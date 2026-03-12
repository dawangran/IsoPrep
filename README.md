# IsoPrep

IsoPrep is a lightweight long-read single-cell preprocessing toolkit.
It orchestrates FASTQ-level preprocessing, barcode/UMI annotation, alignment,
and final BAM/QC materialization for a single sample.

## Features

- Multi-FASTQ single-sample orchestration.
- Modular pipeline stages (`cutadapt`/`tsoclip`/`seqkit`/barcode scan/alignment).
- Unified logging with environment-driven log level and timestamp control.
- QC aggregation from stage summaries and BAM statistics.
- Optional cleanup to keep only final deliverables.

## Project layout

- `IsoPrep/runner.py`: main CLI orchestration logic.
- `IsoPrep/stages.py`: stage-level command wrappers.
- `IsoPrep/config.py`: dataclass-based defaults and executable paths.
- `IsoPrep/utils.py`: shell/filesystem helper utilities.
- `IsoPrep/logging.py`: centralized logger setup.
- `IsoPrep/bin/`: command-line utilities used by the pipeline (barcode scan, sharding helpers, etc.).

## Quick start

```bash
isoprep-runner \
  --fastqs a.fastq.gz b.fastq.gz \
  --sample SAMPLE01 \
  --ref ref.fa \
  --whitelist whitelist.txt \
  --workdir work \
  --threads 16
```

## Key outputs

Under `<workdir>/<sample>/01.data`:

- `<sample>.bam`
- `<sample>.bam.bai`
- `<sample>.qc.tsv`


## CLI parameters

Installed console entry points:

- `isoprep-runner` → `IsoPrep.runner:main`
- `isoprep-run-sharded` → `IsoPrep.bin.run_sharded:main`

### `isoprep-runner`

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--fastqs` | Yes | - | One or more FASTQ files for the same sample. |
| `--sample` | Yes | - | Sample ID used for output naming. |
| `--ref` | Yes | - | Reference FASTA for `minimap2`. |
| `--whitelist` | Yes | - | Cell barcode whitelist used by `scan_cb_umi`. |
| `--valid-list` | No | `""` | Optional valid cell list for `add_cb_umi_db`. |
| `--fulllength` | No | `only` | TSO filter mode: `only` or `all`. |
| `--workdir` | No | `work` | Working directory root. |
| `--threads` | No | `16` | Threads per FASTQ task. |
| `--procs` | No | `1` | Number of FASTQs processed in parallel. |
| `--keep-intermediate` | No | `False` | Keep `tmp/` intermediate files. |
| `--qc-debug` | No | `False` | Print per-FASTQ QC parsing debug logs. |
| `--shards` | No | `32` | Legacy compatibility parameter (unused in no-UMI flow). |
| `--ham` | No | `1` | Legacy compatibility parameter (unused in no-UMI flow). |
| `--ratio` | No | `2.0` | Legacy compatibility parameter (unused in no-UMI flow). |
| `--jitter` | No | `10` | Legacy compatibility parameter (unused in no-UMI flow). |
| `--locus-bin` | No | `1000` | Legacy compatibility parameter (unused in no-UMI flow). |
| `--no-gene` | No | `False` | Legacy compatibility flag (unused in no-UMI flow). |
| `--fast-h1` | No | `False` | Legacy compatibility flag (unused in no-UMI flow). |

### `isoprep-run-sharded`

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--outdir` | Yes | - | Output directory for shard results and merged BAM. |
| `--bam-dir` | Yes | - | Directory containing input BAM files to shard. |
| `--shards` | No | `32` | Number of shard partitions. |
| `--threads` | No | `32` | Threads for `samtools merge/sort`. |
| `--keep-shards` | No | `False` | Keep per-shard directories after final merged output. |

## QC metric definitions

`<sample>.qc.tsv` currently includes:

| Field | Meaning | Source | Parsing rule (priority) | Notes |
|---|---|---|---|---|
| `sample` | Sample identifier | CLI argument | `--sample` | Label only |
| `raw_fastq_reads` | Raw reads before adapter/model filtering | `tmp/<fq>/01.barcode/read1_model.out` | `Total reads processed:\s+([\d,]+)` | Summed across FASTQs |
| `model_reads` | Reads retained by model/TSO filtering | `read1_model.out` and `tso.out` | Parse `Reads written(...)`; if `tso.out` has `trimmed_written=...`, it overrides | `tso.out` has final priority |
| `barcode_corrected_reads` | Reads with successful CB correction and FASTQ emission | `tmp/<fq>/01.barcode/split.tsv.summary.txt` | Prefer `FASTQ_OUTPUT`, fallback `CORR_BOTH_OK`; if summary missing, fallback to log line `FASTQ output:` | Summed across FASTQs |
| `valid_reads` | Reads kept after `add_cb_umi_db` join/filter | `tmp/<fq>/01.barcode/*.summary.txt` | Prefer `model.retain.mask.drc.merge.valid.fq.gz.summary.txt`, otherwise scan `*.summary.txt` for `OUT_KEPT`/`OUT_KEPT(joined)` | Summed across FASTQs |
| `aligned_mapped_reads` | Primary mapped aligned reads | `tmp/<fq>/02.align/*.bam` (or recursive fallback) | `samtools view -c -F 2308` (exclude unmapped/secondary/supplementary) | Summed across FASTQs |

> Note: metrics are stage-specific and not guaranteed to be strictly monotonic across columns.

## Logging

- `SCLR_LOG_LEVEL` (default: `INFO`)
- `SCLR_LOG_TIME` (default: `1`)

Example:

```bash
SCLR_LOG_LEVEL=DEBUG SCLR_LOG_TIME=0 isoprep-runner ...
```
