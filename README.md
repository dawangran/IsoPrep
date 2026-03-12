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
python -m IsoPrep.runner \
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

### `python -m IsoPrep.runner`

| Parameter | Required | Description |
|---|---|---|
| `--fastqs` | Yes | One or more FASTQ files for the same sample. |
| `--sample` | Yes | Sample ID used for output naming. |
| `--ref` | Yes | Reference FASTA for `minimap2`. |
| `--whitelist` | Yes | Cell barcode whitelist used by `scan_cb_umi`. |
| `--valid-list` | No | Optional valid cell list for `add_cb_umi_db`. |
| `--fulllength` | No | TSO filter mode: `only` or `all` (default: `only`). |
| `--workdir` | No | Working directory root (default: `work`). |
| `--threads` | No | Threads per FASTQ task (default: `16`). |
| `--procs` | No | Number of FASTQs processed in parallel (default: `1`). |
| `--keep-intermediate` | No | Keep `tmp/` intermediate files. |
| `--qc-debug` | No | Print per-FASTQ QC parsing debug logs. |
| `--shards`, `--ham`, `--ratio`, `--jitter`, `--locus-bin`, `--no-gene`, `--fast-h1` | No | Legacy compatibility parameters; currently unused in the no-UMI flow. |

### `python -m IsoPrep.bin.run_sharded`

| Parameter | Required | Description |
|---|---|---|
| `--outdir` | Yes | Output directory for shard results and merged BAM. |
| `--bam-dir` | Yes | Directory containing input BAM files to shard. |
| `--shards` | No | Number of shard partitions (default: `32`). |
| `--threads` | No | Threads for `samtools merge/sort` (default: `32`). |
| `--keep-shards` | No | Keep per-shard directories after final merged output. |

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
SCLR_LOG_LEVEL=DEBUG SCLR_LOG_TIME=0 python -m IsoPrep.runner ...
```
