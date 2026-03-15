<h1 align="center">IsoPrep</h1>

<p align="center">
  <strong>A lightweight and production-friendly long-read single-cell preprocessing toolkit.</strong>
</p>

<p align="center">
  FASTQ preprocessing • CB/UMI annotation • alignment • BAM + QC delivery
</p>

---

## Why IsoPrep

IsoPrep is designed for **single-sample long-read scRNA/scISO workflows** where you need a clear, reproducible, and scriptable preprocessing path from raw FASTQ inputs to final BAM/QC artifacts.

It provides:

- **End-to-end orchestration** for multi-FASTQ single-sample processing.
- **Modular stage execution** (`cutadapt`, `tsoclip`, `seqkit`, barcode scan, alignment).
- **Consistent logging** with environment-controlled verbosity and timestamps.
- **QC aggregation** across stage summaries and alignment statistics.
- **Optional cleanup mode** for compact final deliverables.

## Installation

### From source (recommended)

```bash
git clone <your-repo-url>
cd IsoPrep
pip install .
```

### Requirements

- Python `>=3.8`
- Runtime Python dependency:
  - `pysam`
- External tools used by pipeline stages should be available in `PATH` (for example: `cutadapt`, `seqkit`, `samtools`, `minimap2`, and stage-related helper tools).

## Quick start

```bash
isoprep \
  --fastqs a.fastq.gz b.fastq.gz \
  --sample SAMPLE01 \
  --ref ref.fa \
  --whitelist whitelist.txt \
  --workdir work \
  --threads 16
```

## Output structure

Primary outputs are generated under:

```text
<workdir>/<sample>/01.data/
```

Key files:

- `<sample>.bam`
- `<sample>.bam.bai`
- `<sample>.qc.tsv`

## Command-line interface

### Installed entry point

- `isoprep` → `IsoPrep.runner:main`

### `isoprep` parameters

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

> Note: these metrics are stage-specific and are not guaranteed to be strictly monotonic across columns.

## Logging

Environment variables:

- `SCLR_LOG_LEVEL` (default: `INFO`)
- `SCLR_LOG_TIME` (default: `1`)

Example:

```bash
SCLR_LOG_LEVEL=DEBUG SCLR_LOG_TIME=0 isoprep ...
```

## Project layout

- `IsoPrep/runner.py`: CLI orchestration logic.
- `IsoPrep/stages.py`: stage command wrappers.
- `IsoPrep/config.py`: dataclass-based defaults and executable paths.
- `IsoPrep/utils.py`: shell/filesystem helpers.
- `IsoPrep/logging.py`: centralized logger setup.
- `IsoPrep/bin/`: utility scripts used by the pipeline.

## Naming conventions

- Distribution/package name for installation: `isoprep` (lowercase).
- Python source package: `IsoPrep`.
- Main CLI entry point: `isoprep`.

## License

Released under the terms of the `LICENSE` file in this repository.
