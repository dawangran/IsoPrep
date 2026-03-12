# IsoPrep

IsoPrep is a lightweight long-read single-cell preprocessing pipeline package.
It orchestrates FASTQ-level preprocessing, barcode/UMI annotation, alignment,
and final BAM/QC materialization for one sample.

## Features

- Multi-FASTQ single-sample orchestration.
- Modular pipeline stages (`cutadapt`/`tsoclip`/`seqkit`/barcode scan/alignment).
- Unified logging with environment-driven log level and timestamp control.
- Basic QC aggregation from stage summaries and BAM statistics.
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

## Logging

- `SCLR_LOG_LEVEL` (default: `INFO`)
- `SCLR_LOG_TIME` (default: `1`)

Example:

```bash
SCLR_LOG_LEVEL=DEBUG SCLR_LOG_TIME=0 python -m IsoPrep.runner ...
```
