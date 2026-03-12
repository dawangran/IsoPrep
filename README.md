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


## QC 指标口径表

`<sample>.qc.tsv` 当前输出列如下：

| 字段名 | 指标含义 | 统计来源 | 解析规则（优先级） | 备注 |
|---|---|---|---|---|
| `sample` | 样本名 | CLI 参数 | `--sample` | 仅标识 |
| `raw_fastq_reads` | 原始 read 数（进入 cutadapt 前口径） | `tmp/<fq>/01.barcode/read1_model.out` | 匹配 `Total reads processed:\s+([\d,]+)` | FASTQ 级汇总后按样本求和 |
| `model_reads` | 通过模型/TSO 过滤后的 read 数 | `read1_model.out` 与 `tso.out` | 先匹配 `Reads written(...)`，若 `tso.out` 存在 `trimmed_written=...` 则覆盖 | 以 `tso.out` 为最终优先口径 |
| `barcode_corrected_reads` | 完成 CB 纠错并输出到 FASTQ 的 read 数 | `tmp/<fq>/01.barcode/split.tsv.summary.txt` | 先读 `FASTQ_OUTPUT`，回退 `CORR_BOTH_OK`；若 summary 缺失，回退日志文本 `FASTQ output:` | FASTQ 级汇总后按样本求和 |
| `valid_reads` | add_cb_umi_db 保留（joined）read 数 | `tmp/<fq>/01.barcode/*.summary.txt` | 优先默认文件 `model.retain.mask.drc.merge.valid.fq.gz.summary.txt`，否则遍历 `*.summary.txt`，匹配 `OUT_KEPT` 或 `OUT_KEPT(joined)` | FASTQ 级汇总后按样本求和 |
| `aligned_mapped_reads` | 比对成功的主比对 reads | `tmp/<fq>/02.align/*.bam`（或回退搜索） | `samtools view -c -F 2308`（过滤 unmapped/secondary/supplementary） | 统计 primary mapped，FASTQ 级汇总后按样本求和 |
| `umi_dedup_reads` | UMI 去重后 reads | 固定值 | `NA` | 当前流程不做 UMI dedup |

> 说明：不同指标位于不同阶段，存在“口径不一致”的天然现象（例如 `valid_reads` 与 `aligned_mapped_reads` 不一定单调）。

## Logging

- `SCLR_LOG_LEVEL` (default: `INFO`)
- `SCLR_LOG_TIME` (default: `1`)

Example:

```bash
SCLR_LOG_LEVEL=DEBUG SCLR_LOG_TIME=0 python -m IsoPrep.runner ...
```
