
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pipeline stage wrappers.

Each function assembles and executes one stage command and returns the key
artifact path for the next stage.
"""

from pathlib import Path
from typing import Literal

from sclrtoolkit.logging import setup_logger
from .config import ToolPaths, AdapterModel, Defaults
from .utils import run_cmd, safe_mkdir

logger = setup_logger(__name__)

def cutadapt_stage(
    fq: Path,
    outdir: Path,
    tools: ToolPaths,
    model: AdapterModel,
    df: Defaults,
    threads: int,
    fullength: Literal["only", "all"],
) -> Path:
    """Run adapter/TSO processing and return trimmed FASTQ."""
    safe_mkdir(outdir / "01.barcode")
    fq1 = outdir / "01.barcode" / "model.read1_model.fq.gz"
    fq2 = outdir / "01.barcode" / "model.read1_model_trim.fq.gz"
    hits = outdir / "01.barcode" / "hits.tsv"
    out1 = outdir / "01.barcode" / "read1_model.out"
    out2 = outdir / "01.barcode" / "tso.out"
    cmd1 = (
        f'{tools.cutadapt} -j {threads} -g "{model.model}" -o {fq1} {fq} '
        f'--revcomp --action=retain -m {df.min_read1_len} --rename="{{header}}" --discard-untrimmed > {out1}'
    )
    run_cmd(cmd1)
    if fullength == "only":
        cmd2 = (
            f'{tools.tsoclip} --fastq {fq1} --tso {model.tso} --threads 8 --tail-window 100 --tso-max-mm 4 --tso-max-shift 4 --tso-min-overlap 18 --tso-max-mmr 0.20  --tso-max-hits 10 --min-spacing 6 --n-as-match  --batch-size 40000  --out-tsv {hits} --batch-size 40000 --no-json --gzip-level 1 --gzbuf-kb 1024 --emit-only-hit 1   --out-trim-fastq {fq2} 2> {out2}')
    elif fullength == "all":
        cmd2 = (
            f'{tools.tsoclip} --fastq {fq1} --tso {model.tso} --threads 8 --tail-window 100 --tso-max-mm 4 --tso-max-shift 4 --tso-min-overlap 18 --tso-max-mmr 0.20  --tso-max-hits 10 --min-spacing 6 --n-as-match  --batch-size 40000  --out-tsv {hits} --batch-size 40000 --no-json --gzip-level 1 --gzbuf-kb 1024 --out-trim-fastq {fq2} 2> {out2}')
    else:
        raise ValueError(f"Invalid --fulllength value: {fullength}")

    run_cmd(cmd2)
    return fq2

def seqkit_slice(fq_retain: Path, outdir: Path, tools: ToolPaths, model: AdapterModel, threads: int):
    """Slice barcode-containing windows from trimmed reads."""
    r1 = outdir / "01.barcode" / "model.retain.drc_1.fq.gz"
    r2 = outdir / "01.barcode" / "model.retain.drc_2.fq.gz"
    run_cmd(f'{tools.seqkit} subseq -j {threads} -r {model.r1_slice_1} {fq_retain} -o {r1}')
    run_cmd(f'{tools.seqkit} subseq -j {threads} -r {model.r1_slice_2} {fq_retain} -o {r2}')
    return r1, r2

def scan_cb_umi(r1_slice: Path, outdir: Path, whitelist: Path, threads: int, tools: ToolPaths):
    """Run CB/UMI scanner and return summary TSV + filtered FASTQ."""
    out_tsv = outdir / "01.barcode" / "split.tsv"
    r1_flt = outdir / "01.barcode" / "model.retain.drc_1.filter.fq.gz"
    cmd = (
        f'pigz -dc {r1_slice} | '
        f'{tools.python} -m sclrtoolkit.bin.scan_cb_umi - '
        f'--out_tsv {out_tsv} '
        f'--out_fastq {r1_flt} '
        f'--threads {threads} --mp_chunksize 256 --batch 100000 '
        f'--cb_correct --whitelist {whitelist} --cb_allow_ham2_unique'
    )
    run_cmd(cmd)
    return out_tsv, r1_flt

def add_cb_umi(r1_fq: Path, r2_fq: Path, outdir: Path, sample: str, valid_list: Path, model: AdapterModel, tools: ToolPaths):
    """Attach corrected CB/UMI tags and return merged FASTQ."""
    merged = outdir / "01.barcode" / "model.retain.mask.drc.merge.valid.fq.gz"
    cmd = (
        f'{tools.python} -m sclrtoolkit.bin.add_cb_umi_db '
        f'--r1 {r1_fq} --r2 {r2_fq} '
        f'--out {merged} '
        f'--model {model.cb_umi_model} '
        f'{("--valid_list " + str(valid_list)) if valid_list else ""} '
        f'--sample {sample}'
    )
    run_cmd(cmd)
    return merged

def align_and_tag(fq_merged: Path, ref_fasta: Path, outdir: Path, tools: ToolPaths, threads: int) -> Path:
    """Align merged reads and convert to a tagged BAM."""
    bam_dir = outdir / "02.align"
    safe_mkdir(bam_dir)
    sam = bam_dir / "aln.cDNA.sorted.sam"
    bam = bam_dir / "aln.cDNA.sorted.bam"
    bam_tag = bam_dir / "aln.cDNA.sorted.tag.bam"
    run_cmd(
        f'{tools.minimap2} -t {threads} -ax splice -k14 --secondary=no -C5 -p 0.9 -G 200k --MD -Y --cs=long '
        f'{ref_fasta} {fq_merged} | {tools.samtools} sort -@{max(1, threads//2)} -O BAM -o {bam}'
    )
    run_cmd(f'{tools.samtools} view -h {bam} > {sam}')
    run_cmd(f'{tools.pisa} sam2bam {sam} -o {bam_tag}')
    return bam_tag
