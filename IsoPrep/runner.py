#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Single-sample multi-FASTQ orchestrator.

Pipeline behavior is intentionally preserved, while QC statistics are parsed
from stage outputs and final artifacts are written to ``01.data``.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from multiprocessing import Pool
from typing import Optional, Tuple

from .logging import setup_logger
from .utils import safe_mkdir, symlink_force, run_cmd
from .config import ToolPaths, AdapterModel, Defaults, ShardedParams
from . import stages

logger = setup_logger(__name__)

# -------------------- helpers: io / names --------------------
def _read_text(path: Path) -> str:
    if not path or not path.exists():
        return ""
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""

def fastq_basename(fq_path: Path) -> str:
    """Robust FASTQ basename: *.fastq.gz / *.fq.gz / *.fastq / *.fq -> stem without extensions."""
    name = fq_path.name
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suf):
            return name[: -len(suf)]
    # fallback: strip .gz if present, then last extension
    base = os.path.splitext(name)[0]
    base = os.path.splitext(base)[0]
    return base

# -------------------- helpers: parse summaries / counters --------------------
def parse_cutadapt_summary(bc_dir: Path) -> Tuple[int, int]:
    """Parse raw/full-length read counts from cutadapt and tsoclip outputs.

    Priority:
    - `01.barcode/read1_model.out`: `Total reads processed` for raw count.
    - `01.barcode/read1_model.out`: `Reads written (...)` as first full-length fallback.
    - `01.barcode/tso.out`: `trimmed_written=<N>` as preferred retained-read source.
    """
    raw = full_len = 0

    read1_out = bc_dir / "read1_model.out"
    txt = _read_text(read1_out)
    if txt:
        m_raw = re.search(r"Total reads processed:\s+([\d,]+)", txt)
        if m_raw:
            raw = int(m_raw.group(1).replace(",", ""))

        # cutadapt summary fallback, e.g. "Reads written (passing filters): 123"
        m_written = re.search(r"Reads written(?: \([^)]*\))?:\s+([\d,]+)", txt)
        if m_written:
            full_len = int(m_written.group(1).replace(",", ""))

    retain_out = bc_dir / "tso.out"
    t2 = _read_text(retain_out)
    if t2:
        m_done = re.search(r"trimmed_written=([\d,]+)", t2)
        if m_done:
            full_len = int(m_done.group(1).replace(",", ""))

    return raw, full_len
def parse_scan_summary(bc_dir: Path) -> int:
    """Parse corrected barcode read count from scan summary outputs."""
    p = bc_dir / "split.tsv.summary.txt"
    txt = _read_text(p)
    if txt:
        m = re.search(r"^FASTQ_OUTPUT\s+(\d+)", txt, flags=re.M)
        if not m:
            m = re.search(r"^CORR_BOTH_OK\s+(\d+)", txt, flags=re.M)
        if m:
            return int(m.group(1))
    else:
        # Fallback: some runs only preserve stdout logs.
        for f in (bc_dir / "split.tsv", bc_dir / "scan_cb_umi.log"):
            t = _read_text(f)
            if not t:
                continue
            m = re.search(r"FASTQ output:\s*([\d,]+)", t)
            if m:
                return int(m.group(1).replace(",", ""))
    return 0


def parse_addcb_summary(bc_dir: Path) -> int:
    """Parse kept-read count from add_cb_umi summary files."""
    # Match lines like: OUT_KEPT(joined)    399
    pat = re.compile(r"^OUT_KEPT(?:\(joined\))?\s+(\d+)", flags=re.M)

    # Try default naming first.
    candidates = [bc_dir / "model.retain.mask.drc.merge.valid.fq.gz.summary.txt"]
    # Then include wildcard candidates.
    candidates += [p for p in sorted(bc_dir.glob("*.summary.txt")) if p not in candidates]

    for p in candidates:
        txt = _read_text(p)
        if not txt:
            continue
        m = pat.search(txt)
        if m:
            return int(m.group(1))
    return 0


def count_bam_reads(bam: Path, tools: ToolPaths, unique_only: bool = True) -> int:
    """Count mapped reads with samtools.
    If unique_only=True, count only primary alignments (exclude secondary/supplementary/unmapped)."""
    if not bam or not bam.exists():
        logger.warning(f"[QC] BAM not found for counting: {bam}")
        return 0
    if unique_only:
        # 2308 = unmapped(4) + secondary(256) + supplementary(2048)
        flag = "-F 2308"
    else:
        # only drop unmapped
        flag = "-F 4"
    cmd = f'{tools.samtools} view -c {flag} {bam}'.strip()
    try:
        out = subprocess.check_output(cmd, shell=True, text=True)
        return int(out.strip())
    except Exception as e:
        logger.warning(f"[QC] samtools count failed for {bam}: {e}")
        return 0


# -------------------- per-FASTQ processing (unchanged) --------------------
def process_one_fastq(
    fq: Path, sample_dir: Path, sample: str, ref_fasta: Path, whitelist: Path,
    valid_list: Path, tools: ToolPaths, model: AdapterModel,
    df: Defaults, threads: int, final_dir: Path, fullength: bool
) -> Path:
    """
    Process one FASTQ to a tagged BAM under tmp/, and symlink it into final_dir.
    Return the symlinked BAM path for downstream sample-level merge.
    """
    fq = fq.resolve()
    fq_basename = fastq_basename(fq)
    tmpdir = sample_dir / "tmp" / fq_basename
    safe_mkdir(tmpdir)

    logger.info(f"[{sample}/{fq.name}] 1) cutadapt")
    retain = stages.cutadapt_stage(fq, tmpdir, tools, model, df, threads,fullength)

    logger.info(f"[{sample}/{fq.name}] 2) seqkit slice")
    r1, r2 = stages.seqkit_slice(retain, tmpdir, tools, model, threads)

    logger.info(f"[{sample}/{fq.name}] 3) scan_cb_umi")
    out_tsv, r1_flt = stages.scan_cb_umi(r1, tmpdir, whitelist, threads, tools)

    logger.info(f"[{sample}/{fq.name}] 4) add_cb_umi_db (accelerated)")
    merged_fq = stages.add_cb_umi(r1_flt, r2, tmpdir, sample, valid_list, model, tools)

    logger.info(f"[{sample}/{fq.name}] 5) minimap2 + PISA")
    bam_tag = stages.align_and_tag(merged_fq, ref_fasta, tmpdir, tools, threads)

    # link per-FASTQ BAM into final_dir
    target = final_dir / f"{sample}.{fq_basename}.tag.bam"
    symlink_force(bam_tag, target)
    run_cmd(f"{tools.samtools} index -@ {max(1, threads//2)} {target}")
    return target

# -------------------- post-hoc QC aggregator (file parsing, robust paths) --------------------
def aggregate_qc_from_tmp(sample_dir: Path, tools: ToolPaths, debug: bool = False) -> Tuple[int, int, int, int, int]:
    """
    扫描 sample_dir/tmp 下每个 FASTQ 子目录，累加：
      raw（cutadapt #1；read1_model.out）
      full_len（cutadapt #1；read1_model.out；回退 tso.out）
      bc_corrected（split.tsv.summary.txt：FASTQ_OUTPUT 优先，回退 CORR_BOTH_OK）
      valid（*.summary.txt：KEPT）
      aligned_mapped（samtools view -c -F 4 对齐 BAM）
    兼容目录名：01.barcode 或 barcode；02.align 或 align；必要时回退 rglob。
    同时打印每个 FASTQ 子目录的解析结果；当 debug=True 时打印所用文件路径。
    """
    tmp_root = sample_dir / "tmp"
    if not tmp_root.exists():
        logger.warning(f"[QC] tmp directory not found: {tmp_root}")
        return 0, 0, 0, 0, 0

    total_raw = total_full = total_corr = total_valid = total_aln = 0

    for fq_dir in sorted([d for d in tmp_root.iterdir() if d.is_dir()]):
        # --- barcode dir ---
        bc_candidates = [fq_dir / "01.barcode", fq_dir / "barcode"]
        bc_dir = next((p for p in bc_candidates if p.exists()), None)
        if bc_dir is None:
            got = list(fq_dir.rglob("read1_model.out"))
            bc_dir = got[0].parent if got else None

        # --- align dir ---
        aln_candidates = [fq_dir / "02.align", fq_dir / "align"]
        aln_dir = next((p for p in aln_candidates if p.exists()), None)

        raw = full_len = corr = valid = aln = 0

        if bc_dir is None:
            logger.warning(f"[QC] barcode dir not found under: {fq_dir}")
        else:
            # cutadapt
            r, f = parse_cutadapt_summary(bc_dir)
            if debug:
                logger.info(f"[QC-debug] {fq_dir.name} cutadapt: {(bc_dir/'read1_model.out')} -> raw={r}, full_len={f}")
            raw += r; full_len += f

            # scan_cb_umi
            c = parse_scan_summary(bc_dir)
            if debug:
                logger.info(f"[QC-debug] {fq_dir.name} scan: {(bc_dir/'split.tsv.summary.txt')} -> bc_corrected={c}")
            corr += c

            # add_cb_umi
            v = parse_addcb_summary(bc_dir)
            if debug:
                logger.info(f"[QC-debug] {fq_dir.name} add_cb: search in {bc_dir} -> valid={v}")
            valid += v

        # 对齐 BAM
        aln_bam = None
        if aln_dir:
            cand = aln_dir / "aln.cDNA.bam"
            if cand.exists():
                aln_bam = cand
            else:
                cands = sorted(aln_dir.glob("*.bam"))
                if cands:
                    aln_bam = cands[0]
        else:
            cands = sorted(fq_dir.rglob("*.bam"))
            if cands:
                aln_bam = cands[0]

        if aln_bam:
            a = count_bam_reads(aln_bam, tools, unique_only=True)
            if debug:
                logger.info(f"[QC-debug] {fq_dir.name} align: {aln_bam} -> aligned_mapped={a}")
            aln += a
        else:
            logger.warning(f"[QC] aligned BAM not found under: {aln_dir or fq_dir}")

        # 汇总日志：每个 FASTQ 子目录的 QC
        logger.info(
            f"[QC] {fq_dir.name}: raw={raw}, full_len={full_len}, "
            f"bc_corrected={corr}, valid={valid}, aligned_mapped={aln}"
        )

        total_raw   += raw
        total_full  += full_len
        total_corr  += corr
        total_valid += valid
        total_aln   += aln

    return total_raw, total_full, total_corr, total_valid, total_aln

# -------------------- main --------------------
def main():
    ap = argparse.ArgumentParser(
        description="IsoPrep: single-sample runner (multi-FASTQ -> per-FASTQ BAMs -> merged BAM; QC via file parsing; finals in 01.data)"
    )
    ap.add_argument("--fastqs", nargs="+", required=True, help="List of FASTQ files (all belong to the same sample)")
    ap.add_argument("--sample", required=True, help="Sample ID for this run")
    ap.add_argument("--ref", required=True, help="Reference FASTA for minimap2")
    ap.add_argument("--whitelist", required=True, help="CB whitelist for scan_cb_umi")
    ap.add_argument("--valid-list", default="", help="valid.cell.tsv (optional)")
    ap.add_argument("--fulllength", choices=["only", "all"], default="only", help="TSO filtering mode")
    ap.add_argument("--workdir", default="work", help="Working directory")
    ap.add_argument("--threads", type=int, default=16, help="Threads per FASTQ")
    ap.add_argument("--procs", type=int, default=1, help="Parallel FASTQs")
    ap.add_argument("--keep-intermediate", action="store_true",
                    help="Keep tmp/ intermediates. If not set, only final outputs in 01.data are kept.")
    ap.add_argument("--qc-debug", action="store_true",
                    help="Print per-FASTQ QC file paths and parsed numbers")
    # Legacy sharding parameters (currently unused; retained for CLI compatibility).
    ap.add_argument("--shards", type=int, default=ShardedParams.shards)
    ap.add_argument("--ham", type=int, default=ShardedParams.ham)
    ap.add_argument("--ratio", type=float, default=ShardedParams.ratio)
    ap.add_argument("--jitter", type=int, default=ShardedParams.jitter)
    ap.add_argument("--locus-bin", type=int, default=ShardedParams.locus_bin)
    ap.add_argument("--no-gene", action="store_true", default=False)
    ap.add_argument("--fast-h1", action="store_true", default=False)
    args = ap.parse_args()

    fastqs = [Path(f).resolve() for f in args.fastqs]
    workdir = Path(args.workdir).resolve()
    safe_mkdir(workdir)
    sample_dir = workdir / args.sample
    safe_mkdir(sample_dir)
    final_dir = sample_dir / "01.data"
    safe_mkdir(final_dir)

    tools = ToolPaths(); model = AdapterModel(); df = Defaults()
    ref = Path(args.ref).resolve()
    wl = Path(args.whitelist).resolve()
    valid_list = Path(args.valid_list).resolve() if args.valid_list else None

    # 1) Per-FASTQ → final_dir symlinked BAMs (UNCHANGED)
    tasks = [(fq, sample_dir, args.sample, ref, wl, valid_list, tools, model, df, args.threads, final_dir, args.fulllength) for fq in fastqs]
    if args.procs and args.procs > 1:
        with Pool(processes=args.procs) as pool:
            for _ in pool.starmap(process_one_fastq, tasks):
                pass
    else:
        for t in tasks:
            process_one_fastq(*t)

    
    # 2) UMI correction is intentionally skipped: merge per-FASTQ tag BAMs, then sort/index.
    # Collect per-FASTQ BAMs from final_dir
    tag_bams = sorted([p for p in final_dir.glob(f"{args.sample}.*.tag.bam") if p.is_file()])
    if not tag_bams:
        raise SystemExit("No per-FASTQ tag BAMs found under final_dir")
    merged_bam = final_dir / f"{args.sample}.merged.bam"
    run_cmd(f"{tools.samtools} merge -@ {args.threads} -l 1 -o {merged_bam} " + " ".join(map(str, tag_bams)))
    final_dedup_sorted = final_dir / f"{args.sample}.bam"
    run_cmd(f"{tools.samtools} sort -@ {args.threads} -o {final_dedup_sorted} {merged_bam}")
    run_cmd(f"{tools.samtools} index -@ {max(1, args.threads//2)} {final_dedup_sorted}")

    # 3) Aggregate QC metrics across all FASTQs.
    raw, full_len, bc_corr, valid, aligned = aggregate_qc_from_tmp(sample_dir, tools, debug=args.qc_debug)

    bai = final_dir / f"{args.sample}.bam.bai"
    if final_dedup_sorted.exists() and not bai.exists():
        run_cmd(f"{tools.samtools} index -@ {max(1, args.threads//2)} {final_dedup_sorted}")

    qc_path = final_dir / f"{args.sample}.qc.tsv"
    with open(qc_path, "w", encoding="utf-8") as qf:
        qf.write("sample	raw_fastq_reads	model_reads	barcode_corrected_reads	valid_reads	aligned_mapped_reads\n")
        qf.write(f"{args.sample}	{raw}	{full_len}	{bc_corr}	{valid}	{aligned}")
    logger.info(f"QC summary written: {qc_path}")

    # 4) Cleanup (keep only final outputs unless --keep-intermediate).
    if not args.keep_intermediate:
        # Remove tmp/.
        tmp_root = sample_dir / "tmp"
        if tmp_root.exists():
            shutil.rmtree(tmp_root, ignore_errors=True)
        # Keep only final .bam/.bai and QC summary.
        keep = {
            f"{args.sample}.bam",
            f"{args.sample}.bam.bai",
            f"{args.sample}.qc.tsv",
        }
        for p in final_dir.iterdir():
            if p.name not in keep:
                try:
                    shutil.rmtree(p, ignore_errors=True) if p.is_dir() else p.unlink()
                except Exception:
                    pass

    logger.info(f"Final outputs: {final_dir / (args.sample + '.bam')} (.bai); QC: {qc_path}")

if __name__ == "__main__":
    main()
