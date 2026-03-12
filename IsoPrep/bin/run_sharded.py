#!/usr/bin/env python3
import sys
# -*- coding: utf-8 -*-
"""
run_sharded.py — Orchestrate sharding by CB and per-shard UMI correction

Usage:
  python run_sharded.py \
    --outdir work_sharded \
    --shards 32 \
    --bam-dir bam_dir \
    --threads 32 \
    --no-gene \
    --fast-h1
"""
import argparse, os, sys, glob, subprocess, pathlib, multiprocessing as mp

def run(cmd):
    sys.stderr.write(f"[cmd] {cmd}\n")
    ret = subprocess.call(cmd, shell=True)
    if ret != 0:
        raise RuntimeError(f"Command failed: {cmd}")

def per_shard_task(args):
    shard_id, outdir, extra = args
    shard_bam = f"{outdir}/shard_{shard_id}/shard_{shard_id}.bam"
    out_prefix = f"{outdir}/shard_{shard_id}/result"
    cmd = f"{sys.executable} -m sclrtoolkit.bin.umi_correct_longreads --bam {shard_bam} --out {out_prefix} {extra}"
    run(cmd)
    return shard_id

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--shards", type=int, default=32)
    ap.add_argument("--bam-dir", dest="bam_dir", required=True, help="Directory containing input BAM(s)")
    ap.add_argument("--threads", type=int, default=32)
# [NO-UMI]     ap.add_argument("--no-gene", action="store_true", help="propagate to umi_correct_longreads.py")
# [NO-UMI]     ap.add_argument("--fast-h1", action="store_true", help="propagate to umi_correct_longreads.py")
    ap.add_argument("--gtf", help="Optional GTF (used only when --no-gene is OFF)")
    ap.add_argument("--ham", type=int, default=1)
    ap.add_argument("--ratio", type=float, default=2.0)
    ap.add_argument("--jitter", type=int, default=10)
    ap.add_argument("--locus-bin", type=int, default=1000)
    ap.add_argument("--write-all", action="store_true")
    args = ap.parse_args()

    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    bam_paths = sorted(glob.glob(os.path.join(args.bam_dir, "*.bam")))
    if not bam_paths:
        raise SystemExit("No BAM files found in --bam-dir")

    # 1) Shard by CB
    shard_cmd = f"{sys.executable} -m sclrtoolkit.bin.shard_by_cb --outdir {args.outdir} --shards {args.shards} " + " ".join(bam_paths)
    run(shard_cmd)

    # build extra flags for reducer
    extra = f"--ham {args.ham} --ratio {args.ratio} --jitter {args.jitter} --locus-bin {args.locus_bin} "
    if args.no_gene:  extra += "--no-gene "
    if args.fast_h1:  extra += "--fast-h1 "
    if args.write_all: extra += "--write-all "
    if (not args.no_gene) and args.gtf:
        extra += f"--gtf {args.gtf} "

    # 2) Per-shard reduce (parallel)
    tasks = [(i, args.outdir, extra) for i in range(args.shards)]
    with mp.Pool(processes=min(args.threads, args.shards)) as pool:
        for sid in pool.imap_unordered(per_shard_task, tasks):
            sys.stderr.write(f"[reduce] shard {sid} done\n")

    # 3) Merge shard BAMs
    shard_dedups = sorted(glob.glob(f"{args.outdir}/shard_*/result.dedup.bam"))
    if not shard_dedups:
        raise SystemExit("No shard dedup BAMs found")
    merged_bam = f"{args.outdir}/sample.dedup.bam"
    run(f"samtools merge -@ {args.threads} -l 1 -o {merged_bam} " + " ".join(shard_dedups))
    run(f"samtools sort -@ {args.threads} -o {args.outdir}/sample.dedup.sorted.bam {merged_bam}")
    run(f"samtools index {args.outdir}/sample.dedup.sorted.bam")

    # 4) Merge TSVs
    # molecules
    first = True
    with open(f"{args.outdir}/sample.molecules.tsv","w") as out:
        for p in sorted(glob.glob(f"{args.outdir}/shard_*/result.molecules.tsv")):
            with open(p) as fh:
                for i,line in enumerate(fh):
                    if first or i>0:
                        out.write(line)
            first = False
    # assignments
    first = True
    with open(f"{args.outdir}/sample.assignments.tsv","w") as out:
        for p in sorted(glob.glob(f"{args.outdir}/shard_*/result.assignments.tsv")):
            with open(p) as fh:
                for i,line in enumerate(fh):
                    if first or i>0:
                        out.write(line)
            first = False

    sys.stderr.write("[run_sharded] ALL DONE\n")
    sys.stderr.write(f"  BAM : {args.outdir}/sample.dedup.sorted.bam (.bai)\n")
    sys.stderr.write(f"  TSV : {args.outdir}/sample.molecules.tsv\n")
    sys.stderr.write(f"  TSV : {args.outdir}/sample.assignments.tsv\n")

if __name__ == "__main__":
    main()
