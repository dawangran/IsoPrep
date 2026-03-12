#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shard BAMs by cell barcode (CB) and optionally merge shard outputs.

This utility keeps the sharding workflow but removes the deprecated
`umi_correct_longreads` reduction stage.
"""

from __future__ import annotations

import argparse
import glob
import os
import pathlib
import shutil
import subprocess
import sys


def run(cmd: str) -> None:
    """Run a shell command and fail fast on non-zero exit code."""
    sys.stderr.write(f"[cmd] {cmd}\n")
    ret = subprocess.call(cmd, shell=True)
    if ret != 0:
        raise RuntimeError(f"Command failed ({ret}): {cmd}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Shard BAM files by CB tag and merge shard BAMs (no UMI correction)."
    )
    ap.add_argument("--outdir", required=True, help="Output directory for shards and merged BAM")
    ap.add_argument("--shards", type=int, default=32, help="Number of shard partitions")
    ap.add_argument("--bam-dir", required=True, help="Directory containing input BAM(s)")
    ap.add_argument("--threads", type=int, default=32, help="Threads for samtools merge/sort")
    ap.add_argument(
        "--keep-shards",
        action="store_true",
        help="Keep shard BAM files after merged output is generated",
    )
    args = ap.parse_args()

    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)

    bam_paths = sorted(glob.glob(os.path.join(args.bam_dir, "*.bam")))
    if not bam_paths:
        raise SystemExit("No BAM files found in --bam-dir")

    # 1) Shard by CB.
    shard_cmd = (
        f"{sys.executable} -m IsoPrep.bin.shard_by_cb "
        f"--outdir {args.outdir} --shards {args.shards} " + " ".join(bam_paths)
    )
    run(shard_cmd)

    # 2) Merge shard BAMs to a deterministic final output.
    shard_bams = sorted(glob.glob(f"{args.outdir}/shard_*/shard_*.bam"))
    if not shard_bams:
        raise SystemExit("No shard BAMs generated under outdir/shard_*/")

    merged_bam = f"{args.outdir}/sample.sharded.bam"
    sorted_bam = f"{args.outdir}/sample.sharded.sorted.bam"
    run(f"samtools merge -@ {args.threads} -l 1 -o {merged_bam} " + " ".join(shard_bams))
    run(f"samtools sort -@ {args.threads} -o {sorted_bam} {merged_bam}")
    run(f"samtools index {sorted_bam}")

    if not args.keep_shards:
        for shard_dir in glob.glob(f"{args.outdir}/shard_*"):
            if os.path.isdir(shard_dir):
                shutil.rmtree(shard_dir, ignore_errors=True)

    sys.stderr.write("[run_sharded] ALL DONE\n")
    sys.stderr.write(f"  BAM : {sorted_bam} (.bai)\n")


if __name__ == "__main__":
    main()
