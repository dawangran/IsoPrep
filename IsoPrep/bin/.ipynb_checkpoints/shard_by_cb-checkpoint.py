#!/usr/bin/env python3
import sys
# -*- coding: utf-8 -*-
"""
shard_by_cb.py — Split multiple input BAMs into N shards by CB tag (same CB -> same shard)

Usage:
  {sys.executable} -m sclrtoolkit.bin.shard_by_cb --outdir work_sharded --shards 32 bam/*.bam
"""
import argparse, sys, pathlib, pysam, hashlib

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--shards", type=int, required=True)
    ap.add_argument("--keep-all", action="store_true",
                    help="Keep unmapped/secondary/supplementary (default: skip)")
    ap.add_argument("bams", nargs="+", help="Input BAM(s)")
    return ap.parse_args()

def stable_hash(s: str) -> int:
    # stable 32-bit hash
    return int(hashlib.md5(s.encode()).hexdigest()[:8], 16)

def ensure_dirs(base, n):
    basep = pathlib.Path(base)
    basep.mkdir(parents=True, exist_ok=True)
    paths = []
    for i in range(n):
        p = basep / f"shard_{i}"
        p.mkdir(parents=True, exist_ok=True)
        paths.append(p)
    return paths

def open_shard_writers(template_bam: pysam.AlignmentFile, outdirs, n):
    writers = []
    for i in range(n):
        out_bam = str(outdirs[i] / f"shard_{i}.bam")
        writers.append(pysam.AlignmentFile(out_bam, "wb", template=template_bam))
    return writers

def main():
    args = parse_args()
    outdirs = ensure_dirs(args.outdir, args.shards)
    template = pysam.AlignmentFile(args.bams[0], "rb")
    writers = open_shard_writers(template, outdirs, args.shards)

    n_written = [0]*args.shards
    n_total = 0
    n_skipped = 0

    for bam_path in args.bams:
        bam = pysam.AlignmentFile(bam_path, "rb")
        for r in bam:
            n_total += 1
            if not args.keep_all:
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    n_skipped += 1
                    continue
            if not r.has_tag("CB"):
                n_skipped += 1
                continue
            cb = r.get_tag("CB")
            shard = stable_hash(cb) % args.shards
            writers[shard].write(r)
            n_written[shard] += 1
        bam.close()

    for w in writers: w.close()
    template.close()

    sys.stderr.write(f"[shard_by_cb] total={n_total} skipped={n_skipped} "
                     f"kept={sum(n_written)} shards={args.shards} "
                     f"per_shard={n_written}\n")

if __name__ == "__main__":
    main()
