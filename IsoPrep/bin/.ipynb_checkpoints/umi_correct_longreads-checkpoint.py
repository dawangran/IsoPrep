#!/usr/bin/env python3
import sys
# -*- coding: utf-8 -*-
"""
umi_correct_longreads.py — Long-read UMI correction & dedup
- Input  : BAM with CB (cell), UR (raw UMI); optional GX (gene_id)
- Output : {out}.dedup.bam  (UB=corrected UMI, DA=1/0)
           {out}.molecules.tsv
           {out}.assignments.tsv

Key design:
  --no-gene OFF: Key=(CB, GX, strand, SJ)
  --no-gene ON :
     if splice (N>=1): Key=(CB, strand, SJ)
     if no-sj (N==0) : Key=(CB, strand, f"no_sj|{chrom:start'-end'}")  # 1kb bin

Usage:
  {sys.executable} -m sclrtoolkit.bin.umi_correct_longreads \
    --bam shard_0.bam --out shard_0/result \
    --ham 1 --ratio 2.0 --jitter 10 --no-gene --fast-h1
"""
import argparse, sys, gzip, collections
import pysam
from collections import defaultdict

# ---------- helpers ----------
def hamming(a, b):
    if len(a)!=len(b): return 999
    return sum(x!=y for x,y in zip(a,b))

def directional_adjacency(umis_cnt, ham_max=1, ratio=2.0):
    """Generic O(K^2) adjacency with directional merge"""
    order = sorted(umis_cnt.items(), key=lambda x:(-x[1], x[0]))
    rep = {u:u for u,_ in order}
    def root(u):
        while rep[u]!=u:
            rep[u]=rep[rep[u]]
            u=rep[u]
        return u
    for i,(ui,ci) in enumerate(order):
        ri = root(ui)
        for uj,cj in order[i+1:]:
            if hamming(ui,uj) <= ham_max:
                rj = root(uj)
                if ri == rj: continue
                if ci >= ratio * cj:
                    rep[rj] = ri
                elif cj >= ratio * ci:
                    rep[ri] = rj
                    ri = rj
    for u in list(rep.keys()):
        rep[u]=root(u)
    return rep

def directional_adjacency_fast_h1(umis_cnt, ratio=2.0):
    """Near O(K·L) adjacency for ham=1 via neighbor enumeration"""
    bases = "ACGT"
    order = sorted(umis_cnt.items(), key=lambda x:(-x[1], x[0]))
    present = set(umis_cnt.keys())
    rep = {u:u for u,_ in order}
    def root(u):
        while rep[u]!=u:
            rep[u]=rep[rep[u]]
            u=rep[u]
        return u
    for ui, ci in order:
        ri = root(ui)
        s = list(ui)
        for i, ch in enumerate(s):
            for b in bases:
                if b == ch: continue
                s[i] = b
                v = "".join(s)
                if v in present:
                    rj = root(v)
                    if ri == rj: 
                        continue
                    cj = umis_cnt[v]
                    if ci >= ratio * cj:
                        rep[rj] = ri
                    elif cj >= ratio * ci:
                        rep[ri] = rj
                        ri = rj
                s[i] = ch
    for u in list(rep.keys()):
        rep[u] = root(u)
    return rep

def cigar_to_blocks_and_introns(pos0, cigar):
    """Return exon blocks and introns as 0-based half-open intervals"""
    ref = pos0
    blocks = []
    block_start = ref
    introns = []
    for op, ln in cigar or []:
        # 0M 1I 2D 3N 4S 5H 6P 7= 8X
        if op in (0,7,8):  # M,=,X
            ref += ln
        elif op == 2:      # D
            ref += ln
        elif op == 3:      # N (intron)
            blocks.append((block_start, ref))
            s, e = ref, ref + ln
            introns.append((s, e))
            ref = e
            block_start = ref
        elif op in (1,4,5,6):  # I/S/H/P: no ref advance
            pass
    blocks.append((block_start, ref))
    return blocks, introns

def make_sj_signature(chrom, strand, introns, jitter=10):
    """N==0 -> '<strand>:no_sj'; N>=1 -> '<strand>:chr:s1'-e1'|...' using floor binning"""
    prefix = "-" if strand=="-" else "+"
    if not introns:
        return f"{prefix}:no_sj"
    parts = []
    for s,e in introns:
        qs = (s // jitter) * jitter
        qe = (e // jitter) * jitter
        parts.append(f"{chrom}:{qs}-{qe}")
    return f"{prefix}:" + "|".join(parts)

def load_gtf_exons(gtf):
    """chrom -> list of (start,end,gene_id,strand)"""
    idx = defaultdict(list)
    if not gtf: return idx
    fh = gzip.open(gtf,'rt') if gtf.endswith(('.gz','.gzip')) else open(gtf,'r')
    with fh:
        for line in fh:
            if not line or line[0]=='#': continue
            toks = line.rstrip('\n').split('\t')
            if len(toks) < 9: continue
            chrom, _src, ftype, start, end, _score, strand, _phase, attr = toks
            if ftype != 'exon': continue
            d={}
            for kv in attr.rstrip(';').split(';'):
                kv=kv.strip()
                if not kv: continue
                k,v = kv.split(' ',1)
                d[k]=v.strip('"')
            gid = d.get('gene_id')
            if not gid: continue
            idx[chrom].append((int(start)-1, int(end), gid, strand))
    for chrom in idx:
        idx[chrom].sort()
    return idx

def infer_gene_by_overlap(chrom, blocks, exon_idx):
    if chrom not in exon_idx: return None
    best_gid, best_ol = None, 0
    for (gs,ge,gid,_st) in exon_idx[chrom]:
        for (s,e) in blocks:
            ol = max(0, min(e, ge) - max(s, gs))
            if ol > best_ol:
                best_ol = ol
                best_gid = gid
    return best_gid

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True, help="input BAM (coord-sorted)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--gtf", help="optional GTF; used only when --no-gene is OFF and GX is absent")
    ap.add_argument("--no-gene", action="store_true",
                    help="ignore gene info; use (CB,strand,SJ) or no-sj locus bin keys")
    ap.add_argument("--ham", type=int, default=1, help="UMI Hamming threshold")
    ap.add_argument("--ratio", type=float, default=2.0, help="merge ratio threshold (big>=ratio*small)")
    ap.add_argument("--jitter", type=int, default=10, help="SJ bin size (bp)")
    ap.add_argument("--locus-bin", type=int, default=1000, help="bin size for no-sj locus (bp) when --no-gene")
    ap.add_argument("--fast-h1", action="store_true",
                    help="use fast near O(K·L) adjacency when ham==1")
    ap.add_argument("--write-all", action="store_true",
                    help="write all reads with DA=1/0 (default: only representatives)")
    args = ap.parse_args()

    exon_idx = {} if args.no_gene else load_gtf_exons(args.gtf)

    bam = pysam.AlignmentFile(args.bam, "rb")
    out_bam = pysam.AlignmentFile(f"{args.out}.dedup.bam", "wb", template=bam)

    buckets = defaultdict(list)   # key -> [(name, UR, MAPQ, qlen)]
    read2UB  = {}                 # read_name -> corrected UMI (UB)
    keep_names = set()

    # 1) Build buckets
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if not (r.has_tag("CB") and r.has_tag("UR")):
            continue
        cb = r.get_tag("CB")
        ur = r.get_tag("UR")
        strand = "-" if r.is_reverse else "+"
        chrom = bam.get_reference_name(r.reference_id)
        blocks, introns = cigar_to_blocks_and_introns(r.reference_start, r.cigartuples)
        sj = make_sj_signature(chrom, strand, introns, jitter=args.jitter)

        if args.no_gene:
            if sj.endswith("no_sj"):
                bstart = min(s for s, e in blocks)
                bend   = max(e for s, e in blocks)
                q = args.locus_bin
                bstart_q = (bstart // q) * q
                bend_q   = (bend   // q) * q
                locus = f"{chrom}:{bstart_q}-{bend_q}"
                key = (cb, strand, f"no_sj|{locus}")
            else:
                key = (cb, strand, sj)
        else:
            gx = r.get_tag("GX") if r.has_tag("GX") else None
            if gx is None and exon_idx:
                gx = infer_gene_by_overlap(chrom, blocks, exon_idx) or "NA"
                # 可选：写回 GX 标签（不强制）
                # if gx and gx != "NA": r.set_tag("GX", gx, value_type='Z')
            key = (cb, gx or "NA", strand, sj)

        qlen = r.query_alignment_length or 0
        buckets[key].append((r.query_name, ur, r.mapping_quality, qlen))

    # 2) Per-bucket UMI correction & pick representative
    molecules_out = open(f"{args.out}.molecules.tsv","w")
    molecules_out.write("cell\tkey2\tstrand\tSJ_or_noSJ\tUB_corrected\tcount\tmembers\n")
    assign_out = open(f"{args.out}.assignments.tsv","w")
    assign_out.write("read_name\tcell\tkey2\tstrand\tSJ_or_noSJ\tUR_raw\tUB_corrected\tkept\n")

    for key, items in buckets.items():
        # key unpack: for --no-gene: (CB, strand, SJorNoSJ); else: (CB, GX, strand, SJ)
        if args.no_gene:
            cell, strand, sj_or_no = key
            key2 = ""  # 占位
        else:
            cell, key2, strand, sj_or_no = key

        cnt = collections.Counter([ur for (_n, ur, _mq, _ql) in items])

        if args.fast_h1 and args.ham == 1:
            rep = directional_adjacency_fast_h1(cnt, ratio=args.ratio)
        else:
            rep = directional_adjacency(cnt, ham_max=args.ham, ratio=args.ratio)

        # molecules.tsv
        members = defaultdict(list)
        for u,c in cnt.items():
            members[rep[u]].extend([u]*c)
        for ub_corr, mems in members.items():
            molecules_out.write(f"{cell}\t{key2}\t{strand}\t{sj_or_no}\t{ub_corr}\t{len(mems)}\t{','.join(mems)}\n")

        # pick representative per corrected UMI
        per_mol = defaultdict(list)  # ub_corr -> [(name, ur, mq, ql)]
        for (name, ur, mq, ql) in items:
            per_mol[rep[ur]].append((name, ur, mq, ql))
        for ub_corr, rs in per_mol.items():
            rs.sort(key=lambda t: (t[2], t[3]), reverse=True)  # MAPQ, qlen
            keep = rs[0][0]
            keep_names.add(keep)
            for (name, _ur, _mq, _ql) in rs:
                read2UB[name] = ub_corr

        # assignments.tsv
        for (name, ur, _mq, _ql) in items:
            kept = 1 if name in keep_names else 0
            ub_corr = rep[ur]
            assign_out.write(f"{name}\t{cell}\t{key2}\t{strand}\t{sj_or_no}\t{ur}\t{ub_corr}\t{kept}\n")

    molecules_out.close(); assign_out.close()

    # 3) Second pass: write BAM with UB (corrected) and DA
    bam.reset()
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        kept = (r.query_name in keep_names)
        if r.query_name in read2UB and r.has_tag("UR"):
            r.set_tag("UB", read2UB[r.query_name], value_type='Z')  # corrected UMI
        r.set_tag("DA", 1 if kept else 0, value_type='i')
        if args.write_all:
            out_bam.write(r)
        else:
            if kept:
                out_bam.write(r)

    out_bam.close(); bam.close()
    sys.stderr.write(f"[umi_correct_longreads] buckets={len(buckets)} kept_reads={len(keep_names)}\n")

if __name__ == "__main__":
    main()
