#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scan_cb_umi.py
--------------
CB/UMI scanner with optional whitelist-based barcode correction.

Built on the v2 + numba + multiprocessing workflow, with CB correction:
- Under the whitelist assumption (minimum Hamming distance of 3):
  * EXACT / HAM1: unique direct correction.
  * (Optional) HAM2 with globally unique minimum: enable with
    --cb_allow_ham2_unique.
- When correction is enabled (--cb_correct): write FASTQ only when both CB1
  and CB2 are corrected successfully
  (EXACT/HAM1/(optional)HAM2_MIN_UNIQ) and extraction succeeds.
- When correction is disabled: behavior matches previous logic
  (write FASTQ on extraction success).

Examples:
  Extraction only (no correction):
    python -m IsoPrep.bin.scan_cb_umi in.fastq.gz \
      --out_tsv out.tsv --out_fastq out.fastq.gz \
      --threads 8 --mp_chunksize 256 --batch 100000

  Enable correction (EXACT/HAM1):
    python -m IsoPrep.bin.scan_cb_umi in.fastq.gz \
      --out_tsv out.tsv --out_fastq out.fastq.gz \
      --threads 8 --cb_correct --whitelist whitelist.txt

  Enable correction + HAM2 unique-minimum rescue:
    python -m IsoPrep.bin.scan_cb_umi in.fastq.gz \
      --out_tsv out.tsv --out_fastq out.fastq.gz \
      --threads 8 --cb_correct --whitelist whitelist.txt --cb_allow_ham2_unique
"""

import argparse, gzip, sys, multiprocessing as mp
from datetime import datetime

# ---------- Optional Numba ----------
USE_NUMBA = True
try:
    import numba as nb
    import numpy as np
except Exception:
    USE_NUMBA = False

TAG1 = "CTACGATCCGACTTTCTGCG"
TAG2 = "CCTTCC"
TAG3 = "CGATG"

def log(msg: str):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write(f"[{ts}] {msg}\n"); sys.stderr.flush()

# ------------ FASTQ utils ------------
def open_maybe_gzip(path: str):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def parse_fastq(stream):
    while True:
        name = stream.readline()
        if not name:
            break
        seq  = stream.readline()
        plus = stream.readline()
        qual = stream.readline()
        if not qual:
            break
        if not name.startswith('@') or not plus.startswith('+'):
            continue
        yield name[1:].strip(), seq.strip().upper(), qual.strip()

# ------------ Encoding helpers for numba -------------
def s2a(s):
    if not USE_NUMBA:
        return None
    arr = np.empty(len(s), dtype=np.uint8)
    m = { 'A':1, 'C':2, 'G':3, 'T':4, 'N':5, '-':6 }
    for i,ch in enumerate(s):
        arr[i] = m.get(ch, 0)
    return arr

# ------------ Numba kernels ------------
if USE_NUMBA:
    @nb.njit(cache=True)
    def t_match_score_n(tchar, rchar, score_match, score_mismatch):
        # Template base 'N'(5) is wildcard; score as a match.
        return score_match if (tchar == 5 or tchar == rchar) else score_mismatch

    @nb.njit(cache=True)
    def nw_semiglobal_template_numba(read_sub, templ, score_match, score_mismatch, score_gap):
        n = read_sub.size
        m = templ.size
        NEG_INF = -10_000_000
        dp = np.empty((n+1, m+1), dtype=np.int32)
        tb = np.empty((n+1, m+1), dtype=np.int8)   # 0=diag,1=up,2=left

        # init
        for j in range(m+1):
            dp[0,j] = -10_000_000; tb[0,j] = 2
        dp[0,0] = 0
        for i in range(1, n+1):
            dp[i,0] = 0; tb[i,0] = 1
        for j in range(1, m+1):
            dp[0,j] = dp[0,j-1] + score_gap

        # fill
        for i in range(1, n+1):
            rc = read_sub[i-1]
            for j in range(1, m+1):
                tc = templ[j-1]
                s_diag = dp[i-1,j-1] + t_match_score_n(tc, rc, score_match, score_mismatch)
                s_up   = dp[i-1,j] + score_gap
                s_left = dp[i,j-1] + score_gap
                best = s_diag; bt = 0
                if s_up > best: best = s_up; bt = 1
                if s_left > best: best = s_left; bt = 2
                dp[i,j] = best; tb[i,j] = bt

        # best end on last column
        best = -10_000_000; best_i = 0
        for i in range(n+1):
            if dp[i,m] > best:
                best = dp[i,m]; best_i = i

        # traceback
        i = best_i; j = m
        path_i = []; path_j = []; path_bt = []
        while j>0 or i>0:
            bt = tb[i,j]
            path_i.append(i); path_j.append(j); path_bt.append(bt)
            if bt == 0: i -= 1; j -= 1
            elif bt == 1: i -= 1
            else: j -= 1
            if i==0 and j==0: break
        path_i = path_i[::-1]; path_j = path_j[::-1]; path_bt = path_bt[::-1]
        return best_i, np.array(path_i, np.int32), np.array(path_j, np.int32), np.array(path_bt, np.int8)

    @nb.njit(cache=True)
    def sw_local_best_numba(seq, ref, max_k):
        n = seq.size; m = ref.size
        H = np.zeros((n+1, m+1), dtype=np.int32)
        ptr = np.zeros((n+1, m+1), dtype=np.int8)  # 0=stop,1=diag,2=up,3=left
        best_score = 0; bi=0; bj=0
        for i in range(1, n+1):
            si = seq[i-1]
            for j in range(1, m+1):
                rj = ref[j-1]
                diag = H[i-1,j-1] + (2 if si==rj else -1)
                up   = H[i-1,j] + (-2)
                left = H[i,j-1] + (-2)
                sc   = 0
                if diag > sc: sc = diag; ptr[i,j] = 1
                if up   > sc: sc = up;   ptr[i,j] = 2
                if left > sc: sc = left; ptr[i,j] = 3
                H[i,j] = sc
                if sc > best_score:
                    best_score = sc; bi=i; bj=j
        if best_score == 0:
            return -1, -1, -1
        i = bi; j = bj; edits=0
        while i>0 and j>0 and ptr[i,j] != 0:
            p = ptr[i,j]
            if p == 1:
                if seq[i-1] != ref[j-1]: edits += 1
                i -= 1; j -= 1
            elif p == 2:
                edits += 1; i -= 1
            else:
                edits += 1; j -= 1
        if edits > max_k:
            return -1, -1, -1
        return i, bi, edits

# ------------ Pure-Python fallback ------------
def sw_local_best_py(seq: str, ref: str, max_k: int):
    n, m = len(seq), len(ref)
    H = [[0]*(m+1) for _ in range(n+1)]
    ptr = [[0]*(m+1) for _ in range(n+1)]
    best = (0, 0, 0)
    for i in range(1, n+1):
        si = seq[i-1]
        for j in range(1, m+1):
            rj = ref[j-1]
            diag = H[i-1][j-1] + (2 if si==rj else -1)
            up   = H[i-1][j] + (-2)
            left = H[i][j-1] + (-2)
            sc   = max(0, diag, up, left)
            H[i][j] = sc
            if sc == 0: ptr[i][j] = 0
            elif sc == diag: ptr[i][j] = 1
            elif sc == up:   ptr[i][j] = 2
            else:            ptr[i][j] = 3
            if sc > best[0]: best = (sc, i, j)
    score, i, j = best
    if score == 0: return None
    i_end = i; edits = 0
    while i>0 and j>0 and ptr[i][j] != 0:
        p = ptr[i][j]
        if p == 1:
            if seq[i-1] != ref[j-1]: edits += 1
            i -= 1; j -= 1
        elif p == 2:
            edits += 1; i -= 1
        else:
            edits += 1; j -= 1
    i_start = i
    if edits > max_k: return None
    return (i_start, i_end, edits)

def nw_semiglobal_template_py(read_sub: str, template: str, score_gap=-2):
    n = len(read_sub); m = len(template)
    NEG_INF = -10**9
    dp = [[NEG_INF]*(m+1) for _ in range(n+1)]
    tb = [[0]*(m+1) for _ in range(n+1)]
    dp[0][0] = 0
    for i in range(1, n+1):
        dp[i][0] = 0; tb[i][0] = 1
    for j in range(1, m+1):
        dp[0][j] = dp[0][j-1] + score_gap; tb[0][j] = 2
    def t_match_score(tc, rc):
        return 2 if (tc=='N' or tc==rc) else -2
    for i in range(1, n+1):
        rc = read_sub[i-1]
        for j in range(1, m+1):
            tc = template[j-1]
            s_diag = dp[i-1][j-1] + t_match_score(tc, rc)
            s_up   = dp[i-1][j]   + score_gap
            s_left = dp[i][j-1]   + score_gap
            best = s_diag; bt=0
            if s_up > best: best = s_up; bt=1
            if s_left > best: best = s_left; bt=2
            dp[i][j] = best; tb[i][j] = bt
    best_i = 0; best = -10**9
    for i in range(0, n+1):
        if dp[i][m] > best:
            best = dp[i][m]; best_i = i
    i = best_i; j = m
    align_read = []
    align_temp = []
    while j > 0 or i > 0:
        bt = tb[i][j]
        if bt == 0:
            align_read.append(read_sub[i-1] if i>0 else '-')
            align_temp.append(template[j-1] if j>0 else '-')
            i -= 1; j -= 1
        elif bt == 1:
            align_read.append(read_sub[i-1] if i>0 else '-')
            align_temp.append('-')
            i -= 1
        else:
            align_read.append('-')
            align_temp.append(template[j-1] if j>0 else '-')
            j -= 1
        if i==0 and j==0:
            break
    align_read = align_read[::-1]; align_temp = align_temp[::-1]
    return align_read, align_temp, best_i

# ------------ Edit counting & N-block ------------
def edits_in_range(a_read, a_temp, j_start, j_end):
    edits = 0
    temp_idx = -1
    for r,t in zip(a_read, a_temp):
        if t != '-':
            temp_idx += 1
            in_seg = (j_start <= temp_idx < j_end)
        else:
            in_seg = False
        if in_seg:
            if r == '-' or t == '-':
                edits += 1
            else:
                if (t != 'N') and (r != t):
                    edits += 1
    return edits

def collect_segment_nblock(a_read, a_temp, j_start, j_end, read_offset, raw_seq):
    seq = []
    pos_first = None
    temp_idx = -1
    in_block = False
    consumed_read = 0
    for k, (r, t) in enumerate(zip(a_read, a_temp)):
        if r != '-':
            consumed_read += 1
        if t != '-':
            temp_idx += 1
            if temp_idx == j_start:
                in_block = True
            elif temp_idx == j_end:
                in_block = False
        if in_block and r != '-':
            if pos_first is None:
                pos_first = read_offset + consumed_read - 1
            seq.append(r)
            if len(seq) == 10:
                return ''.join(seq), pos_first
    # Fallback: if alignment yields <10 bp, append downstream read bases after N block.
    block_cols = []
    temp_idx = -1
    for k, t in enumerate(a_temp):
        if t != '-':
            temp_idx += 1
        if j_start <= temp_idx < j_end:
            block_cols.append(k)
    if block_cols:
        last_k = block_cols[-1]
        for k in range(last_k+1, len(a_read)):
            r = a_read[k]
            if r != '-':
                consumed_read += 1
                if pos_first is None:
                    pos_first = read_offset + consumed_read - 1
                seq.append(r)
                if len(seq) == 10:
                    return ''.join(seq), pos_first
    # Final fallback: if pos_first is known, slice 10 bp directly from raw sequence.
    if pos_first is not None and pos_first + 10 <= len(raw_seq):
        start = pos_first
        return raw_seq[start:start+10], start
    return None, None

# ------------ Whitelist correction ------------
def load_whitelist(path):
    wl = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip().upper()
            if not s: continue
            if len(s) != 10: continue
            wl.append(s)
    return wl

def hamdist(a, b):
    d = 0
    for i in range(10):
        if a[i] != b[i]:
            d += 1
    return d

def correct_cb_10mer(obs, wl, allow_ham2_unique=False):
    """
    Return (corr, status, dmin).
    status ∈ {"EXACT","HAM1","HAM2_MIN_UNIQ","UNCORR"}
    """
    if obs is None or len(obs) != 10:
        return None, "UNCORR", None
    best = None
    best_d = 11
    best_count = 0
    for w in wl:
        d = hamdist(obs, w)
        if d < best_d:
            best_d = d; best = w; best_count = 1
        elif d == best_d:
            best_count += 1
        if d == 0:
            return w, "EXACT", 0
    if best_d == 1 and best_count == 1:
        return best, "HAM1", 1
    if allow_ham2_unique and best_d == 2 and best_count == 1:
        return best, "HAM2_MIN_UNIQ", 2
    return None, "UNCORR", best_d

def is_corr_success(status: str) -> bool:
    return status in ("EXACT", "HAM1", "HAM2_MIN_UNIQ")

# ------------ Worker ------------
def worker(task):
    idx, rid, seq, qual, conf = task
    maxk1 = conf['maxk1']; score_match = conf['score_match']
    score_mismatch = conf['score_mismatch']; score_gap = conf['score_gap']
    T = conf['T']; T_arr = conf['T_arr']
    j0=conf['j0']; j1=conf['j1']; j2=conf['j2']; j3=conf['j3']; j4=conf['j4']; j5=conf['j5']; j6=conf['j6']; j7=conf['j7']
    do_corr = conf['do_corr']; wl = conf['whitelist']; allow_ham2 = conf['allow_ham2']

    # 1) anchor TAG1
    if USE_NUMBA:
        s_arr = s2a(seq); t1_arr = s2a(TAG1)
        s1,e1,ed = sw_local_best_numba(s_arr, t1_arr, maxk1)
        if s1 == -1:
            return (idx, rid, "NO_TAG1",
                    None, None, None, None, None,   # cb1_raw,pos,corr,status,dmin
                    None, None, None, None, None,   # cb2_...
                    None, None,                     # umi, umi_pos
                    None, None, None, None,         # ed1..ed4
                    None)                           # fq_rec
        start = s1
        read_sub = seq[start:]
        rs_arr = s2a(read_sub)
        _best_i, path_i, path_j, path_bt = nw_semiglobal_template_numba(rs_arr, T_arr, score_match, score_mismatch, score_gap)
        # reconstruct for downstream simplicity
        ar = []; at = []
        i = 0; j = 0
        for k in range(path_bt.size):
            bt = path_bt[k]
            if bt == 0:
                ar.append(read_sub[i]); at.append(T[j]); i+=1; j+=1
            elif bt == 1:
                ar.append(read_sub[i]); at.append('-'); i+=1
            else:
                ar.append('-'); at.append(T[j]); j+=1
        ar = "".join(ar); at = "".join(at)
    else:
        hit1 = sw_local_best_py(seq, TAG1, maxk1)
        if not hit1:
            return (idx, rid, "NO_TAG1",
                    None, None, None, None, None,
                    None, None, None, None, None,
                    None, None,
                    None, None, None, None,
                    None)
        s1,e1,_ = hit1
        start = s1
        read_sub = seq[start:]
        ar, at, _ = nw_semiglobal_template_py(read_sub, T, score_gap)

    # 2) edits (for QC)
    ed1 = edits_in_range(ar, at, j0, j1)
    ed2 = edits_in_range(ar, at, j2, j3)
    ed3 = edits_in_range(ar, at, j4, j5)
    ed4 = edits_in_range(ar, at, j6, j7)

    # 3) extract segments
    cb1_raw, cb1_pos = collect_segment_nblock(ar, at, j1, j2, read_offset=start, raw_seq=seq)
    cb2_raw, cb2_pos = collect_segment_nblock(ar, at, j3, j4, read_offset=start, raw_seq=seq)
    umi,     umi_pos = collect_segment_nblock(ar, at, j5, j6, read_offset=start, raw_seq=seq)
    if (cb1_raw is None) or (cb2_raw is None) or (umi is None):
        return (idx, rid, "EXTRACTION_FALLBACK_FAILED",
                cb1_raw, cb1_pos, None, "UNCORR", None,
                cb2_raw, cb2_pos, None, "UNCORR", None,
                umi, umi_pos,
                ed1, ed2, ed3, ed4,
                None)

    # 4) correction (optional)
    if do_corr and wl is not None:
        cb1_corr, cb1_status, cb1_dmin = correct_cb_10mer(cb1_raw, wl, allow_ham2)
        cb2_corr, cb2_status, cb2_dmin = correct_cb_10mer(cb2_raw, wl, allow_ham2)
    else:
        cb1_corr, cb1_status, cb1_dmin = cb1_raw, "BYPASS", 0
        cb2_corr, cb2_status, cb2_dmin = cb2_raw, "BYPASS", 0

    # 5) QUAL slicing
    def slice_q(pos):
        if pos is None: return None
        end = pos + 10
        return None if (pos < 0 or end > len(qual)) else qual[pos:end]
    q1 = slice_q(cb1_pos); q2 = slice_q(cb2_pos); q3 = slice_q(umi_pos)
    if (q1 is None) or (q2 is None) or (q3 is None):
        return (idx, rid, "QUAL_FAIL",
                cb1_raw, cb1_pos, cb1_corr, cb1_status, cb1_dmin,
                cb2_raw, cb2_pos, cb2_corr, cb2_status, cb2_dmin,
                umi, umi_pos,
                ed1, ed2, ed3, ed4,
                None)

    # 6) Build FASTQ SEQ (prefer corrected CB; fallback to raw when correction is unavailable).
    c1 = cb1_corr if cb1_corr is not None else cb1_raw
    c2 = cb2_corr if cb2_corr is not None else cb2_raw
    concat = f"{c1}{c2}{umi}"
    fq_rec = f"@{rid}\n{concat}\n+\n{q1+q2+q3}\n"

    return (idx, rid, "OK",
            cb1_raw, cb1_pos, cb1_corr, cb1_status, cb1_dmin,
            cb2_raw, cb2_pos, cb2_corr, cb2_status, cb2_dmin,
            umi, umi_pos,
            ed1, ed2, ed3, ed4,
            fq_rec)

# ------------ Main ------------
def main():
    ap = argparse.ArgumentParser(description="Template-align CB/UMI extractor v2 + Numba + MP + optional CB correction.")
    ap.add_argument("fastq", help="FASTQ path or '-' for stdin")
    ap.add_argument("--out_tsv", required=True, help="TSV output path")
    ap.add_argument("--out_fastq", required=True, help="FASTQ output path or '-' for stdout")
    ap.add_argument("--polyT_min", type=int, default=4)
    ap.add_argument("--maxk1", type=int, default=3)
    ap.add_argument("--score_match", type=int, default=2)
    ap.add_argument("--score_mismatch", type=int, default=-2)
    ap.add_argument("--score_gap", type=int, default=-2)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--mp_chunksize", type=int, default=256)
    ap.add_argument("--batch", type=int, default=100000)

    # correction options
    ap.add_argument("--cb_correct", action="store_true", help="Enable CB1/CB2 correction.")
    ap.add_argument("--whitelist", type=str, default=None, help="Whitelist path; one 10bp barcode per line (A/C/G/T).")
    ap.add_argument("--cb_allow_ham2_unique", action="store_true",
                    help="Allow Hamming=2 rescue when the minimum is unique.")

    args = ap.parse_args()

    if args.cb_correct and not args.whitelist:
        log("ERROR: --whitelist is required when --cb_correct is enabled"); sys.exit(2)

    whitelist = None
    if args.cb_correct:
        whitelist = load_whitelist(args.whitelist)
        if len(whitelist) == 0:
            log("ERROR: whitelist is empty or has no valid 10bp entries"); sys.exit(2)
        log(f"Whitelist loaded: {len(whitelist)} entries.")

    poly = "T"*args.polyT_min
    T = TAG1 + "N"*10 + TAG2 + "N"*10 + TAG3 + "N"*10 + poly
    T_arr = s2a(T) if USE_NUMBA else None

    # template indices
    j0=0; j1=len(TAG1); j2=j1+10
    j3=j2+len(TAG2); j4=j3+10
    j5=j4+len(TAG3); j6=j5+10
    j7=j6+len(poly)

    conf = dict(
        maxk1=args.maxk1,
        score_match=args.score_match,
        score_mismatch=args.score_mismatch,
        score_gap=args.score_gap,
        T=T, T_arr=T_arr,
        j0=j0,j1=j1,j2=j2,j3=j3,j4=j4,j5=j5,j6=j6,j7=j7,
        do_corr=args.cb_correct,
        whitelist=whitelist,
        allow_ham2=args.cb_allow_ham2_unique
    )

    # TSV header (includes correction fields).
    header = "\t".join([
        "read_id",
        "cb1_raw","cb1_pos","cb1_corr","cb1_status","cb1_dmin",
        "cb2_raw","cb2_pos","cb2_corr","cb2_status","cb2_dmin",
        "umi","umi_pos",
        "ed_tag1","ed_tag2","ed_tag3","ed_tag4","total_edits",
        "proc_status"
    ])

    out_tsv = open(args.out_tsv, "w"); out_tsv.write(header + "\n")
    if args.out_fastq == "-":
        out_fq = sys.stdout
    else:
        out_fq = gzip.open(args.out_fastq, "wt") if args.out_fastq.endswith(".gz") else open(args.out_fastq, "w")

    total = kept = 0
    stat = {"OK":0, "NO_TAG1":0, "EXTRACTION_FALLBACK_FAILED":0, "QUAL_FAIL":0}

    # Correction stats (used only when --cb_correct is enabled).
    if args.cb_correct:
        for k in ["CORR_CB1_EXACT","CORR_CB1_HAM1","CORR_CB1_HAM2","CORR_CB1_UNCORR",
                  "CORR_CB2_EXACT","CORR_CB2_HAM1","CORR_CB2_HAM2","CORR_CB2_UNCORR",
                  "CORR_BOTH_OK","CORR_ANY_FAIL"]:
            stat[k] = 0

    log(f"Numba={'ON' if USE_NUMBA else 'OFF'}, threads={args.threads}, chunksize={args.mp_chunksize}, batch={args.batch}")
    with open_maybe_gzip(args.fastq) as fh, mp.get_context("fork").Pool(processes=args.threads) as pool:
        batch = []; idx0 = 0
        def handle_result(res):
            nonlocal kept
            (_idx, rid, status,
             cb1_raw, cb1_pos, cb1_corr, cb1_status, cb1_dmin,
             cb2_raw, cb2_pos, cb2_corr, cb2_status, cb2_dmin,
             umi, umi_pos,
             ed1, ed2, ed3, ed4,
             fq_rec) = res

            total_ed = 0
            for x in (ed1,ed2,ed3,ed4):
                total_ed += (x or 0)

            # Count correction-status distribution (only when correction is enabled).
            if args.cb_correct:
                # CB1
                if cb1_status == "EXACT": stat["CORR_CB1_EXACT"] += 1
                elif cb1_status == "HAM1": stat["CORR_CB1_HAM1"] += 1
                elif cb1_status == "HAM2_MIN_UNIQ": stat["CORR_CB1_HAM2"] += 1
                else: stat["CORR_CB1_UNCORR"] += 1
                # CB2
                if cb2_status == "EXACT": stat["CORR_CB2_EXACT"] += 1
                elif cb2_status == "HAM1": stat["CORR_CB2_HAM1"] += 1
                elif cb2_status == "HAM2_MIN_UNIQ": stat["CORR_CB2_HAM2"] += 1
                else: stat["CORR_CB2_UNCORR"] += 1

            # TSV row
            out_tsv.write("\t".join([
                rid,
                cb1_raw or "", str(cb1_pos) if cb1_pos is not None else "",
                cb1_corr or "", cb1_status or "", "" if cb1_dmin is None else str(cb1_dmin),
                cb2_raw or "", str(cb2_pos) if cb2_pos is not None else "",
                cb2_corr or "", cb2_status or "", "" if cb2_dmin is None else str(cb2_dmin),
                umi or "", str(umi_pos) if umi_pos is not None else "",
                str(ed1) if ed1 is not None else "", str(ed2) if ed2 is not None else "",
                str(ed3) if ed3 is not None else "", str(ed4) if ed4 is not None else "",
                str(total_ed),
                status
            ]) + "\n")

            stat[status] = stat.get(status, 0) + 1

            # Whether to write FASTQ
            write_fastq = False
            if status == "OK" and fq_rec is not None:
                if args.cb_correct:
                    if is_corr_success(cb1_status) and is_corr_success(cb2_status):
                        write_fastq = True
                        stat["CORR_BOTH_OK"] += 1
                    else:
                        stat["CORR_ANY_FAIL"] += 1
                else:
                    write_fastq = True

            if write_fastq:
                out_fq.write(fq_rec); kept += 1

        # Main loop
        for rid, seq, qual in parse_fastq(fh):
            batch.append((idx0, rid, seq, qual, conf)); idx0 += 1; total += 1
            if len(batch) >= args.batch:
                for res in pool.imap(worker, batch, chunksize=args.mp_chunksize):
                    handle_result(res)
                log(f"Processed {idx0} reads... (FASTQ kept {kept})")
                batch.clear()
        if batch:
            for res in pool.imap(worker, batch, chunksize=args.mp_chunksize):
                handle_result(res)

    ratio = (kept/total) if total>0 else 0.0
    log(f"Done. Total reads: {total}, FASTQ output: {kept} ({ratio:.4f}).")
    log("Status counts: " + ", ".join([f"{k}={stat.get(k,0)}" for k in ["OK","QUAL_FAIL","EXTRACTION_FALLBACK_FAILED","NO_TAG1"]]))
    if args.cb_correct:
        log("Correction stats: " +
            ", ".join([f"{k}={stat.get(k,0)}" for k in
                       ["CORR_CB1_EXACT","CORR_CB1_HAM1","CORR_CB1_HAM2","CORR_CB1_UNCORR",
                        "CORR_CB2_EXACT","CORR_CB2_HAM1","CORR_CB2_HAM2","CORR_CB2_UNCORR",
                        "CORR_BOTH_OK","CORR_ANY_FAIL"]]))

    # SUMMARY
    summ_path = args.out_tsv + ".summary.txt"
    with open(summ_path, "w") as sf:
        sf.write("SUMMARY\n")
        sf.write(f"INPUT\t{total}\n")
        sf.write(f"FASTQ_OUTPUT\t{kept}\n")
        sf.write(f"RATIO\t{ratio:.4f}\n")
        for k in ["OK","QUAL_FAIL","EXTRACTION_FALLBACK_FAILED","NO_TAG1"]:
            sf.write(f"{k}\t{stat.get(k,0)}\n")
        if args.cb_correct:
            sf.write("CORRECTION_STATS\n")
            sf.write(f"CB1_EXACT\t{stat.get('CORR_CB1_EXACT',0)}\n")
            sf.write(f"CB1_HAM1\t{stat.get('CORR_CB1_HAM1',0)}\n")
            sf.write(f"CB1_HAM2\t{stat.get('CORR_CB1_HAM2',0)}\n")
            sf.write(f"CB1_UNCORR\t{stat.get('CORR_CB1_UNCORR',0)}\n")
            sf.write(f"CB2_EXACT\t{stat.get('CORR_CB2_EXACT',0)}\n")
            sf.write(f"CB2_HAM1\t{stat.get('CORR_CB2_HAM1',0)}\n")
            sf.write(f"CB2_HAM2\t{stat.get('CORR_CB2_HAM2',0)}\n")
            sf.write(f"CB2_UNCORR\t{stat.get('CORR_CB2_UNCORR',0)}\n")
            sf.write(f"CORR_BOTH_OK\t{stat.get('CORR_BOTH_OK',0)}\n")
            sf.write(f"CORR_ANY_FAIL\t{stat.get('CORR_ANY_FAIL',0)}\n")
    log(f"Summary written to {summ_path}")

if __name__ == "__main__":
    main()
