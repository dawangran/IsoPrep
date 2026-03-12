#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, gzip, io, os, re, sys, tempfile, shutil, multiprocessing as mp
from datetime import datetime
from typing import Tuple

# ---------- logging ----------
def log(msg: str):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write(f"[{ts}] {msg}\n"); sys.stderr.flush()

# ---------- file I/O ----------
def open_in_text(path: str, buffer_mb: int):
    bufsize = max(1, buffer_mb) * 1024 * 1024
    if path == "-":
        return io.TextIOWrapper(sys.stdin.buffer, encoding="utf-8", newline="")
    if path.endswith(".gz"):
        raw = open(path, "rb", buffering=bufsize)
        gz  = gzip.GzipFile(fileobj=raw, mode="rb")
        return io.TextIOWrapper(gz, encoding="utf-8", newline="")
    return open(path, "rt", buffering=bufsize)

def open_out_text(path: str, buffer_mb: int, gzip_level: int):
    bufsize = max(1, buffer_mb) * 1024 * 1024
    if path == "-":
        return io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", newline="")
    if path.endswith(".gz"):
        raw = open(path, "wb", buffering=bufsize)
        gz  = gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=gzip_level, mtime=0)
        return io.TextIOWrapper(gz, encoding="utf-8", newline="")
    return open(path, "wt", buffering=bufsize)

# ---------- read id helpers ----------
_id_strip_re = re.compile(r'(/1|/2)$')

def core_read_id(name_line: str) -> str:
    """
    提取 FASTQ 名称的核心ID：@后到首个空格；去掉 /1 或 /2。
    """
    s = name_line.strip()
    if s.startswith("@"):
        s = s[1:]
    s = s.split()[0]
    s = _id_strip_re.sub("", s)
    return s

# ---------- fastq iter ----------
def iter_fastq(path: str, buffer_mb: int):
    """
    逐条读 FASTQ，yield (name_line, seq, plus_line, qual)
    """
    with open_in_text(path, buffer_mb) as fh:
        while True:
            n = fh.readline()
            if not n:
                break
            s = fh.readline(); p = fh.readline(); q = fh.readline()
            if not q:
                break
            if not n.startswith("@") or not p.startswith("+"):
                continue
            yield n.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n")

# ---------- model parsing ----------
def parse_model_literal(model_str: str):
    m = model_str.strip().upper()
    if not m or any(c not in "BU" for c in m):
        raise ValueError("字面模型仅允许 'B' 与 'U'")
    blocks_B = []; blocks_U = []
    i = 0; n = len(m)
    while i < n:
        ch = m[i]; j = i + 1
        while j < n and m[j] == ch: j += 1
        ln = j - i
        (blocks_B if ch == 'B' else blocks_U).append((i, ln))
        i = j
    return blocks_B, blocks_U, n

def parse_model_compact(model_str: str):
    s = model_str.strip().upper()
    if not re.fullmatch(r'(?:[BU]\d+)+', s):
        raise ValueError("简写模型格式错误，示例：B20U10 或 B20U10B5")
    pos = 0; blocks_B = []; blocks_U = []
    for m in re.finditer(r'([BU])(\d+)', s):
        ch = m.group(1); ln = int(m.group(2))
        if ln <= 0: raise ValueError("模型长度必须为正整数")
        (blocks_B if ch == 'B' else blocks_U).append((pos, ln))
        pos += ln
    return blocks_B, blocks_U, pos

def parse_model(model_str: str):
    return parse_model_compact(model_str) if any(c.isdigit() for c in model_str) else parse_model_literal(model_str)

def slice_blocks(seq: str, blocks):
    if not blocks: return ""
    if len(blocks) == 1:
        st, ln = blocks[0]; return seq[st:st+ln]
    return "".join(seq[st:st+ln] for st, ln in blocks)

# ---------- valid list ----------
def load_valid_map(path: str, cb_len: int):
    m = {}
    with open(path, "r") as f:
        for ln in f:
            ln = ln.strip()
            if not ln: continue
            parts = ln.split("\t")
            if len(parts) < 2:
                parts = ln.split()
            if len(parts) < 2:
                continue
            cb = parts[0].strip().upper()
            cell = parts[1].strip()
            if len(cb) != cb_len:
                continue
            m[cb] = cell
    return m

# ---------- hashing ----------
def bucket_index(read_id: str, n_shards: int) -> int:
    # 简单稳定哈希；避免 Python hash 随机化：使用自定义FNV-like
    h = 2166136261
    for ch in read_id.encode('utf-8'):
        h ^= ch
        h = (h * 16777619) & 0xFFFFFFFF
    return h % n_shards

# ---------- pass 1: write R1 buckets (id -> cb,umi) ----------
def pass1_r1_to_buckets(r1_path: str, buffer_mb: int, cb_blocks, umi_blocks, model_len: int,
                        uppercase: bool, n_shards: int, tmpdir: str) -> Tuple[int,int]:
    # 打开每个桶的 TSV 文件： id \t cb \t umi
    writers = []
    for i in range(n_shards):
        fp = os.path.join(tmpdir, f"r1.{i:04d}.tsv")
        writers.append(open(fp, "wt", buffering=1024*1024))
    total = kept = 0
    for n, s, p, q in iter_fastq(r1_path, buffer_mb):
        total += 1
        if len(s) < model_len:
            continue
        window = s[:model_len]
        if uppercase: window = window.upper()
        cb  = slice_blocks(window, cb_blocks)
        umi = slice_blocks(window, umi_blocks)
        rid = core_read_id(n)
        b = bucket_index(rid, n_shards)
        writers[b].write(f"{rid}\t{cb}\t{umi}\n")
        kept += 1
        if kept % 1_000_000 == 0:
            log(f"R1 pass1 kept ~{kept} (total read ~{total})")
    for w in writers:
        w.close()
    return total, kept

# ---------- pass 2: write R2 buckets (id, name_wo_at, seq, qual) ----------
def pass2_r2_to_buckets(r2_path: str, buffer_mb: int, n_shards: int, tmpdir: str) -> int:
    writers = []
    for i in range(n_shards):
        fp = os.path.join(tmpdir, f"r2.{i:04d}.tsv")
        writers.append(open(fp, "wt", buffering=1024*1024))
    total = 0
    for n, s, p, q in iter_fastq(r2_path, buffer_mb):
        rid = core_read_id(n)
        name_wo_at = n[1:] if n.startswith("@") else n
        b = bucket_index(rid, n_shards)
        writers[b].write(f"{rid}\t{name_wo_at}\t{s}\t{q}\n")
        total += 1
        if total % 1_000_000 == 0:
            log(f"R2 pass2 written ~{total}")
    for w in writers:
        w.close()
    return total

# ---------- pass 3: join buckets and output ----------
def join_one_bucket(i: int, tmpdir: str, out_part: str, valid_map, sample: str,
                    gzip_level: int, buffer_mb: int) -> Tuple[str, int, int]:
    """
    读 r1.i.tsv -> dict; 流式扫描 r2.i.tsv，命中的输出到 out_part（FASTQ）。
    返回 (out_part_path, r1_ids, out_kept)
    """
    r1_path = os.path.join(tmpdir, f"r1.{i:04d}.tsv")
    r2_path = os.path.join(tmpdir, f"r2.{i:04d}.tsv")
    if not os.path.exists(r1_path):
        # 该桶没有 R1，直接返回空文件
        open(out_part, "wb").close()
        return out_part, 0, 0

    # 1) 加载 R1 映射
    r1_map = {}
    with open(r1_path, "rt", buffering=1024*1024) as fh:
        for ln in fh:
            rid, cb, umi = ln.rstrip("\n").split("\t")
            r1_map[rid] = (cb, umi)
    r1_ids = len(r1_map)

    # 2) 扫描 R2 并输出
    kept = 0
    with open(out_part, "wt", buffering=1024*1024) as out, open(r2_path, "rt", buffering=1024*1024) as fh:
        for ln in fh:
            rid, name_wo_at, s2, q2 = ln.rstrip("\n").split("\t")
            x = r1_map.get(rid)
            if x is None:
                continue
            cb, umi = x
            if valid_map is not None:
                cell = valid_map.get(cb)
                if cell is None:
                    continue
                name = f"@{name_wo_at}|||CB:Z:{cb}|||UR:Z:{umi}|||DB:Z:{sample}_{cell}"
            else:
                name = f"@{name_wo_at}|||CB:Z:{cb}|||UR:Z:{umi}"
            out.write(name); out.write("\n")
            out.write(s2);   out.write("\n+\n")
            out.write(q2);   out.write("\n")
            kept += 1
    return out_part, r1_ids, kept

# 顶层包装：供 multiprocessing 使用（lambda 不可被 pickle）
def _join_bucket_star(args):
    return join_one_bucket(*args)

def join_buckets_parallel(n_shards: int, tmpdir: str, out_path: str,
                          threads: int, valid_map, sample: str,
                          buffer_mb: int, gzip_level: int):
    # 为每个桶准备一个 part 输出文件
    part_files = [os.path.join(tmpdir, f"part.{i:04d}.fq") for i in range(n_shards)]
    args_list = [(i, tmpdir, part_files[i], valid_map, sample, gzip_level, buffer_mb) for i in range(n_shards)]

    total_r1_ids = 0
    total_kept   = 0

    if threads <= 1 or n_shards == 1:
        # 串行回退（便于调试 / 低开销）
        for a in args_list:
            _, r1_ids, kept = join_one_bucket(*a)
            total_r1_ids += r1_ids
            total_kept   += kept
    else:
        ctx = mp.get_context("fork") if hasattr(mp, "get_context") else mp
        pool = ctx.Pool(processes=threads)
        stats = []
        try:
            chunksz = max(1, len(args_list) // max(1, threads * 2))
            for res in pool.imap_unordered(_join_bucket_star, args_list, chunksize=chunksz):
                stats.append(res)
        finally:
            pool.close(); pool.join()
        total_r1_ids = sum(x[1] for x in stats)
        total_kept   = sum(x[2] for x in stats)

    # 合并 part → 最终输出
    out = open_out_text(out_path, buffer_mb, gzip_level)
    try:
        for pf in part_files:  # 固定按桶序写入，保证确定性
            with open(pf, "rt", buffering=1024*1024) as fh:
                shutil.copyfileobj(fh, out, length=1024*1024*8)
    finally:
        if out is not sys.stdout:
            out.close()

    # 清理 part
    for pf in part_files:
        try: os.remove(pf)
        except Exception: pass

    return total_r1_ids, total_kept

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(
        description="以 R1 为准做 read-id 连接：从 R1 切 CB/UMI，贴到匹配的 R2 readname，仅输出 R1 存在的 read id。"
    )
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--out", required=True, help="输出 FASTQ（可 .gz 或 '-'）")
    ap.add_argument("--model", required=True, help="如 B20U10 或 字面 BBBB...UUU...")
    ap.add_argument("--valid_list", help="CB<tab>细胞身份；提供则仅保留名单内 CB（需配合 --sample）")
    ap.add_argument("--sample", help="样本ID；与 --valid_list 同用，在 readname 写 DB:Z:<sample>_<细胞身份>")

    ap.add_argument("--threads", type=int, default=8, help="并行桶连接线程数（pass3）")
    ap.add_argument("--shards",  type=int, default=256, help="哈希分桶数（越大越省内存，占用更多小文件）")
    ap.add_argument("--buffer_mb", type=int, default=32, help="读写缓冲(MB)")
    ap.add_argument("--gzip_level", type=int, default=1, help="输出 .gz 压缩等级（1=最快）")
    ap.add_argument("--uppercase", type=int, default=1, help="是否将模型窗口转大写(1/0)")
    ap.add_argument("--tmpdir", default="./", help="临时目录")
    args = ap.parse_args()

    # parse model
    try:
        cb_blocks, umi_blocks, model_len = parse_model(args.model)
    except Exception as e:
        log(f"ERROR: 解析模型失败：{e}"); sys.exit(2)
    cb_len = sum(ln for _, ln in cb_blocks)
    umi_len = sum(ln for _, ln in umi_blocks)
    if cb_len == 0 or umi_len == 0:
        log("ERROR: 模型中必须同时包含 B（CB）与 U（UMI）区段。"); sys.exit(2)

    # valid list
    valid_map = None
    if args.valid_list:
        if not args.sample:
            log("ERROR: 提供 --valid_list 时必须同时提供 --sample"); sys.exit(2)
        valid_map = load_valid_map(args.valid_list, cb_len)
        if not valid_map:
            log("WARN: 有效名单为空或没有匹配长度的 CB，读段将全部被过滤。")

    # temp dir
    if args.tmpdir:
        base_tmpdir = os.path.abspath(args.tmpdir)
        os.makedirs(base_tmpdir, exist_ok=True)
        if not os.path.isdir(base_tmpdir):
            log(f"ERROR: tmpdir 不是目录：{base_tmpdir}"); sys.exit(2)
        tmpdir = tempfile.mkdtemp(prefix="r1join_tmp_", dir=base_tmpdir)
    else:
        tmpdir = tempfile.mkdtemp(prefix="r1join_tmp_")
    log(f"Temp dir: {tmpdir}")
    n_shards = max(1, int(args.shards))

    # pass1: R1 → buckets
    log("PASS1: Scanning R1 and writing hash buckets (id -> CB,UMI)...")
    r1_total, r1_kept = pass1_r1_to_buckets(
        args.r1, args.buffer_mb, cb_blocks, umi_blocks, model_len, bool(args.uppercase),
        n_shards, tmpdir
    )
    log(f"PASS1 done. R1 total={r1_total}, usable={r1_kept} (len(seq)>=model_len)")

    # pass2: R2 → buckets
    log("PASS2: Scanning R2 and writing hash buckets (id -> rec)...")
    r2_total = pass2_r2_to_buckets(args.r2, args.buffer_mb, n_shards, tmpdir)
    log(f"PASS2 done. R2 total={r2_total}")

    # pass3: join per-bucket
    log("PASS3: Joining buckets and writing output...")
    total_r1_ids, total_kept = join_buckets_parallel(
        n_shards, tmpdir, args.out, args.threads, valid_map, args.sample,
        args.buffer_mb, args.gzip_level
    )
    log(f"PASS3 done. Buckets joined: r1_ids={total_r1_ids}, kept_out={total_kept}")

    # cleanup buckets
    for i in range(n_shards):
        for kind in ("r1", "r2"):
            fp = os.path.join(tmpdir, f"{kind}.{i:04d}.tsv")
            try: os.remove(fp)
            except Exception: pass
    try:
        if args.tmpdir is None:
            os.rmdir(tmpdir)
    except Exception:
        pass

    ratio = (total_kept / r1_kept) if r1_kept else 0.0
    log(f"Done. R1_usable={r1_kept}, Output_kept={total_kept} (join_ratio={ratio:.4f}).")
    if args.out != "-":
        with open(args.out + ".summary.txt", "w") as sf:
            sf.write("SUMMARY\n")
            sf.write(f"R1_TOTAL\t{r1_total}\n")
            sf.write(f"R1_USABLE(model_ok)\t{r1_kept}\n")
            sf.write(f"R2_TOTAL\t{r2_total}\n")
            sf.write(f"R1_IDS_IN_BUCKETS\t{total_r1_ids}\n")
            sf.write(f"OUT_KEPT(joined)\t{total_kept}\n")
            sf.write(f"JOIN_RATIO(out/r1_usable)\t{ratio:.4f}\n")
            if valid_map is not None:
                sf.write(f"SAMPLE\t{args.sample}\n")
                sf.write(f"MODEL\t{args.model}\n")

if __name__ == "__main__":
    main()
