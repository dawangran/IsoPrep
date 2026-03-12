#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NO-UMI build: umi_correct_longreads.py is neutralized (keeps CLI stable)."""
import sys, argparse
def main():
    ap = argparse.ArgumentParser(description="NO-UMI passthrough (does nothing)")
    ap.add_argument("--in", dest="inp", default="-")
    ap.add_argument("--out", dest="out", default="-")
    ap.add_argument("--whitelist"); ap.add_argument("--allow-ham2-unique", action="store_true")
    ap.add_argument("--threads", type=int, default=1)
    _ = ap.parse_args()
    sys.stdout.buffer.write(sys.stdin.buffer.read())
if __name__ == "__main__":
    main()
