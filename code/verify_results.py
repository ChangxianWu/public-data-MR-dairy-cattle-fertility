#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify reproduced result tables against the published checksums.

Usage: python code/verify_results.py
Exits 0 with ALL CHECKSUMS MATCH when every table listed in results/CHECKSUMS.txt
reproduces byte-for-byte from the public raw data.
"""
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
LIST = HERE / "results" / "CHECKSUMS.txt"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def main():
    if not LIST.exists():
        print("results/CHECKSUMS.txt not found")
        return 2
    bad = total = 0
    for line in LIST.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        h, name = line.split(maxsplit=1)
        fp = HERE / "results" / name.strip()
        total += 1
        if not fp.exists():
            print(f"  MISSING  {name.strip()}")
            bad += 1
        elif sha256(fp) != h:
            print(f"  MISMATCH {name.strip()}")
            bad += 1
    print(f"{total - bad}/{total} tables match the published checksums")
    if bad == 0 and total:
        print("ALL CHECKSUMS MATCH")
        return 0
    print("MISMATCHES PRESENT - check data/raw snapshots against supplementary/data_manifest.json")
    return 1


if __name__ == "__main__":
    sys.exit(main())
