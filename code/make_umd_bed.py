#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert the Dryad EMMAX outcome CSV into the BED coordinate file used by the pipeline.

The IndexP60 EMMAX release (Dryad doi:10.5061/dryad.0gb5mkm04, file
0gb5mkm04__File_S1_EMMAX_RI_extracted.csv) lists one probe per row with
Chromosome/Position on UMD3.1.1; unplaced probes (chromosome 0) carry
placeholder coordinates. This script keeps the genome-placed probes in file
order and writes the 0-based BED that liftover_chain.py consumes:

    chr{chromosome}    {position-1}    {position}    {running_index}

The fourth column is the 0-based index of the row among the kept (placed)
rows; run_mr_main.py uses it to join BED positions back to the CSV statistics.

Usage:
    python code/make_umd_bed.py <emmax_csv> <output.bed>

On Windows the output is written in text mode, matching the CRLF line endings
of the original analysis (the file is consumed internally, so the exact
line-ending convention only matters for byte-level verification).
"""
import csv
import sys

UCSC_MIN, UCSC_MAX = 1, 29          # cattle autosomes in the source release


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    src, dst = argv[1], argv[2]
    kept = 0
    with open(src, encoding="utf-8", errors="replace") as fi, open(dst, "w") as fo:
        reader = csv.reader(fi)
        header = next(reader)
        i_c = header.index("Chromosome")
        i_p = header.index("Position")
        for row in reader:
            try:
                c = int(row[i_c])
                p = int(row[i_p])
            except (ValueError, IndexError):
                continue
            if not (UCSC_MIN <= c <= UCSC_MAX and p >= 1):
                continue
            fo.write(f"chr{c}\t{p - 1}\t{p}\t{kept}\n")
            kept += 1
    print(f"wrote {kept} placed probes -> {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
