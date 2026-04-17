#!/usr/bin/env python3
"""Concatenate per-sample *.qcsum.txt files into one qcsum_mqc.csv for MultiQC.

Mirrors the aggregation block from ucladx/bcl-qc bclqc.py qcsum() (lines 312-318):
write the canonical header once, then every sample's data row (skipping each
file's own header).
"""
import argparse
import sys
from pathlib import Path

HEADER = (
    "Sample,Sequencing_Platform,Pipeline_version,Alignment_QC,Coverage_QC,"
    "Total_Reads,%Reads_Aligned,Capture,Avg_Capture_Coverage,%On/Near_Bait_Bases,"
    "%On_Bait_Bases,FOLD_80_BASE_PENALTY,Avg_ROI_Coverage,MEDIAN_ROI_COVERAGE,"
    "MAX_ROI_COVERAGE,%ROI_1x,%ROI_20x,%ROI_100x,%ROI_250x,%ROI_500x"
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="Output qcsum_mqc.csv path")
    ap.add_argument("qcsum_files", nargs="+", help="Per-sample *.qcsum.txt files")
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(HEADER + "\n")
        for f in sorted(args.qcsum_files):
            lines = Path(f).read_text().splitlines(keepends=True)
            if len(lines) < 2:
                print(f"warning: {f} has no data row, skipping", file=sys.stderr)
                continue
            fh.writelines(lines[1:])
        fh.write("\n")


if __name__ == "__main__":
    main()
