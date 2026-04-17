#!/usr/bin/env bash
# Offline parity test: run bin/qcsum_metrics.pl against a frozen Picard
# HsMetrics fixture and diff the emitted qcsum.txt against a golden copy.
#
# This is the only end-to-end check we can run without Dragen/Picard/MultiQC
# installed. Since qcsum_metrics.pl is copied verbatim from ucladx/bcl-qc
# HEAD, byte-identical output here proves per-sample QC parity with upstream.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
FIX="$REPO/tests/fixtures"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

cp "$FIX/sample.hsm.txt" "$WORK/sample.hsm.txt"

perl "$REPO/bin/qcsum_metrics.pl" \
    --prefix sample \
    --qcfolder "$WORK" \
    --pipeline_version PCPTNDx-v2.0 \
    --platform NovaSeq6000 \
    --pass_min_align_pct 49 --fail_min_align_pct 49 \
    --covered 250 \
    --pass_min_roi_pct 39 --fail_min_roi_pct 39 \
    --pass_min_avgcov 399 --fail_min_avgcov 399 \
    --pass_min_reads 99999999 --fail_min_reads 99999999 \
    --capture twist --capture_version v1.0

if ! diff -u "$FIX/sample.qcsum.expected.txt" "$WORK/sample.qcsum.txt"; then
    echo "FAIL: qcsum output diverged from golden fixture" >&2
    exit 1
fi

# Also exercise the aggregator on a single input
cp "$WORK/sample.qcsum.txt" "$WORK/sample1.qcsum.txt"
cp "$WORK/sample.qcsum.txt" "$WORK/sample2.qcsum.txt"
python3 "$REPO/bin/aggregate_qcsum.py" --out "$WORK/qcsum_mqc.csv" \
    "$WORK/sample1.qcsum.txt" "$WORK/sample2.qcsum.txt"

# qcsum_mqc.csv must have: header, N data rows, trailing newline
lines=$(wc -l < "$WORK/qcsum_mqc.csv")
[ "$lines" -eq 4 ] || { echo "FAIL: expected 4 lines in qcsum_mqc.csv, got $lines" >&2; cat "$WORK/qcsum_mqc.csv" >&2; exit 1; }
head -1 "$WORK/qcsum_mqc.csv" | diff -q - "$REPO/assets/qcsum_mqc_header.txt"

echo "PASS: qcsum parity + aggregator"
