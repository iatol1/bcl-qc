"""Helpers for filtering Illumina v1/v2 SampleSheet CSV files by Sample_ID."""

import csv
import io
import os
from collections.abc import Collection, Iterable

DATA_SECTIONS = {"[BCLConvert_Data]", "[Data]"}


def _section_name(line: str) -> str | None:
    """Return the normalized section header (e.g. '[BCLConvert_Data]'), or None."""
    first_cell = line.strip().split(",", 1)[0].strip()
    if first_cell.startswith("[") and first_cell.endswith("]"):
        return first_cell
    return None


def filter_samplesheet_lines(
    lines: Iterable[str],
    samples: Collection[str],
) -> tuple[list[str], set[str]]:
    """Filter SampleSheet CSV lines to only rows whose Sample_ID is in ``samples``.

    Preserves all non-data sections and header rows verbatim. Supports both
    Illumina SampleSheet v2 (``[BCLConvert_Data]``) and v1 (``[Data]``).
    Returns ``(filtered_lines, matched_sample_ids)``.
    """
    requested = {s.strip() for s in samples if s and s.strip()}
    output: list[str] = []
    matched: set[str] = set()

    in_data_section = False
    sample_id_col: int | None = None

    for line in lines:
        section = _section_name(line)
        if section is not None:
            in_data_section = section in DATA_SECTIONS
            sample_id_col = None
            output.append(line)
            continue

        if not in_data_section:
            output.append(line)
            continue

        stripped = line.strip()
        if not stripped or all(ch == "," for ch in stripped):
            output.append(line)
            continue

        row = next(csv.reader(io.StringIO(stripped)))
        if sample_id_col is None:
            # First non-empty row inside the data section is the column header.
            lowered = [col.strip().lower() for col in row]
            if "sample_id" in lowered:
                sample_id_col = lowered.index("sample_id")
            else:
                sample_id_col = 0
            output.append(line)
            continue

        if sample_id_col < len(row):
            sample_id = row[sample_id_col].strip()
            if sample_id in requested:
                matched.add(sample_id)
                output.append(line)

    return output, matched


def write_filtered_samplesheet(
    src_path: str,
    dst_path: str,
    samples: Collection[str],
) -> set[str]:
    """Filter ``src_path`` to ``samples``; write ``dst_path`` only if any match."""
    with open(src_path, "r", newline="") as fin:
        filtered_lines, matched = filter_samplesheet_lines(fin, samples)
    if matched:
        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
        with open(dst_path, "w", newline="") as fout:
            fout.writelines(filtered_lines)
    return matched
