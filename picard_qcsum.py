#!/usr/bin/env python3
"""
Standalone CLI and library for running Picard CollectHsMetrics and generating
panel-aware QC summary (qcsum) reports.

Usage:
  # Run on a PCP tumor sample
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel pcp_tumor

  # Run on a Heme sample with custom sample ID
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel heme_comp --sample-id MySample

  # Run on a Heme subpanel (MPN)
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel heme_mpn

  # Output to a specific directory
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel pcp_normal --output-dir /path/to/output

  # List available panels
  python3 picard_qcsum.py --list-panels

  # Merge multiple qcsum files
  python3 picard_qcsum.py --merge sample1.qcsum.txt sample2.qcsum.txt -o /path/to/output

  # Batch process multiple CRAMs and merge results
  python3 picard_qcsum.py --batch sample1.cram sample2.cram --panel pcp_tumor -o /path/to/output
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import logging
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
from typing import Iterable, Optional
import yaml

SCRIPT_DIR = Path(__file__).parent.resolve()
QCSUM_CONFIG_YAML = SCRIPT_DIR / "config" / "qcsum_config.yaml"
PICARD_REF = "/mnt/pns/tracks/ref/hg38.fa"

PANEL_TO_SAMPLE_TYPE = {
    "pcp_tumor": "Tumor",
    "pcp_normal": "Normal",
}

QC_SUM_COLUMNS = (
    "Sample",
    "Tumor_Normal",
    "Sequencing_Platform",
    "Pipeline_version",
    "Alignment_QC",
    "Coverage_QC",
    "Total_Reads",
    "%Reads_Aligned",
    "Capture",
    "Avg_Capture_Coverage",
    "%On/Near_Bait_Bases",
    "%On_Bait_Bases",
    "%On_Target_Bases",
    "FOLD_80_BASE_PENALTY",
    "AT_DROPOUT",
    "GC_DROPOUT",
    "Avg_ROI_Coverage",
    "MEDIAN_ROI_COVERAGE",
    "MAX_ROI_COVERAGE",
    "%ROI_1x",
    "%ROI_20x",
    "%ROI_100x",
    "%ROI_250x",
    "%ROI_500x",
)
QC_SUM_HEADER = ",".join(QC_SUM_COLUMNS)

logger = logging.getLogger(__name__)


def _load_raw_config(config_path: Path | str = QCSUM_CONFIG_YAML) -> dict:
    """Load and validate the raw panel YAML configuration."""
    path = Path(config_path)
    with open(path) as f:
        yaml_obj = yaml.safe_load(f)
    if not isinstance(yaml_obj, dict) or not yaml_obj:
        raise ValueError(f"Failed to load {path}. Ensure it is a valid non-empty YAML mapping.")
    return yaml_obj


class QCSumConfig:
    """Manages QC summary configuration for a given panel."""

    def __init__(self, panel: str, config_path: Path | str = QCSUM_CONFIG_YAML):
        self.panel = panel
        self.config_path = Path(config_path)
        self.config = self._load_config()

    @classmethod
    def available_panels(cls, config_path: Path | str = QCSUM_CONFIG_YAML) -> list[str]:
        """Return the list of configured panel names from YAML."""
        return list(_load_raw_config(config_path).keys())

    def _load_config(self) -> dict[str, str]:
        """Load and validate panel configuration from YAML."""
        yaml_obj = _load_raw_config(self.config_path)

        if self.panel not in yaml_obj:
            available = ", ".join(yaml_obj.keys())
            raise ValueError(
                f"Panel '{self.panel}' not found in {self.config_path}. Available panels: {available}"
            )

        # Heme subpanels inherit defaults from heme_comp and override any specified keys
        if self.panel.startswith("heme_") and self.panel != "heme_comp":
            base = yaml_obj.get("heme_comp") or {}
            yaml_dict = dict(base)
            subpanel_dict = yaml_obj.get(self.panel) or {}
            yaml_dict.update({k: v for k, v in subpanel_dict.items() if v is not None})
        else:
            yaml_dict = dict(yaml_obj.get(self.panel) or {})

        return {k: str(v) for k, v in yaml_dict.items()}

    def get(self, key: str, default: str = "") -> str:
        """Get a config value by key."""
        return self.config.get(key, default)


@dataclass(frozen=True)
class QCSumSummary:
    """Evaluated QC summary metrics and PASS/FAIL status for a single sample."""

    sample_id: str
    sample_type: str
    alignment_qc: str
    coverage_qc: str
    total_reads: int
    pct_aligned: float
    mean_target_coverage: float
    covered_threshold: int
    pct_roi_at_threshold: float
    csv_line: str


def _write_text_atomic(path: Path, content: str) -> None:
    """Atomically write text content to path via a temporary file in the same directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as tmp_file:
            tmp_file.write(content)
        os.replace(tmp_name, path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise


def _parse_float(metrics: dict, key: str, default: float = 0.0) -> float:
    """Parse a float metric from Picard output, normalizing missing/NaN values to default."""
    raw = metrics.get(key)
    if raw is None or raw in ("", "?"):
        return default
    val = float(raw)
    return default if math.isnan(val) else val


def run_command(cmd: list[str], dry_run: bool = False) -> subprocess.CompletedProcess:
    """Execute a command with logging and error handling."""
    cmd_str = shlex.join(cmd)
    logger.info("Running: %s", cmd_str)

    if dry_run:
        logger.info("[DRY RUN] Command not executed")
        return subprocess.CompletedProcess(cmd, 0, b"", b"")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.stderr:
            logger.warning("stderr: %s", result.stderr.decode(errors="replace"))
        return result
    except subprocess.CalledProcessError as e:
        logger.error("Command failed (exit %s): %s", e.returncode, cmd_str)
        if e.stderr:
            logger.error("stderr: %s", e.stderr.decode(errors="replace"))
        raise


def run_picard_hs_metrics(
    cram_file: str,
    output_file: str,
    bait_intervals: str,
    target_intervals: str,
    reference: str = PICARD_REF,
    dry_run: bool = False,
) -> None:
    """Run Picard CollectHsMetrics on a CRAM/BAM file."""
    picard_cmd = [
        "picard",
        "CollectHsMetrics",
        f"I={cram_file}",
        f"O={output_file}",
        f"R={reference}",
        f"BAIT_INTERVALS={bait_intervals}",
        f"TARGET_INTERVALS={target_intervals}",
    ]
    run_command(picard_cmd, dry_run)


def parse_hs_metrics(hsm_file: str | Path) -> dict[str, str]:
    """
    Parse Picard HsMetrics output file.

    Locates the METRICS header row starting with ``BAIT_SET`` and maps each
    column name to its value from the first (aggregate) metrics data row.
    Missing values (``?`` or empty) are normalized to ``"0"``.
    """
    with open(hsm_file) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith("BAIT_SET\t"):
            if i + 1 >= len(lines) or not lines[i + 1].strip():
                raise ValueError(f"Missing metrics data row in {hsm_file}")
            headers = line.rstrip("\r\n").split("\t")
            values = lines[i + 1].rstrip("\r\n").split("\t")
            normalized = ["0" if x in ("?", "") else x for x in values]
            return dict(zip(headers, normalized))

    raise ValueError(f"Could not parse metrics from {hsm_file}")


def evaluate_qcsum(
    sample_id: str,
    metrics: dict[str, str],
    config: QCSumConfig,
    sample_type: Optional[str] = None,
) -> QCSumSummary:
    """
    Evaluate alignment/coverage QC PASS/FAIL status and build a CSV line matching QC_SUM_HEADER.

    Args:
        sample_id: Sample identifier.
        metrics: Parsed Picard HsMetrics dictionary.
        config: Loaded QCSumConfig for the sample's panel.
        sample_type: Optional override for the Tumor_Normal column. When None,
            defaults to "Tumor"/"Normal" for PCP panels and "" for Heme panels.
    """
    total_reads_float = _parse_float(metrics, "TOTAL_READS")
    total_reads = int(total_reads_float)
    mean_target_coverage = _parse_float(metrics, "MEAN_TARGET_COVERAGE")
    pct_aligned = _parse_float(metrics, "PCT_PF_UQ_READS_ALIGNED") * 100

    covered_threshold = int(config.get("covered", "100"))
    roi_coverage_col = f"PCT_TARGET_BASES_{covered_threshold}X"
    if roi_coverage_col not in metrics:
        raise ValueError(
            f"Configured coverage column '{roi_coverage_col}' not found in Picard HsMetrics"
        )
    pct_roi_at_threshold = _parse_float(metrics, roi_coverage_col) * 100

    # Determine alignment QC (PASS/FAIL)
    fail_min_align = float(config.get("fail_min_align_pct", "0"))
    fail_min_reads = float(config.get("fail_min_reads", "0"))
    alignment_qc = "FAIL" if (pct_aligned < fail_min_align or total_reads_float < fail_min_reads) else "PASS"

    # Determine coverage QC (PASS/FAIL)
    fail_min_roi = float(config.get("fail_min_roi_pct", "0"))
    fail_min_avgcov = float(config.get("fail_min_avgcov", "0"))
    coverage_qc = "FAIL" if (pct_roi_at_threshold < fail_min_roi or mean_target_coverage < fail_min_avgcov) else "PASS"

    def get_pct(key: str) -> str:
        """Convert Picard fraction to percentage string rounded to 2 decimals."""
        return f"{_parse_float(metrics, key) * 100:.2f}"

    # ON_BAIT_BASES / PF_BASES_ALIGNED (and matching ON_TARGET_BASES / PF_BASES_ALIGNED):
    # excludes near-bait bases without penalizing PCR duplicates like PCT_USABLE_BASES_*
    pf_bases_aligned = _parse_float(metrics, "PF_BASES_ALIGNED")
    pct_on_bait = (_parse_float(metrics, "ON_BAIT_BASES") / pf_bases_aligned * 100) if pf_bases_aligned else 0.0
    pct_on_target = (_parse_float(metrics, "ON_TARGET_BASES") / pf_bases_aligned * 100) if pf_bases_aligned else 0.0

    resolved_sample_type = (
        PANEL_TO_SAMPLE_TYPE.get(config.panel, "") if sample_type is None else sample_type
    )

    fields = [
        sample_id,
        resolved_sample_type,
        config.get("platform", "NovaSeq6000"),
        config.get("pipeline_version", "Unknown"),
        alignment_qc,
        coverage_qc,
        str(metrics.get("TOTAL_READS", "0")),
        f"{pct_aligned:.2f}",
        config.get("capture", "unknown"),
        f"{_parse_float(metrics, 'MEAN_BAIT_COVERAGE'):.2f}",
        get_pct("PCT_SELECTED_BASES"),  # % on/near bait bases
        f"{pct_on_bait:.2f}",           # % on bait bases
        f"{pct_on_target:.2f}",         # % on target bases
        f"{_parse_float(metrics, 'FOLD_80_BASE_PENALTY'):.2f}",
        f"{_parse_float(metrics, 'AT_DROPOUT'):.2f}",
        f"{_parse_float(metrics, 'GC_DROPOUT'):.2f}",
        f"{mean_target_coverage:.2f}",
        str(metrics.get("MEDIAN_TARGET_COVERAGE", "0")),
        str(metrics.get("MAX_TARGET_COVERAGE", "0")),
        get_pct("PCT_TARGET_BASES_1X"),
        get_pct("PCT_TARGET_BASES_20X"),
        get_pct("PCT_TARGET_BASES_100X"),
        get_pct("PCT_TARGET_BASES_250X"),
        get_pct("PCT_TARGET_BASES_500X"),
    ]
    assert len(fields) == len(QC_SUM_COLUMNS), (
        f"Field count ({len(fields)}) does not match QC_SUM_COLUMNS ({len(QC_SUM_COLUMNS)})"
    )

    return QCSumSummary(
        sample_id=sample_id,
        sample_type=resolved_sample_type,
        alignment_qc=alignment_qc,
        coverage_qc=coverage_qc,
        total_reads=total_reads,
        pct_aligned=pct_aligned,
        mean_target_coverage=mean_target_coverage,
        covered_threshold=covered_threshold,
        pct_roi_at_threshold=pct_roi_at_threshold,
        csv_line=",".join(fields),
    )


def generate_qcsum_line(
    sample_id: str,
    metrics: dict[str, str],
    config: QCSumConfig,
    sample_type: Optional[str] = None,
) -> str:
    """Generate a single CSV line for qcsum output matching QC_SUM_HEADER."""
    return evaluate_qcsum(sample_id, metrics, config, sample_type=sample_type).csv_line


def derive_sample_id(cram_file: str | Path) -> str:
    """Derive a sample identifier from a CRAM/BAM path."""
    return Path(cram_file).stem


def run_qcsum(
    cram_file: str,
    panel: str,
    sample_id: Optional[str] = None,
    output_dir: Optional[str] = None,
    reference: str = PICARD_REF,
    dry_run: bool = False,
    sample_type: Optional[str] = None,
    config_path: Path | str = QCSUM_CONFIG_YAML,
) -> dict:
    """
    Run the complete Picard CollectHsMetrics + qcsum workflow on a single sample.

    Args:
        cram_file: Path to input CRAM/BAM file.
        panel: Panel name for selecting cutoffs (e.g., 'pcp_tumor', 'heme_comp').
        sample_id: Sample identifier. If not provided, derived from filename.
        output_dir: Directory for output files. Defaults to same directory as CRAM.
        reference: Path to reference FASTA for Picard.
        dry_run: If True, print commands without executing.
        sample_type: Optional override for Tumor_Normal column.
        config_path: Path to qcsum_config.yaml.

    Returns:
        dict with QC results and output file paths.
    """
    cram_path = Path(cram_file).resolve()

    if sample_id is None:
        sample_id = derive_sample_id(cram_path)

    if output_dir is None:
        out_path = cram_path.parent
    else:
        out_path = Path(output_dir).resolve()
        if not dry_run:
            out_path.mkdir(parents=True, exist_ok=True)

    logger.info("Processing sample: %s", sample_id)
    logger.info("Input CRAM: %s", cram_path)
    logger.info("Panel: %s", panel)
    logger.info("Output directory: %s", out_path)

    config = QCSumConfig(panel, config_path=config_path)

    bait_intervals = config.get("bait_intervals")
    target_intervals = config.get("target_intervals")
    if not bait_intervals or not target_intervals:
        raise ValueError(f"Panel '{panel}' missing bait_intervals or target_intervals")

    hsm_file = out_path / f"{sample_id}.hsm.txt"
    qcsum_file = out_path / f"{sample_id}.qcsum.txt"

    # Run Picard CollectHsMetrics (skip if a non-empty output already exists)
    if hsm_file.is_file() and hsm_file.stat().st_size > 0:
        logger.info("Picard output exists, skipping: %s", hsm_file)
    else:
        if not cram_path.exists() and not dry_run:
            raise FileNotFoundError(f"CRAM/BAM file not found: {cram_path}")
        run_picard_hs_metrics(
            cram_file=str(cram_path),
            output_file=str(hsm_file),
            bait_intervals=bait_intervals,
            target_intervals=target_intervals,
            reference=reference,
            dry_run=dry_run,
        )

    if dry_run:
        logger.info("[DRY RUN] Skipping metrics parsing")
        return {"sample_id": sample_id, "dry_run": True}

    logger.info("Parsing metrics from: %s", hsm_file)
    metrics = parse_hs_metrics(hsm_file)
    summary = evaluate_qcsum(sample_id, metrics, config, sample_type=sample_type)

    logger.info("Alignment QC: %s", summary.alignment_qc)
    logger.info("  Total reads: %s", f"{summary.total_reads:,}")
    logger.info("  %% Aligned: %.2f%%", summary.pct_aligned)
    logger.info("Coverage QC: %s", summary.coverage_qc)
    logger.info("  Mean target coverage: %.2fx", summary.mean_target_coverage)
    logger.info("  %% ROI at %dx: %.2f%%", summary.covered_threshold, summary.pct_roi_at_threshold)

    _write_text_atomic(qcsum_file, f"{QC_SUM_HEADER}\n{summary.csv_line}\n")
    logger.info("QC summary written to: %s", qcsum_file)

    return {
        "sample_id": sample_id,
        "panel": panel,
        "sample_type": summary.sample_type,
        "hsm_file": str(hsm_file),
        "qcsum_file": str(qcsum_file),
        "alignment_qc": summary.alignment_qc,
        "coverage_qc": summary.coverage_qc,
    }


def merge_qcsum_files(
    qcsum_files: Iterable[str | Path],
    output_dir: Optional[str | Path] = None,
    output_name: str = "merged_qcsum.csv",
) -> str:
    """
    Merge multiple qcsum files into a single CSV file.

    Args:
        qcsum_files: Paths to qcsum files to merge.
        output_dir: Directory for output file. Defaults to current directory.
        output_name: Name of the merged output file.

    Returns:
        Path to the merged output file.
    """
    qcsum_file_list = list(qcsum_files)
    if not qcsum_file_list:
        raise ValueError("No qcsum files provided to merge")

    out_dir = Path.cwd() if output_dir is None else Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    output_file = out_dir / output_name

    all_lines: list[str] = []
    header: Optional[str] = None
    files_merged = 0

    for qcsum_file in qcsum_file_list:
        qcsum_path = Path(qcsum_file)
        if not qcsum_path.exists():
            logger.warning("Skipping missing file: %s", qcsum_path)
            continue

        with open(qcsum_path) as f:
            lines = f.readlines()

        if not lines:
            logger.warning("Skipping empty file: %s", qcsum_path)
            continue

        file_header = lines[0].strip()
        if header is None:
            header = file_header
        elif file_header != header:
            raise ValueError(
                f"Header mismatch in {qcsum_path}: expected '{header}', got '{file_header}'"
            )

        files_merged += 1
        for line in lines[1:]:
            stripped = line.strip()
            if stripped:
                all_lines.append(stripped)

    if header is None:
        raise ValueError("No valid qcsum files found to merge")

    content = header + "\n" + "".join(f"{line}\n" for line in all_lines)
    _write_text_atomic(output_file, content)

    logger.info("Merged %d samples from %d files", len(all_lines), files_merged)
    logger.info("Output: %s", output_file)

    return str(output_file)


def list_panels(config_path: Path | str = QCSUM_CONFIG_YAML) -> list[str]:
    """Print and return available panels from the config file."""
    yaml_obj = _load_raw_config(config_path)

    print("Available panels:")
    for panel, config in yaml_obj.items():
        if config and "pipeline_version" in config:
            print(f"  {panel:15s} - {config.get('pipeline_version', '')}")
            print(f"    Coverage threshold: {config.get('covered', 'N/A')}x")
            print(f"    Min avg coverage: {config.get('pass_min_avgcov', 'N/A')}x")
        else:
            print(f"  {panel:15s} - (subpanel, inherits from heme_comp)")
    return list(yaml_obj.keys())


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command line arguments."""
    panel_choices = (
        QCSumConfig.available_panels() if QCSUM_CONFIG_YAML.exists() else None
    )

    parser = argparse.ArgumentParser(
        description="Run Picard CollectHsMetrics and QC summary on a CRAM file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on a PCP tumor sample
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel pcp_tumor

  # Run on a Heme sample with custom sample ID
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel heme_comp --sample-id MySample

  # Run on a Heme subpanel (MPN)
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel heme_mpn

  # Output to a specific directory
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel pcp_normal --output-dir /path/to/output

  # List available panels
  python3 picard_qcsum.py --list-panels

  # Merge multiple qcsum files
  python3 picard_qcsum.py --merge sample1.qcsum.txt sample2.qcsum.txt -o /path/to/output

  # Batch process multiple CRAMs and merge results
  python3 picard_qcsum.py --batch sample1.cram sample2.cram --panel pcp_tumor -o /path/to/output
        """,
    )

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--cram", "-c", help="Path to input CRAM/BAM file")
    mode_group.add_argument("--batch", "-b", nargs="+", metavar="CRAM", help="Batch process multiple CRAMs")
    mode_group.add_argument("--merge", "-m", nargs="+", metavar="FILE", help="Merge qcsum files")
    mode_group.add_argument("--list-panels", action="store_true", help="List available panels")

    parser.add_argument("--panel", "-p", choices=panel_choices, help="Panel for QC cutoffs")
    parser.add_argument("--sample-id", "-s", help="Sample identifier (default: from filename)")
    parser.add_argument("--output-dir", "-o", help="Output directory")
    parser.add_argument("--reference", "-r", default=PICARD_REF, help="Reference FASTA")
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=None,
        help="Max parallel Picard jobs in batch mode (default: min(num_crams, cpu_count))",
    )
    parser.add_argument("--dry-run", "-n", action="store_true", help="Print commands only")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument("--merge-output", default="merged_qcsum.csv", help="Merged output filename")

    return parser.parse_args(argv)


def _check_batch_output_collisions(cram_files: list[str], output_dir: Optional[str]) -> None:
    """Raise ValueError if two batch inputs would write to the same output qcsum file."""
    seen: dict[Path, str] = {}
    for cram_file in cram_files:
        cram_path = Path(cram_file).resolve()
        sample_id = derive_sample_id(cram_path)
        out_dir = Path(output_dir).resolve() if output_dir else cram_path.parent
        target = out_dir / f"{sample_id}.qcsum.txt"
        if target in seen:
            raise ValueError(
                f"Batch output collision for sample '{sample_id}' ({target}): "
                f"'{seen[target]}' and '{cram_file}'"
            )
        seen[target] = cram_file


def main(argv: Optional[list[str]] = None) -> int:
    """Main entry point."""
    args = parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    if args.threads is not None and args.threads < 1:
        logger.error("--threads must be >= 1, got %d", args.threads)
        return 1

    if args.list_panels:
        list_panels()
        return 0

    if args.merge:
        try:
            output_file = merge_qcsum_files(
                qcsum_files=args.merge,
                output_dir=args.output_dir,
                output_name=args.merge_output,
            )
            print(f"Merged qcsum written to: {output_file}")
            return 0
        except ValueError as e:
            logger.error("%s", e)
            return 1

    if args.batch:
        if not args.panel:
            print("Error: --panel is required for batch mode", file=sys.stderr)
            return 1
        if args.sample_id:
            print("Error: --sample-id cannot be used with --batch", file=sys.stderr)
            return 1
        try:
            _check_batch_output_collisions(args.batch, args.output_dir)
        except ValueError as e:
            logger.error("%s", e)
            return 1

        results_by_cram: dict[str, str] = {}
        failed: list[str] = []

        max_workers = (
            args.threads if args.threads is not None else min(len(args.batch), os.cpu_count() or 1)
        )
        logger.info(
            "Batch mode: processing %d CRAMs with %d parallel workers",
            len(args.batch),
            max_workers,
        )

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_cram = {
                executor.submit(
                    run_qcsum,
                    cram_file=cram_file,
                    panel=args.panel,
                    output_dir=args.output_dir,
                    reference=args.reference,
                    dry_run=args.dry_run,
                ): cram_file
                for cram_file in args.batch
            }

            for future in as_completed(future_to_cram):
                cram_file = future_to_cram[future]
                try:
                    result = future.result()
                    if not args.dry_run:
                        results_by_cram[cram_file] = result["qcsum_file"]
                        logger.info(
                            "Processed: %s - Align:%s Cov:%s",
                            result["sample_id"],
                            result["alignment_qc"],
                            result["coverage_qc"],
                        )
                except Exception as e:
                    logger.error("Failed to process %s: %s", cram_file, e)
                    failed.append(cram_file)

        if args.dry_run:
            return 0

        # Preserve input order when merging batch outputs
        qcsum_files = [results_by_cram[c] for c in args.batch if c in results_by_cram]
        if qcsum_files:
            merged_file = merge_qcsum_files(
                qcsum_files=qcsum_files,
                output_dir=args.output_dir,
                output_name=args.merge_output,
            )
            print(f"\nProcessed {len(qcsum_files)} samples, {len(failed)} failed")
            print(f"Merged qcsum: {merged_file}")

        return 1 if failed else 0

    if not args.cram:
        print(
            "Error: --cram, --batch, --merge, or --list-panels is required",
            file=sys.stderr,
        )
        return 1

    if not args.panel:
        print("Error: --panel is required (use --list-panels to see available panels)", file=sys.stderr)
        return 1

    try:
        result = run_qcsum(
            cram_file=args.cram,
            panel=args.panel,
            sample_id=args.sample_id,
            output_dir=args.output_dir,
            reference=args.reference,
            dry_run=args.dry_run,
        )

        if not args.dry_run:
            print("\n" + "=" * 60)
            print(f"Sample:       {result['sample_id']}")
            print(f"Panel:        {result['panel']}")
            print(f"Alignment QC: {result['alignment_qc']}")
            print(f"Coverage QC:  {result['coverage_qc']}")
            print(f"HSM file:     {result['hsm_file']}")
            print(f"QCSum file:   {result['qcsum_file']}")
            print("=" * 60)

        return 0

    except (FileNotFoundError, ValueError) as e:
        logger.error("%s", e)
        return 1
    except subprocess.CalledProcessError as e:
        logger.error("Command failed with exit code %s", e.returncode)
        return e.returncode


if __name__ == "__main__":
    sys.exit(main())
