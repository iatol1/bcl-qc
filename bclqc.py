#!/usr/bin/env python3
"""
Unified CLI and library for BCL demux, DRAGEN alignment, Picard CollectHsMetrics
QC summaries, and MultiQC reporting.

Subcommands:
  bclqc.py run    --run-dir RUN [--sampleinfo FILE] [--steps demux align qc]
  bclqc.py demux  --run-dir RUN [--fastqs-dir DIR]
  bclqc.py align  --run-dir RUN [--fastqs-dir DIR] [--bams-dir DIR]
  bclqc.py qc     --run-dir RUN --sampleinfo FILE [--fastqs-dir DIR] [--bams-dir DIR]
  bclqc.py qcsum  CRAM [CRAM ...] --panel PANEL [-o DIR] [-s SAMPLE_ID]
  bclqc.py merge  FILE [FILE ...] [-o DIR] [--output NAME]
  bclqc.py panels
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import shlex
import subprocess
import sys
import tempfile
from collections import deque
from collections.abc import Callable, Iterable
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).parent.resolve()
PANELS_YAML = SCRIPT_DIR / "config" / "panels.yaml"
MULTIQC_CONFIG_YAML = SCRIPT_DIR / "config" / "multiqc_config.yaml"

DEF_STEPS = ("demux", "align", "qc")
PICARD_XMX_RE = re.compile(r"^(-Xmx)?([1-9][0-9]*[kKmMgG])$")

QC_SUM_COLUMNS = (
    "Sample", "Tumor_Normal", "Sequencing_Platform", "Pipeline_version",
    "Alignment_QC", "Coverage_QC", "Total_Reads", "%Reads_Aligned", "Capture",
    "Avg_Capture_Coverage", "%On/Near_Bait_Bases", "%On_Bait_Bases",
    "%On_Target_Bases", "FOLD_80_BASE_PENALTY", "AT_DROPOUT", "GC_DROPOUT",
    "Avg_ROI_Coverage", "MEDIAN_ROI_COVERAGE", "MAX_ROI_COVERAGE",
    "%ROI_1x", "%ROI_20x", "%ROI_100x", "%ROI_250x", "%ROI_500x",
)
QC_SUM_HEADER = ",".join(QC_SUM_COLUMNS)

LANE_COLORS = (
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
)

logger = logging.getLogger(__name__)


class PanelsConfig:
    """Declarative pipeline and QC configuration loaded from config/panels.yaml."""

    def __init__(self, config_path: Path | str = PANELS_YAML):
        self.config_path = Path(config_path)
        with self.config_path.open() as f:
            raw = yaml.safe_load(f)
        if not isinstance(raw, dict) or not raw:
            raise ValueError(f"Failed to load {self.config_path}: expected non-empty YAML mapping")

        self.defaults: dict[str, str] = {k: str(v) for k, v in (raw.get("defaults") or {}).items()}
        self.align_groups: dict[str, dict] = dict(raw.get("align_groups") or {})
        self._raw_qc_panels: dict[str, dict] = dict(raw.get("qc_panels") or {})

        self.qc_panels: dict[str, dict[str, str]] = {}
        self.aliases: dict[str, str] = {}
        for name, raw_entry in self._raw_qc_panels.items():
            self.qc_panels[name] = self._resolve_qc_panel(name, stack=())
            aliases = (raw_entry or {}).get("aliases")
            if aliases is None:
                continue
            if not isinstance(aliases, list):
                raise ValueError(f"Panel '{name}' aliases must be a list, got {type(aliases).__name__}")
            for alias in aliases:
                alias_str = str(alias)
                if alias_str in self.aliases:
                    raise ValueError(
                        f"Duplicate alias '{alias_str}' for panels '{self.aliases[alias_str]}' and '{name}'"
                    )
                self.aliases[alias_str] = name

    def _resolve_qc_panel(self, name: str, stack: tuple[str, ...]) -> dict[str, str]:
        if name in stack:
            raise ValueError(f"Circular panel inheritance in {self.config_path}: {' -> '.join([*stack, name])}")
        if name not in self._raw_qc_panels:
            raise ValueError(f"Panel '{name}' not found in {self.config_path}. Available: {', '.join(self._raw_qc_panels)}")

        entry = dict(self._raw_qc_panels[name] or {})
        parent_name = entry.pop("inherits", None)
        entry.pop("aliases", None)

        merged: dict[str, object] = {}
        if "platform" in self.defaults:
            merged["platform"] = self.defaults["platform"]
        if parent_name:
            merged.update(self._resolve_qc_panel(str(parent_name), stack=(*stack, name)))
        merged.update({k: v for k, v in entry.items() if v is not None})
        return {k: str(v) for k, v in merged.items()}

    def resolve_panel_name(self, panel_or_alias: str) -> str:
        """Resolve a sampleinfo display name or canonical key to a configured QC panel."""
        canonical = self.aliases.get(panel_or_alias, panel_or_alias)
        if canonical not in self.qc_panels:
            raise ValueError(f"Unknown panel '{panel_or_alias}'. Available: {', '.join(self.qc_panels)}")
        return canonical

    def get_qc_panel(self, panel_or_alias: str) -> tuple[str, dict[str, str]]:
        """Return (canonical_panel_name, resolved_settings)."""
        name = self.resolve_panel_name(panel_or_alias)
        return name, self.qc_panels[name]

    def get_align_group(self, fastq_list: str | Path, fastqs_dir: str | Path | None = None) -> tuple[str, dict]:
        """Match a fastq_list.csv path to its align_groups entry using the innermost index directory."""
        fq_path = Path(fastq_list)
        if fastqs_dir is not None:
            try:
                rel_parts = fq_path.relative_to(Path(fastqs_dir)).parts
                if rel_parts and rel_parts[0] in self.align_groups:
                    return rel_parts[0], self.align_groups[rel_parts[0]]
            except ValueError:
                pass
        for part in reversed(fq_path.parent.parts):
            if part in self.align_groups:
                return part, self.align_groups[part]
        raise ValueError(f"Could not determine alignment group for fastq list: {fastq_list}")


@dataclass(frozen=True)
class QCSumSummary:
    """Evaluated QC summary metrics and PASS/FAIL status for a single sample."""

    sample_id: str
    panel: str
    sample_type: str
    alignment_qc: str
    coverage_qc: str
    total_reads: int
    pct_aligned: float
    mean_target_coverage: float
    covered_threshold: int
    pct_roi_at_threshold: float
    csv_line: str


def write_text_atomic(path: Path, content: str) -> None:
    """Atomically write text content to path via a temp file in the same directory."""
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


def normalize_picard_xmx(xmx: str) -> str:
    """Validate and normalize a JVM heap flag such as '4g' or '-Xmx4g' into '-Xmx4g'."""
    match = PICARD_XMX_RE.fullmatch(xmx.strip())
    if not match:
        raise ValueError(f"Invalid --picard-xmx value {xmx!r}; expected a JVM heap size like '4g' or '-Xmx8g'")
    return f"-Xmx{match.group(2)}"


def resolve_worker_count(num_tasks: int, requested: int | None, config: PanelsConfig) -> int:
    """Return a bounded concurrent Picard worker count (default capped by panels.yaml qc_workers)."""
    if requested is not None:
        if requested < 1:
            raise ValueError(f"--workers must be >= 1, got {requested}")
        return max(1, min(num_tasks, requested))
    default_cap = int(config.defaults.get("qc_workers", "4"))
    return max(1, min(num_tasks, default_cap, os.cpu_count() or 1))


def run_cmd(
    cmd: list[str],
    dry_run: bool = False,
    log_file: Path | str | None = None,
) -> subprocess.CompletedProcess:
    """Execute a subprocess command, optionally streaming stdout/stderr directly to log_file."""
    cmd_str = shlex.join(cmd)
    logger.info("Running: %s", cmd_str)

    if dry_run:
        logger.info("[DRY RUN] Command not executed")
        return subprocess.CompletedProcess(cmd, 0, b"", b"")

    if log_file is not None:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with log_path.open("w") as lf:
                return subprocess.run(cmd, check=True, stdout=lf, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            logger.error("Command failed (exit %s): %s (log: %s)", e.returncode, cmd_str, log_path)
            try:
                with log_path.open(errors="replace") as lf:
                    tail = "".join(deque(lf, maxlen=20)).strip()
                if tail:
                    logger.error("Last lines of %s:\n%s", log_path, tail)
            except OSError:
                pass
            raise

    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.stderr:
            logger.warning("stderr: %s", result.stderr.decode(errors="replace"))
        return result
    except subprocess.CalledProcessError as e:
        logger.error("Command failed (exit %s): %s", e.returncode, cmd_str)
        if e.stdout:
            logger.error("stdout: %s", e.stdout.decode(errors="replace"))
        if e.stderr:
            logger.error("stderr: %s", e.stderr.decode(errors="replace"))
        raise


def _parse_float(metrics: dict[str, str], key: str, default: float = 0.0) -> float:
    raw = metrics.get(key)
    if raw is None or raw in ("", "?"):
        return default
    val = float(raw)
    return default if math.isnan(val) else val


def parse_hs_metrics(hsm_file: str | Path) -> dict[str, str]:
    """Parse Picard CollectHsMetrics output into a {column: value} mapping."""
    with open(hsm_file) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith("BAIT_SET\t"):
            if i + 1 >= len(lines) or not lines[i + 1].strip():
                raise ValueError(f"Missing metrics data row in {hsm_file}")
            headers = line.rstrip("\r\n").split("\t")
            values = lines[i + 1].rstrip("\r\n").split("\t")
            return dict(zip(headers, ["0" if x in ("?", "") else x for x in values]))

    raise ValueError(f"Could not parse metrics from {hsm_file}")


def evaluate_qcsum(
    sample_id: str,
    panel: str,
    metrics: dict[str, str],
    panel_cfg: dict[str, str],
    sample_type: str | None = None,
) -> QCSumSummary:
    """Evaluate alignment/coverage QC PASS/FAIL status and format the CSV row."""
    total_reads_float = _parse_float(metrics, "TOTAL_READS")
    mean_target_coverage = _parse_float(metrics, "MEAN_TARGET_COVERAGE")
    pct_aligned = _parse_float(metrics, "PCT_PF_UQ_READS_ALIGNED") * 100

    covered_threshold = int(panel_cfg.get("covered", "100"))
    roi_col = f"PCT_TARGET_BASES_{covered_threshold}X"
    if roi_col not in metrics:
        raise ValueError(f"Configured coverage column '{roi_col}' not found in Picard HsMetrics")
    pct_roi = _parse_float(metrics, roi_col) * 100

    alignment_qc = (
        "FAIL"
        if pct_aligned < float(panel_cfg.get("fail_min_align_pct", "0"))
        or total_reads_float < float(panel_cfg.get("fail_min_reads", "0"))
        else "PASS"
    )
    coverage_qc = (
        "FAIL"
        if pct_roi < float(panel_cfg.get("fail_min_roi_pct", "0"))
        or mean_target_coverage < float(panel_cfg.get("fail_min_avgcov", "0"))
        else "PASS"
    )

    def get_pct(key: str) -> str:
        return f"{_parse_float(metrics, key) * 100:.2f}"

    # ON_BAIT_BASES / PF_BASES_ALIGNED (and matching ON_TARGET_BASES / PF_BASES_ALIGNED)
    # excludes near-bait bases without penalizing PCR duplicates like PCT_USABLE_BASES_*.
    pf_bases_aligned = _parse_float(metrics, "PF_BASES_ALIGNED")
    pct_on_bait = (_parse_float(metrics, "ON_BAIT_BASES") / pf_bases_aligned * 100) if pf_bases_aligned else 0.0
    pct_on_target = (_parse_float(metrics, "ON_TARGET_BASES") / pf_bases_aligned * 100) if pf_bases_aligned else 0.0
    resolved_sample_type = panel_cfg.get("sample_type", "") if sample_type is None else sample_type

    fields = [
        sample_id,
        resolved_sample_type,
        panel_cfg.get("platform", "NovaSeq6000"),
        panel_cfg.get("pipeline_version", "Unknown"),
        alignment_qc,
        coverage_qc,
        str(metrics.get("TOTAL_READS", "0")),
        f"{pct_aligned:.2f}",
        panel_cfg.get("capture", "unknown"),
        f"{_parse_float(metrics, 'MEAN_BAIT_COVERAGE'):.2f}",
        get_pct("PCT_SELECTED_BASES"),
        f"{pct_on_bait:.2f}",
        f"{pct_on_target:.2f}",
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
    assert len(fields) == len(QC_SUM_COLUMNS)

    return QCSumSummary(
        sample_id=sample_id,
        panel=panel,
        sample_type=resolved_sample_type,
        alignment_qc=alignment_qc,
        coverage_qc=coverage_qc,
        total_reads=int(total_reads_float),
        pct_aligned=pct_aligned,
        mean_target_coverage=mean_target_coverage,
        covered_threshold=covered_threshold,
        pct_roi_at_threshold=pct_roi,
        csv_line=",".join(fields),
    )


def run_picard_hs_metrics(
    cram_file: str | Path,
    output_file: str | Path,
    bait_intervals: str,
    target_intervals: str,
    reference: str,
    picard_xmx: str = "4g",
    log_file: Path | str | None = None,
    dry_run: bool = False,
) -> None:
    """Run Picard CollectHsMetrics with bounded JVM heap and GC thread count."""
    run_cmd(
        [
            "picard",
            normalize_picard_xmx(picard_xmx),
            "-XX:ActiveProcessorCount=2",
            "CollectHsMetrics",
            f"I={cram_file}",
            f"O={output_file}",
            f"R={reference}",
            f"BAIT_INTERVALS={bait_intervals}",
            f"TARGET_INTERVALS={target_intervals}",
        ],
        dry_run=dry_run,
        log_file=log_file,
    )


def run_qcsum(
    cram_file: str | Path,
    panel: str,
    config: PanelsConfig,
    sample_id: str | None = None,
    output_dir: str | Path | None = None,
    reference: str | None = None,
    picard_xmx: str | None = None,
    sample_type: str | None = None,
    dry_run: bool = False,
) -> tuple[Path, QCSumSummary | None]:
    """Run Picard CollectHsMetrics (skipping if non-empty .hsm.txt exists) and write .qcsum.txt."""
    if not str(cram_file).strip():
        raise ValueError("CRAM/BAM path must not be empty")

    panel_name, panel_cfg = config.get_qc_panel(panel)
    bait_intervals = panel_cfg.get("bait_intervals")
    target_intervals = panel_cfg.get("target_intervals")
    if not bait_intervals or not target_intervals:
        raise ValueError(f"Panel '{panel_name}' missing bait_intervals or target_intervals")

    picard_ref = reference or config.defaults.get("picard_ref")
    if not picard_ref:
        raise ValueError("No Picard reference FASTA configured (pass --picard-ref or set defaults.picard_ref)")
    xmx = picard_xmx or config.defaults.get("picard_xmx", "4g")

    cram_path = Path(cram_file).resolve()
    sample_id = sample_id or cram_path.stem
    out_dir = cram_path.parent if output_dir is None else Path(output_dir).resolve()
    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    hsm_file = out_dir / f"{sample_id}.hsm.txt"
    qcsum_file = out_dir / f"{sample_id}.qcsum.txt"

    if hsm_file.is_file() and hsm_file.stat().st_size > 0:
        logger.info("Picard output exists, skipping: %s", hsm_file)
    else:
        if not dry_run and not cram_path.is_file():
            raise FileNotFoundError(f"CRAM/BAM file not found: {cram_path}")
        run_picard_hs_metrics(
            cram_file=cram_path,
            output_file=hsm_file,
            bait_intervals=bait_intervals,
            target_intervals=target_intervals,
            reference=picard_ref,
            picard_xmx=xmx,
            log_file=out_dir / f"{sample_id}.picard.log",
            dry_run=dry_run,
        )

    if dry_run:
        logger.info("[DRY RUN] Skipping metrics parsing for %s", sample_id)
        return qcsum_file, None

    summary = evaluate_qcsum(sample_id, panel_name, parse_hs_metrics(hsm_file), panel_cfg, sample_type=sample_type)
    write_text_atomic(qcsum_file, f"{QC_SUM_HEADER}\n{summary.csv_line}\n")
    logger.info(
        "QC summary written: %s (Align=%s, Cov=%s, Mean=%.2fx, ROI@%dx=%.2f%%)",
        qcsum_file, summary.alignment_qc, summary.coverage_qc,
        summary.mean_target_coverage, summary.covered_threshold, summary.pct_roi_at_threshold,
    )
    return qcsum_file, summary


def merge_qcsum_files(
    qcsum_files: Iterable[str | Path],
    output_dir: str | Path | None = None,
    output_name: str = "merged_qcsum.csv",
) -> Path:
    """Merge per-sample .qcsum.txt files into one CSV, validating matching headers."""
    file_list = [Path(p) for p in qcsum_files]
    if not file_list:
        raise ValueError("No qcsum files provided to merge")

    output_file = (Path.cwd() if output_dir is None else Path(output_dir).resolve()) / output_name
    header: str | None = None
    data_lines: list[str] = []
    files_merged = 0

    for path in file_list:
        if not path.exists():
            logger.warning("Skipping missing file: %s", path)
            continue
        lines = path.read_text().splitlines()
        if not lines:
            logger.warning("Skipping empty file: %s", path)
            continue
        file_header = lines[0].strip()
        if header is None:
            header = file_header
        elif file_header != header:
            raise ValueError(f"Header mismatch in {path}: expected '{header}', got '{file_header}'")
        files_merged += 1
        data_lines.extend(line.strip() for line in lines[1:] if line.strip())

    if header is None:
        raise ValueError("No valid qcsum files found to merge")

    write_text_atomic(output_file, header + "\n" + "".join(f"{line}\n" for line in data_lines))
    logger.info("Merged %d samples from %d files into %s", len(data_lines), files_merged, output_file)
    return output_file


def format_summary_table(summaries: list[QCSumSummary]) -> str:
    """Render a compact human-readable PASS/FAIL QC summary table."""
    if not summaries:
        return "No samples summarized."
    headers = ("Sample", "Panel", "Type", "Align", "Cov", "Reads", "%Aln", "MeanCov", "%ROI")
    rows = [
        (
            s.sample_id, s.panel, s.sample_type or "-", s.alignment_qc, s.coverage_qc,
            f"{s.total_reads:,}", f"{s.pct_aligned:.2f}%",
            f"{s.mean_target_coverage:.1f}x", f"{s.pct_roi_at_threshold:.2f}%@{s.covered_threshold}x",
        )
        for s in summaries
    ]
    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    return "\n".join([fmt.format(*headers), "  ".join("-" * w for w in widths), *(fmt.format(*r) for r in rows)])


# --- Pipeline steps ---


def derive_run_name(run_dir: str | Path) -> str:
    """Extract the run directory basename cleanly regardless of path depth or trailing slashes."""
    name = os.path.basename(os.path.normpath(str(run_dir)))
    if not name or name in (".", "..", os.sep):
        raise ValueError(f"Could not derive run name from --run-dir: {run_dir!r}")
    return name


def get_samplesheets(run_dir: str | Path) -> list[Path]:
    """Return sorted SampleSheet_*.csv paths from a sequencer run directory."""
    run_path = Path(run_dir)
    if not run_path.is_dir():
        raise FileNotFoundError(f"Run directory not found: {run_path}")
    sheets = sorted(
        p for p in run_path.iterdir()
        if p.is_file() and p.name.startswith("SampleSheet_") and p.suffix == ".csv"
    )
    if not sheets:
        raise FileNotFoundError(f"No SampleSheet_*.csv files found in {run_path}")
    return sheets


def get_index(samplesheet: str | Path) -> str:
    """Extract the index group label from SampleSheet_<index>.csv."""
    return Path(samplesheet).stem.removeprefix("SampleSheet_")


def get_sample_ids(fastq_list: str | Path) -> list[str]:
    """Return sorted unique RGSM sample IDs from a DRAGEN fastq_list.csv."""
    with open(fastq_list, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "RGSM" not in reader.fieldnames:
            raise ValueError(f"Missing required 'RGSM' column in {fastq_list}")
        return sorted({row["RGSM"].strip() for row in reader if row.get("RGSM") and row["RGSM"].strip()})


def read_sampleinfo(sampleinfo: str | Path) -> list[dict[str, str]]:
    """Read non-blank sampleinfo TSV rows using stdlib csv (handling UTF-8 BOM)."""
    with open(sampleinfo, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "Samples" not in reader.fieldnames:
            raise ValueError(f"Missing required 'Samples' column in {sampleinfo}")
        return [row for row in reader if any(v and v.strip() for v in row.values())]


def demux_samples(run_dir: str | Path, fastqs_dir: str | Path, dry_run: bool = False) -> None:
    """Run DRAGEN BCLConvert for each SampleSheet_*.csv in run_dir."""
    fastqs_path = Path(fastqs_dir)
    if fastqs_path.exists():
        raise FileExistsError(f"Demux output directory already exists: {fastqs_path}")
    # Validate samplesheets before creating fastqs_dir so typos never poison retries
    sheets = get_samplesheets(run_dir)
    if not dry_run:
        fastqs_path.mkdir(parents=True, exist_ok=True)

    for sheet in sheets:
        run_cmd(
            [
                "dragen",
                "--bcl-conversion-only", "true",
                "--bcl-use-hw", "false",
                "--bcl-only-matched-reads", "true",
                "--bcl-input-directory", str(run_dir),
                "--sample-sheet", str(sheet),
                "--output-directory", str(fastqs_path / get_index(sheet)),
            ],
            dry_run=dry_run,
        )


def build_dragen_align_cmd(
    fastq_list: str | Path,
    bam_output: str | Path,
    sample_id: str,
    group_cfg: dict,
    dragen_tmp_dir: str,
) -> list[str]:
    """Build the DRAGEN alignment + variant calling command from an align_groups entry."""
    bed_file = str(group_cfg["bed"])
    cmd = [
        "dragen",
        "--intermediate-results-dir", dragen_tmp_dir,
        "--enable-map-align", "true",
        "--enable-map-align-output", "true",
        "--output-format", "CRAM",
        "--generate-sa-tags", "true",
        "--enable-sort", "true",
        "--soft-read-trimmers", "polyg,quality",
        "--trim-min-quality", "2",
        "--ref-dir", str(group_cfg["human_ref"]),
        "--qc-coverage-tag-1", "target_bed",
        "--qc-coverage-region-1", bed_file,
        "--qc-coverage-reports-1", "cov_report",
        "--qc-coverage-ignore-overlaps", "true",
        "--enable-variant-caller", "true",
        "--vc-combine-phased-variants-distance", "6",
        "--vc-emit-ref-confidence", "GVCF",
        "--enable-hla", "true",
        "--fastq-list", str(fastq_list),
        "--fastq-list-sample-id", sample_id,
        "--output-directory", str(bam_output),
        "--output-file-prefix", sample_id,
    ]
    cmd.extend(str(token).format(bed=bed_file) for token in (group_cfg.get("extra_dragen_args") or ()))
    return cmd


def align_samples(
    fastqs_dir: str | Path,
    bams_dir: str | Path,
    config: PanelsConfig,
    dry_run: bool = False,
    on_sample_aligned: Callable[[str], None] | None = None,
) -> None:
    """Align all non-skipped samples discovered under fastqs_dir."""
    fastqs_path = Path(fastqs_dir)
    if not dry_run and not fastqs_path.is_dir():
        raise FileNotFoundError(f"FASTQ directory not found: {fastqs_path}")

    fastq_lists = sorted(fastqs_path.rglob("fastq_list.csv")) if fastqs_path.is_dir() else []
    if not fastq_lists:
        logger.warning("No fastq_list.csv files found under %s", fastqs_path)
    dragen_tmp_dir = config.defaults.get("dragen_tmp_dir", "/staging/tmp")

    for fastq_list in fastq_lists:
        group_name, group_cfg = config.get_align_group(fastq_list, fastqs_dir=fastqs_path)
        if group_cfg.get("skip"):
            logger.info("Skipping alignment for group '%s': %s", group_name, fastq_list)
            continue

        for sample_id in get_sample_ids(fastq_list):
            bam_output = Path(bams_dir) / sample_id
            if bam_output.exists():
                logger.info("Skipping alignment for %s (output dir exists): %s", sample_id, bam_output)
            else:
                if not dry_run:
                    bam_output.mkdir(parents=True, exist_ok=True)
                cmd = build_dragen_align_cmd(fastq_list, bam_output, sample_id, group_cfg, dragen_tmp_dir)
                run_cmd(cmd, dry_run=dry_run, log_file=bam_output / f"{sample_id}.dragen.log")
            if on_sample_aligned is not None:
                on_sample_aligned(sample_id)


def build_occ_pf_mqc_payload(headers: list[str], rows: Iterable[Iterable[float]]) -> dict | None:
    """Build a MultiQC v1.23 custom-content scatter JSON payload for % Occupied vs % Pass Filter."""
    required = ("% Occupied", "% Pass Filter", "Lane")
    if any(col not in headers for col in required):
        return None

    x_idx = headers.index("% Occupied")
    y_idx = headers.index("% Pass Filter")
    lane_idx = headers.index("Lane")
    tile_idx = headers.index("Tile") if "Tile" in headers else None

    seen_tiles: set[tuple[int, int]] = set()
    lane_colors_used: dict[int, str] = {}
    points: dict[str, dict[str, object]] = {}
    for row_num, row in enumerate(rows):
        vals = list(row)
        occ, pf, lane_raw = float(vals[x_idx]), float(vals[y_idx]), float(vals[lane_idx])
        tile_raw = float(vals[tile_idx]) if tile_idx is not None else float(row_num)
        if not (math.isfinite(occ) and math.isfinite(pf) and math.isfinite(lane_raw) and math.isfinite(tile_raw)):
            continue
        lane, tile = int(lane_raw), int(tile_raw)
        if (lane, tile) in seen_tiles:
            continue
        seen_tiles.add((lane, tile))
        color = LANE_COLORS[(lane - 1) % len(LANE_COLORS)]
        lane_colors_used[lane] = color
        points[f"Lane {lane} Tile {tile}"] = {
            "x": round(occ, 3),
            "y": round(pf, 3),
            "color": color,
            "opacity": 0.5,
            "marker_size": 5,
        }

    if not points:
        return None

    legend = ", ".join(f"Lane {lane} ({lane_colors_used[lane]})" for lane in sorted(lane_colors_used))
    return {
        "id": "occ_pf_lane",
        "section_name": "Occupied vs Pass Filter",
        "description": f"Illumina InterOp % Occupied vs % Pass Filter per tile. Lane colors: {legend}.",
        "plot_type": "scatter",
        "pconfig": {
            "id": "occ_pf_lane_plot",
            "title": "% Occupied vs % Pass Filter",
            "xlab": "% Occupied",
            "ylab": "% Pass Filter",
            "xmin": 0,
            "xmax": 100,
            "ymin": 50,
            "ymax": 100,
        },
        "data": points,
    }


def save_occ_pf_plot(run_dir: str | Path, output_dir: str | Path) -> Path | None:
    """Emit native MultiQC scatter JSON (occ_pf_lane_mqc.json) from Illumina InterOp metrics."""
    from interop import py_interop_run, py_interop_run_metrics, py_interop_table
    from numpy import float32, zeros

    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    # Only Tile (% Pass Filter) and ExtendedTile (% Occupied) are needed; skip per-cycle Extraction
    valid_to_load[py_interop_run.ExtendedTile] = 1
    valid_to_load[py_interop_run.Tile] = 1

    try:
        run_metrics.read(str(run_dir), valid_to_load)
    except Exception as e:
        logger.warning("Unable to read InterOp files in %s (%s); skipping plot", run_dir, e)
        return None

    columns = py_interop_table.imaging_column_vector()
    py_interop_table.create_imaging_table_columns(run_metrics, columns)
    if columns.size() == 0:
        logger.warning("No InterOp imaging data found in %s; skipping plot", run_dir)
        return None

    headers: list[str] = []
    for i in range(columns.size()):
        col = columns[i]
        headers.extend(f"{col.name()} ({sub})" for sub in col.subcolumns()) if col.has_children() else headers.append(col.name())

    row_offsets = py_interop_table.map_id_offset()
    py_interop_table.count_table_rows(run_metrics, row_offsets)
    data = zeros((row_offsets.size(), py_interop_table.count_table_columns(columns)), dtype=float32)
    py_interop_table.populate_imaging_table_data(run_metrics, columns, row_offsets, data.ravel())

    payload = build_occ_pf_mqc_payload(headers, data)
    if payload is None:
        logger.warning("No valid InterOp Occupied/PF data found in %s; skipping plot", run_dir)
        return None

    json_path = Path(output_dir) / "occ_pf_lane_mqc.json"
    write_text_atomic(json_path, json.dumps(payload, allow_nan=False) + "\n")
    logger.info("Saved native MultiQC InterOp plot: %s (%d tiles)", json_path, len(payload["data"]))
    return json_path


def multiqc(fastqs_dir: str | Path, bams_dir: str | Path, dry_run: bool = False) -> None:
    """Remove irrelevant WGS DRAGEN coverage CSVs and run MultiQC."""
    for wgs_csv in sorted(Path(bams_dir).glob("*/*.wgs_*.csv")):
        if dry_run:
            logger.info("[DRY RUN] Would remove %s", wgs_csv)
            continue
        try:
            wgs_csv.unlink()
            logger.info("Removed irrelevant WGS coverage report: %s", wgs_csv)
        except FileNotFoundError:
            pass

    run_cmd(
        [
            "multiqc",
            "--force",
            "--config", str(MULTIQC_CONFIG_YAML),
            "--outdir", str(bams_dir),
            str(bams_dir),
            str(fastqs_dir),
        ],
        dry_run=dry_run,
    )


def _sampleinfo_qc_tasks(
    bams_dir: Path, sampleinfo: str | Path, config: PanelsConfig
) -> list[tuple[str, str, str, Path]]:
    """Parse and pre-validate all sampleinfo rows before starting alignment or QC workers."""
    tasks: list[tuple[str, str, str, Path]] = []
    seen_ids: set[str] = set()
    for row_num, row in enumerate(read_sampleinfo(sampleinfo), start=2):
        sample_id = (row.get("Samples") or "").strip()
        if not sample_id:
            raise ValueError(f"Missing sample ID on row {row_num} of {sampleinfo}")
        if sample_id in seen_ids:
            raise ValueError(f"Duplicate sample ID '{sample_id}' in {sampleinfo}")
        seen_ids.add(sample_id)

        panel = (row.get("Panel") or row.get("Tumor") or "").strip()
        panel_name, panel_cfg = config.get_qc_panel(panel)
        if not panel_cfg.get("bait_intervals") or not panel_cfg.get("target_intervals"):
            raise ValueError(f"Panel '{panel_name}' missing bait_intervals or target_intervals")

        bam_path = (row.get("BAM Path") or "").strip()
        if not bam_path:
            raise ValueError(f"Missing 'BAM Path' for sample '{sample_id}' in {sampleinfo}")
        tasks.append((bam_path, panel_name, sample_id, bams_dir / sample_id))
    return tasks


def run_qc_stage(
    run_dir: str | Path,
    fastqs_dir: str | Path,
    bams_dir: str | Path,
    sampleinfo: str | Path,
    config: PanelsConfig,
    align_first: bool = False,
    workers: int | None = None,
    reference: str | None = None,
    picard_xmx: str | None = None,
    dry_run: bool = False,
) -> list[QCSumSummary]:
    """Run QC (optionally overlapped with DRAGEN alignment), write qcsum_mqc.csv, and run MultiQC."""
    bams_path = Path(bams_dir)
    tasks = _sampleinfo_qc_tasks(bams_path, sampleinfo, config)
    if not tasks:
        if align_first:
            align_samples(fastqs_dir, bams_path, config, dry_run=dry_run)
        if not dry_run:
            save_occ_pf_plot(run_dir, bams_path)
            write_text_atomic(bams_path / "qcsum_mqc.csv", QC_SUM_HEADER + "\n")
        multiqc(fastqs_dir, bams_path, dry_run=dry_run)
        return []

    sample_to_idx = {sample_id: idx for idx, (_, _, sample_id, _) in enumerate(tasks)}
    max_workers = resolve_worker_count(len(tasks), workers, config)
    submitted: dict[int, Future[tuple[Path, QCSumSummary | None]]] = {}
    results: dict[int, tuple[Path, QCSumSummary | None]] = {}

    pool = ThreadPoolExecutor(max_workers=max_workers)
    try:
        def submit_idx(idx: int) -> None:
            if idx in submitted:
                return
            bam_path, panel, sample_id, sample_dir = tasks[idx]
            submitted[idx] = pool.submit(
                run_qcsum,
                cram_file=bam_path,
                panel=panel,
                config=config,
                sample_id=sample_id,
                output_dir=sample_dir,
                reference=reference,
                picard_xmx=picard_xmx,
                dry_run=dry_run,
            )

        if align_first:
            align_samples(
                fastqs_dir,
                bams_path,
                config,
                dry_run=dry_run,
                on_sample_aligned=lambda sid: submit_idx(sample_to_idx[sid]) if sid in sample_to_idx else None,
            )

        for idx in range(len(tasks)):
            submit_idx(idx)

        if not dry_run:
            save_occ_pf_plot(run_dir, bams_path)

        future_to_idx = {fut: idx for idx, fut in submitted.items()}
        for fut in as_completed(future_to_idx):
            results[future_to_idx[fut]] = fut.result()
    except BaseException:
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    else:
        pool.shutdown(wait=True)

    if not dry_run:
        merge_qcsum_files([results[i][0] for i in range(len(tasks))], output_dir=bams_path, output_name="qcsum_mqc.csv")

    multiqc(fastqs_dir, bams_path, dry_run=dry_run)
    if dry_run:
        return []
    return [results[i][1] for i in range(len(tasks)) if results[i][1] is not None]


# --- CLI ---


def _check_batch_collisions(crams: list[str], output_dir: str | None) -> None:
    seen: dict[Path, str] = {}
    for cram in crams:
        cram_path = Path(cram).resolve()
        out_dir = Path(output_dir).resolve() if output_dir else cram_path.parent
        target = out_dir / f"{cram_path.stem}.qcsum.txt"
        if target in seen:
            raise ValueError(f"Output collision for '{cram_path.stem}' ({target}): '{seen[target]}' and '{cram}'")
        seen[target] = cram


def build_parser(config: PanelsConfig) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="BCL-QC: demux, align, Picard QCSum, and MultiQC CLI.")
    parser.add_argument("--config", default=str(PANELS_YAML), help="Path to panels.yaml")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging")

    sub = parser.add_subparsers(dest="command", required=True)
    default_fastqs = config.defaults.get("fastqs_dir", "/staging/hot/reads")
    default_bams = config.defaults.get("bams_dir", "/mnt/pns/bams")
    default_ref = config.defaults.get("picard_ref")
    default_xmx = config.defaults.get("picard_xmx", "4g")
    default_workers = config.defaults.get("qc_workers", "4")

    def add_run_io(p: argparse.ArgumentParser, include_bams: bool = True) -> None:
        p.add_argument("--run-dir", required=True, help="Sequencer run directory")
        p.add_argument("--fastqs-dir", default=default_fastqs, help="Parent FASTQ output directory")
        if include_bams:
            p.add_argument("--bams-dir", default=default_bams, help="Parent BAM/CRAM output directory")
        p.add_argument("-n", "--dry-run", action="store_true", help="Log commands without executing")

    def add_picard_opts(p: argparse.ArgumentParser, short_ref: bool = False) -> None:
        p.add_argument("-t", "--workers", type=int, default=None, help=f"Max Picard workers (default: min(tasks, {default_workers}))")
        p.add_argument(*("-r", "--picard-ref") if short_ref else ("--picard-ref",), default=default_ref, help="Reference FASTA for Picard")
        p.add_argument("--picard-xmx", default=default_xmx, help=f"JVM heap size for Picard (default: {default_xmx})")

    p_run = sub.add_parser("run", help="Run selected pipeline steps (demux, align, qc)")
    add_run_io(p_run)
    p_run.add_argument("--sampleinfo", help="Path to sampleinfo TSV (required when 'qc' is in --steps)")
    p_run.add_argument("--steps", nargs="+", choices=DEF_STEPS, default=list(DEF_STEPS), help="Steps to execute")
    add_picard_opts(p_run)

    add_run_io(sub.add_parser("demux", help="Run DRAGEN BCLConvert demultiplexing"), include_bams=False)
    add_run_io(sub.add_parser("align", help="Run DRAGEN alignment and variant calling"))

    p_qc = sub.add_parser("qc", help="Run InterOp plot, Picard QCSum, and MultiQC")
    add_run_io(p_qc)
    p_qc.add_argument("--sampleinfo", required=True, help="Path to sampleinfo TSV")
    add_picard_opts(p_qc)

    p_qcsum = sub.add_parser("qcsum", help="Run standalone Picard CollectHsMetrics + QCSum on 1+ CRAM/BAMs")
    p_qcsum.add_argument("crams", nargs="+", metavar="CRAM", help="Input CRAM/BAM path(s)")
    p_qcsum.add_argument("-p", "--panel", required=True, choices=list(config.qc_panels), help="QC panel")
    p_qcsum.add_argument("-s", "--sample-id", help="Sample ID (only valid for a single CRAM)")
    p_qcsum.add_argument("-o", "--output-dir", help="Output directory (defaults to each CRAM's parent)")
    add_picard_opts(p_qcsum, short_ref=True)
    p_qcsum.add_argument("--merge-output", default="merged_qcsum.csv", help="Merged CSV filename for multi-CRAM runs")
    p_qcsum.add_argument("-n", "--dry-run", action="store_true", help="Log commands without executing")

    p_merge = sub.add_parser("merge", help="Merge multiple .qcsum.txt files into one CSV")
    p_merge.add_argument("files", nargs="+", metavar="FILE", help="Input .qcsum.txt files")
    p_merge.add_argument("-o", "--output-dir", help="Output directory (default: cwd)")
    p_merge.add_argument("--output", default="merged_qcsum.csv", help="Merged output filename")

    sub.add_parser("panels", help="List configured alignment groups and QC panels")
    return parser


def main(argv: list[str] | None = None) -> int:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=str(PANELS_YAML))
    pre.add_argument("-v", "--verbose", action="store_true")
    pre_args, _ = pre.parse_known_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if pre_args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    try:
        config = PanelsConfig(pre_args.config)
    except (FileNotFoundError, ValueError, yaml.YAMLError) as e:
        logger.error("%s", e)
        return 1

    parser = build_parser(config)
    args = parser.parse_args(argv)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if getattr(args, "workers", None) is not None and args.workers < 1:
        parser.error("--workers must be >= 1")
    if getattr(args, "picard_xmx", None) is not None:
        try:
            normalize_picard_xmx(args.picard_xmx)
        except ValueError as e:
            parser.error(str(e))

    try:
        if args.command == "panels":
            print("Alignment groups:")
            for name, grp in config.align_groups.items():
                detail = "skipped" if grp.get("skip") else f"{grp['panel']} ({grp['human_ref']})"
                print(f"  {name:12s} - {detail}")
            print("\nQC panels:")
            for name, pcfg in config.qc_panels.items():
                print(f"  {name:12s} - {pcfg.get('pipeline_version', '')} (cov>={pcfg.get('covered')}x, avg>={pcfg.get('fail_min_avgcov')}x)")
            return 0

        if args.command == "merge":
            print(f"Merged qcsum written to: {merge_qcsum_files(args.files, output_dir=args.output_dir, output_name=args.output)}")
            return 0

        if args.command == "qcsum":
            if args.sample_id and len(args.crams) > 1:
                parser.error("--sample-id can only be used with a single CRAM")
            if len(args.crams) > 1:
                _check_batch_collisions(args.crams, args.output_dir)

            max_workers = resolve_worker_count(len(args.crams), args.workers, config)
            results: dict[int, tuple[Path, QCSumSummary | None]] = {}
            failed: list[str] = []

            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futs = {
                    pool.submit(
                        run_qcsum,
                        cram_file=cram,
                        panel=args.panel,
                        config=config,
                        sample_id=args.sample_id,
                        output_dir=args.output_dir,
                        reference=args.picard_ref,
                        picard_xmx=args.picard_xmx,
                        dry_run=args.dry_run,
                    ): (idx, cram)
                    for idx, cram in enumerate(args.crams)
                }
                for fut in as_completed(futs):
                    idx, cram = futs[fut]
                    try:
                        results[idx] = fut.result()
                    except Exception as e:
                        logger.error("Failed to process %s: %s", cram, e)
                        failed.append(cram)

            if args.dry_run:
                return 1 if failed else 0

            ordered = [i for i in range(len(args.crams)) if i in results]
            summaries = [results[i][1] for i in ordered if results[i][1] is not None]
            if summaries:
                print("\n" + format_summary_table(summaries))
            if len(args.crams) > 1 and ordered:
                merged = merge_qcsum_files([results[i][0] for i in ordered], output_dir=args.output_dir, output_name=args.merge_output)
                print(f"\nMerged qcsum: {merged}")
            if failed:
                print(f"\nFailed ({len(failed)}): {', '.join(failed)}", file=sys.stderr)
            return 1 if failed else 0

        run_name = derive_run_name(args.run_dir)
        fastqs_dir = Path(args.fastqs_dir) / run_name
        steps = [args.command] if args.command in DEF_STEPS else list(args.steps)
        if "qc" in steps and not args.sampleinfo:
            parser.error("--sampleinfo is required when 'qc' is in --steps")

        if "demux" in steps:
            demux_samples(args.run_dir, fastqs_dir, dry_run=args.dry_run)
        if "align" in steps or "qc" in steps:
            bams_dir = Path(args.bams_dir) / run_name
            if "qc" in steps:
                summaries = run_qc_stage(
                    args.run_dir, fastqs_dir, bams_dir, args.sampleinfo, config,
                    align_first="align" in steps,
                    workers=args.workers, reference=args.picard_ref,
                    picard_xmx=args.picard_xmx, dry_run=args.dry_run,
                )
                if summaries:
                    print("\n" + format_summary_table(summaries))
            else:
                align_samples(fastqs_dir, bams_dir, config, dry_run=args.dry_run)
        return 0
    except (FileNotFoundError, FileExistsError, ValueError, yaml.YAMLError) as e:
        logger.error("%s", e)
        return 1
    except subprocess.CalledProcessError as e:
        return e.returncode


if __name__ == "__main__":
    sys.exit(main())
