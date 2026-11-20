"""Unit and integration tests for unified bclqc.py and config/panels.yaml."""

import contextlib
import csv
import io
import json
import subprocess
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock

import yaml

import bclqc

PICARD_HSM_HEADERS = [
    "BAIT_SET", "BAIT_TERRITORY", "BAIT_DESIGN_EFFICIENCY", "ON_BAIT_BASES",
    "NEAR_BAIT_BASES", "OFF_BAIT_BASES", "PCT_SELECTED_BASES", "PCT_OFF_BAIT",
    "ON_BAIT_VS_SELECTED", "MEAN_BAIT_COVERAGE", "PCT_USABLE_BASES_ON_BAIT",
    "PCT_USABLE_BASES_ON_TARGET", "FOLD_ENRICHMENT", "HS_LIBRARY_SIZE",
    "HS_PENALTY_10X", "HS_PENALTY_20X", "HS_PENALTY_30X", "HS_PENALTY_40X",
    "HS_PENALTY_50X", "HS_PENALTY_100X", "TARGET_TERRITORY", "GENOME_SIZE",
    "TOTAL_READS", "PF_READS", "PF_BASES", "PF_UNIQUE_READS",
    "PF_UQ_READS_ALIGNED", "PF_BASES_ALIGNED", "PF_UQ_BASES_ALIGNED",
    "ON_TARGET_BASES", "PCT_PF_READS", "PCT_PF_UQ_READS",
    "PCT_PF_UQ_READS_ALIGNED", "MEAN_TARGET_COVERAGE", "MEDIAN_TARGET_COVERAGE",
    "MAX_TARGET_COVERAGE", "MIN_TARGET_COVERAGE", "ZERO_CVG_TARGETS_PCT",
    "PCT_EXC_DUPE", "PCT_EXC_ADAPTER", "PCT_EXC_MAPQ", "PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP", "PCT_EXC_OFF_TARGET", "FOLD_80_BASE_PENALTY",
    "PCT_TARGET_BASES_1X", "PCT_TARGET_BASES_2X", "PCT_TARGET_BASES_10X",
    "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_40X",
    "PCT_TARGET_BASES_50X", "PCT_TARGET_BASES_100X", "PCT_TARGET_BASES_250X",
    "PCT_TARGET_BASES_500X", "PCT_TARGET_BASES_1000X", "AT_DROPOUT",
    "GC_DROPOUT", "HET_SNP_SENSITIVITY", "HET_SNP_Q", "SAMPLE", "LIBRARY",
    "READ_GROUP",
]


def make_hsm_text(**overrides: object) -> str:
    defaults: dict[str, object] = {
        "BAIT_SET": "test_baits",
        "BAIT_TERRITORY": "1000000",
        "BAIT_DESIGN_EFFICIENCY": "0.85",
        "ON_BAIT_BASES": "7200000000",       # 72.00% of PF_BASES_ALIGNED
        "NEAR_BAIT_BASES": "1600000000",
        "OFF_BAIT_BASES": "1200000000",
        "PCT_SELECTED_BASES": "0.880000",    # 88.00%
        "PCT_OFF_BAIT": "0.120000",
        "ON_BAIT_VS_SELECTED": "0.818182",
        "MEAN_BAIT_COVERAGE": "1523.456789",
        # Deliberately distinct from ON_BAIT/ON_TARGET / PF_BASES_ALIGNED
        "PCT_USABLE_BASES_ON_BAIT": "0.310000",
        "PCT_USABLE_BASES_ON_TARGET": "0.240000",
        "FOLD_ENRICHMENT": "450.12",
        "HS_LIBRARY_SIZE": "12345678",
        "HS_PENALTY_10X": "1.1", "HS_PENALTY_20X": "1.2", "HS_PENALTY_30X": "1.3",
        "HS_PENALTY_40X": "1.4", "HS_PENALTY_50X": "1.5", "HS_PENALTY_100X": "1.6",
        "TARGET_TERRITORY": "850000",
        "GENOME_SIZE": "3088269832",
        "TOTAL_READS": "80000000",
        "PF_READS": "80000000",
        "PF_BASES": "12000000000",
        "PF_UNIQUE_READS": "50000000",
        "PF_UQ_READS_ALIGNED": "49750000",
        "PF_BASES_ALIGNED": "10000000000",
        "PF_UQ_BASES_ALIGNED": "6000000000",
        "ON_TARGET_BASES": "6150000000",     # 61.50% of PF_BASES_ALIGNED
        "PCT_PF_READS": "1.000000",
        "PCT_PF_UQ_READS": "0.625000",
        "PCT_PF_UQ_READS_ALIGNED": "0.995000",
        "MEAN_TARGET_COVERAGE": "1345.678901",
        "MEDIAN_TARGET_COVERAGE": "1310",
        "MAX_TARGET_COVERAGE": "4200",
        "MIN_TARGET_COVERAGE": "12",
        "ZERO_CVG_TARGETS_PCT": "0.0001",
        "PCT_EXC_DUPE": "0.35", "PCT_EXC_ADAPTER": "0.001", "PCT_EXC_MAPQ": "0.01",
        "PCT_EXC_BASEQ": "0.01", "PCT_EXC_OVERLAP": "0.05", "PCT_EXC_OFF_TARGET": "0.10",
        "FOLD_80_BASE_PENALTY": "1.345678",
        "PCT_TARGET_BASES_1X": "0.999900",
        "PCT_TARGET_BASES_2X": "0.999800",
        "PCT_TARGET_BASES_10X": "0.999500",
        "PCT_TARGET_BASES_20X": "0.998000",
        "PCT_TARGET_BASES_30X": "0.997000",
        "PCT_TARGET_BASES_40X": "0.996000",
        "PCT_TARGET_BASES_50X": "0.995000",
        "PCT_TARGET_BASES_100X": "0.991234",
        "PCT_TARGET_BASES_250X": "0.965432",
        "PCT_TARGET_BASES_500X": "0.912345",
        "PCT_TARGET_BASES_1000X": "0.750000",
        "AT_DROPOUT": "1.234567",
        "GC_DROPOUT": "2.345678",
        "HET_SNP_SENSITIVITY": "0.999",
        "HET_SNP_Q": "30",
        "SAMPLE": "", "LIBRARY": "", "READ_GROUP": "",
    }
    defaults.update(overrides)
    header_line = "\t".join(PICARD_HSM_HEADERS)
    value_line = "\t".join(str(defaults[col]) for col in PICARD_HSM_HEADERS)
    return "\n".join([
        "## htsjdk.samtools.metrics.StringHeader",
        "# CollectHsMetrics INPUT=sample.cram OUTPUT=sample.hsm.txt",
        "## htsjdk.samtools.metrics.StringHeader",
        "# Started on: Fri Nov 20 11:00:00 PST 2026",
        "",
        "## METRICS CLASS\tpicard.analysis.directed.HsMetrics",
        header_line,
        value_line,
        "",
    ])


def parse_qcsum_row(csv_line: str) -> dict[str, str]:
    return next(csv.DictReader(io.StringIO(f"{bclqc.QC_SUM_HEADER}\n{csv_line}\n")))


def run_cli(argv: list[str]) -> tuple[int, str, str]:
    out, err = io.StringIO(), io.StringIO()
    with (
        mock.patch("logging.basicConfig"),
        contextlib.redirect_stdout(out),
        contextlib.redirect_stderr(err),
    ):
        rc = bclqc.main(argv)
    return rc, out.getvalue(), err.getvalue()


class TestPanelsConfig(unittest.TestCase):
    def setUp(self):
        self.cfg = bclqc.PanelsConfig()

    def test_panels_aliases_and_inheritance(self):
        self.assertEqual(
            set(self.cfg.qc_panels),
            {"pcp_tumor", "pcp_normal", "heme_comp", "heme_mpn", "heme_jak2", "heme_pblp"},
        )
        self.assertEqual(self.cfg.resolve_panel_name("Tumor"), "pcp_tumor")
        self.assertEqual(self.cfg.resolve_panel_name("Normal"), "pcp_normal")
        self.assertEqual(self.cfg.resolve_panel_name("Comprehensive Heme Panel"), "heme_comp")
        self.assertEqual(self.cfg.resolve_panel_name("MPN Screen Panel (JAK2, CALR, MPL)"), "heme_mpn")

        _, comp = self.cfg.get_qc_panel("heme_comp")
        _, mpn = self.cfg.get_qc_panel("heme_mpn")
        self.assertEqual(mpn["bait_intervals"], comp["bait_intervals"])
        self.assertIn("MPN.cds10bp.interval_list", mpn["target_intervals"])

        _, tumor = self.cfg.get_qc_panel("pcp_tumor")
        _, normal = self.cfg.get_qc_panel("pcp_normal")
        self.assertEqual(normal["bait_intervals"], tumor["bait_intervals"])
        self.assertEqual(normal["sample_type"], "Normal")
        self.assertEqual(normal["covered"], "20")

    def test_config_validation_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            # Circular inheritance
            (p / "cycle.yaml").write_text("qc_panels:\n  a:\n    inherits: b\n  b:\n    inherits: a\n")
            with self.assertRaises(ValueError):
                bclqc.PanelsConfig(p / "cycle.yaml")

            # Scalar alias rejected
            (p / "scalar.yaml").write_text("qc_panels:\n  a:\n    aliases: Tumor\n")
            with self.assertRaises(ValueError):
                bclqc.PanelsConfig(p / "scalar.yaml")

            # Duplicate aliases rejected
            (p / "dup.yaml").write_text("qc_panels:\n  a:\n    aliases: [X]\n  b:\n    aliases: [X]\n")
            with self.assertRaises(ValueError):
                bclqc.PanelsConfig(p / "dup.yaml")

            # Null panel entry handled cleanly
            (p / "null_entry.yaml").write_text("qc_panels:\n  empty_panel:\n")
            cfg = bclqc.PanelsConfig(p / "null_entry.yaml")
            self.assertEqual(cfg.qc_panels["empty_panel"], {})

    def test_align_group_innermost_match_ignores_ancestor_collision(self):
        # Ancestor directory named 'I10' must NOT hijack inner 'U5N2_I10'
        name, grp = self.cfg.get_align_group(
            "/mnt/I10/reads/RUN/U5N2_I10/Reports/fastq_list.csv",
            fastqs_dir="/mnt/I10/reads/RUN",
        )
        self.assertEqual(name, "U5N2_I10")
        self.assertEqual(grp["panel"], "HEME")

        name2, grp2 = self.cfg.get_align_group("/mnt/I10/reads/RUN/U5N2_I10/Reports/fastq_list.csv")
        self.assertEqual(name2, "U5N2_I10")
        self.assertEqual(grp2["panel"], "HEME")

        name_exome, grp_exome = self.cfg.get_align_group("/reads/RUN/I8N2_N10/Reports/fastq_list.csv")
        self.assertEqual(name_exome, "I8N2_N10")
        self.assertTrue(grp_exome["skip"])

        with self.assertRaises(ValueError):
            self.cfg.get_align_group("/reads/RUN/UNKNOWN/fastq_list.csv")

    def test_multiqc_and_conda_env_configs(self):
        with bclqc.MULTIQC_CONFIG_YAML.open() as f:
            mc = yaml.safe_load(f)
        self.assertNotIn("general_stats", mc["module_order"])
        self.assertEqual(mc["exclude_modules"], ["snippy"])
        self.assertEqual(mc["custom_content"]["order"], ["qcsum-module", "occ_pf_lane"])

        env_text = (bclqc.SCRIPT_DIR / "config" / "conda_env.yaml").read_text()
        self.assertIn("numpy", env_text)
        for dropped in ("r-base", "pandas", "seaborn", "matplotlib", "kaleido"):
            self.assertNotIn(dropped, env_text)

    def test_build_occ_pf_mqc_payload_and_edge_cases(self):
        headers = ["Lane", "Tile", "% Occupied", "% Pass Filter"]
        rows = [
            [1, 1101, 85.1234, 92.5678],
            [1, 1101, 85.1234, 92.5678],   # duplicate tile -> deduplicated
            [1, 1102, float("nan"), 90.0],  # NaN occ -> skipped
            [float("nan"), 1103, 80.0, 90.0],  # NaN lane -> skipped
            [2, float("inf"), 80.0, 90.0],  # Inf tile -> skipped
            [2, 2101, 80.5, 88.25],
        ]
        payload = bclqc.build_occ_pf_mqc_payload(headers, rows)
        assert payload is not None
        self.assertEqual(payload["id"], "occ_pf_lane")
        self.assertEqual(payload["plot_type"], "scatter")
        self.assertEqual(payload["pconfig"]["ylab"], "% Pass Filter")
        self.assertIn("Lane 1 (#1f77b4)", payload["description"])
        self.assertIn("Lane 2 (#ff7f0e)", payload["description"])
        self.assertEqual(set(payload["data"].keys()), {"Lane 1 Tile 1101", "Lane 2 Tile 2101"})
        for pt in payload["data"].values():
            self.assertIsInstance(pt, dict)
        json.dumps(payload, allow_nan=False)

        # Missing % Occupied header -> None (no crash)
        self.assertIsNone(bclqc.build_occ_pf_mqc_payload(["Lane", "% Pass Filter"], [[1, 90.0]]))
        # All-NaN rows -> None
        self.assertIsNone(bclqc.build_occ_pf_mqc_payload(headers, [[1, 1101, float("nan"), 90.0]]))
        # Missing Tile column falls back cleanly to row index
        no_tile = bclqc.build_occ_pf_mqc_payload(["Lane", "% Occupied", "% Pass Filter"], [[1, 85.0, 91.0]])
        assert no_tile is not None
        self.assertIn("Lane 1 Tile 0", no_tile["data"])


class TestHsMetricsAndQCSum(unittest.TestCase):
    def setUp(self):
        self.cfg = bclqc.PanelsConfig()

    def _metrics(self, **overrides: object) -> dict[str, str]:
        with tempfile.TemporaryDirectory() as tmp:
            hsm = Path(tmp) / "S1.hsm.txt"
            hsm.write_text(make_hsm_text(**overrides))
            return bclqc.parse_hs_metrics(hsm)

    def test_parse_question_marks_and_errors(self):
        metrics = self._metrics(
            TOTAL_READS="?", FOLD_80_BASE_PENALTY="?",
            MEDIAN_TARGET_COVERAGE="?", MAX_TARGET_COVERAGE="?",
        )
        _, pcfg = self.cfg.get_qc_panel("pcp_tumor")
        row = parse_qcsum_row(bclqc.evaluate_qcsum("S1", "pcp_tumor", metrics, pcfg).csv_line)
        self.assertEqual(row["Total_Reads"], "0")
        self.assertEqual(row["FOLD_80_BASE_PENALTY"], "0.00")
        self.assertEqual(row["MEDIAN_ROI_COVERAGE"], "0")

        with tempfile.TemporaryDirectory() as tmp:
            trunc = Path(tmp) / "trunc.hsm.txt"
            trunc.write_text("\t".join(PICARD_HSM_HEADERS) + "\n")
            with self.assertRaises(ValueError):
                bclqc.parse_hs_metrics(trunc)

    def test_formulas_formatting_thresholds_and_missing_roi_col(self):
        metrics = self._metrics()
        _, pcfg = self.cfg.get_qc_panel("pcp_tumor")
        summary = bclqc.evaluate_qcsum("SAMPLE_T", "pcp_tumor", metrics, pcfg)
        row = parse_qcsum_row(summary.csv_line)

        self.assertEqual(row["Tumor_Normal"], "Tumor")
        self.assertEqual(row["Capture"], "twist")
        self.assertEqual(row["%On/Near_Bait_Bases"], "88.00")
        self.assertEqual(row["%On_Bait_Bases"], "72.00")
        self.assertEqual(row["%On_Target_Bases"], "61.50")
        self.assertEqual(row["AT_DROPOUT"], "1.23")
        self.assertEqual(row["GC_DROPOUT"], "2.35")
        self.assertEqual((row["Alignment_QC"], row["Coverage_QC"]), ("PASS", "PASS"))

        zero_row = parse_qcsum_row(
            bclqc.evaluate_qcsum(
                "S0", "pcp_tumor",
                self._metrics(PF_BASES_ALIGNED="0", ON_BAIT_BASES="0", ON_TARGET_BASES="0"),
                pcfg,
            ).csv_line
        )
        self.assertEqual((zero_row["%On_Bait_Bases"], zero_row["%On_Target_Bases"]), ("0.00", "0.00"))

        at_bound = self._metrics(
            PCT_PF_UQ_READS_ALIGNED="0.970000", TOTAL_READS="66000000",
            PCT_TARGET_BASES_250X="0.900000", MEAN_TARGET_COVERAGE="400.0",
        )
        res = bclqc.evaluate_qcsum("S1", "pcp_tumor", at_bound, pcfg)
        self.assertEqual((res.alignment_qc, res.coverage_qc), ("PASS", "PASS"))

        low_cov = bclqc.evaluate_qcsum("S1", "pcp_tumor", self._metrics(MEAN_TARGET_COVERAGE="399.9"), pcfg)
        self.assertEqual(low_cov.coverage_qc, "FAIL")

        missing_col = dict(metrics)
        del missing_col["PCT_TARGET_BASES_250X"]
        with self.assertRaises(ValueError):
            bclqc.evaluate_qcsum("S1", "pcp_tumor", missing_col, pcfg)

    def test_run_qcsum_edge_cases(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            # Empty string CRAM rejected
            with self.assertRaises(ValueError):
                bclqc.run_qcsum("", "pcp_tumor", self.cfg, output_dir=p)
            # Directory as CRAM rejected
            with self.assertRaises(FileNotFoundError):
                bclqc.run_qcsum(p, "pcp_tumor", self.cfg, output_dir=p)
            # Missing CRAM file rejected
            with self.assertRaises(FileNotFoundError):
                bclqc.run_qcsum(p / "missing.cram", "pcp_tumor", self.cfg, output_dir=p)
            # Invalid panel does not create output_dir
            bad_out = p / "should_not_exist"
            with self.assertRaises(ValueError):
                bclqc.run_qcsum(p / "S1.cram", "no_such_panel", self.cfg, output_dir=bad_out)
            self.assertFalse(bad_out.exists())

            # Dry-run does not create output_dir
            dry_out = p / "dry_out"
            bclqc.run_qcsum(p / "S1.cram", "pcp_tumor", self.cfg, output_dir=dry_out, dry_run=True)
            self.assertFalse(dry_out.exists())

            # 0-byte .hsm.txt triggers Picard rerun
            cram = p / "S1.cram"
            cram.write_bytes(b"cram")
            (p / "S1.hsm.txt").write_text("")

            def fake_picard(cram_file, output_file, **kwargs):
                Path(output_file).write_text(make_hsm_text())

            with mock.patch.object(bclqc, "run_picard_hs_metrics", side_effect=fake_picard) as m_pic:
                _, summary = bclqc.run_qcsum(cram, "heme_comp", self.cfg)
                m_pic.assert_called_once()
                assert summary is not None
                self.assertEqual(summary.alignment_qc, "PASS")

    def test_merge_and_header_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            f1, f2, bad = p / "S1.qcsum.txt", p / "S2.qcsum.txt", p / "BAD.qcsum.txt"
            f1.write_text(f"{bclqc.QC_SUM_HEADER}\nS1,Tumor,rest\n")
            f2.write_text(f"{bclqc.QC_SUM_HEADER}\nS2,Normal,rest\n")
            bad.write_text("Sample,Wrong\nS3,x\n")

            merged = bclqc.merge_qcsum_files([f1, f2], output_dir=p)
            self.assertEqual(merged.read_text().splitlines(), [bclqc.QC_SUM_HEADER, "S1,Tumor,rest", "S2,Normal,rest"])
            with self.assertRaises(ValueError):
                bclqc.merge_qcsum_files([f1, bad], output_dir=p)


class TestPipelineAndCLI(unittest.TestCase):
    def setUp(self):
        self.cfg = bclqc.PanelsConfig()

    def test_derive_run_name_and_tsv_bom_handling(self):
        self.assertEqual(bclqc.derive_run_name("/mnt/pns/runs/RUN_01/"), "RUN_01")
        self.assertEqual(bclqc.derive_run_name("RUN_02"), "RUN_02")
        for bad in ("/", ".", "..", ""):
            with self.assertRaises(ValueError):
                bclqc.derive_run_name(bad)

        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            # UTF-8 BOM + trailing blank line
            si = p / "bom.tsv"
            si.write_bytes(b"\xef\xbb\xbfSamples\tTumor\tBAM Path\nS1\tTumor\t/a/S1.cram\n\t\t\n")
            rows = bclqc.read_sampleinfo(si)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["Samples"], "S1")

            # Missing Samples header
            bad_si = p / "bad.tsv"
            bad_si.write_text("ID\tTumor\nS1\tTumor\n")
            with self.assertRaises(ValueError):
                bclqc.read_sampleinfo(bad_si)

            # Missing RGSM column in fastq_list.csv
            bad_fq = p / "fq.csv"
            bad_fq.write_text("Sample,Read1File\nS1,r1.fq.gz\n")
            with self.assertRaises(ValueError):
                bclqc.get_sample_ids(bad_fq)

    def test_demux_samples_validation_before_mkdir(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            empty_run = p / "empty_run"
            empty_run.mkdir()
            fastqs = p / "reads" / "empty_run"

            # Missing samplesheets must NOT create fastqs_dir
            with self.assertRaises(FileNotFoundError):
                bclqc.demux_samples(empty_run, fastqs)
            self.assertFalse(fastqs.exists())

            # Strict SampleSheet_*.csv matching & sorting
            (empty_run / "SampleSheet_U5N2_I10.csv").write_text("x")
            (empty_run / "SampleSheet_I10.csv").write_text("x")
            (empty_run / "SampleSheet_I10.csv.bak").write_text("ignore")
            (empty_run / "other_SampleSheet_X.csv").write_text("ignore")

            cmds = []
            with mock.patch.object(bclqc, "run_cmd", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                bclqc.demux_samples(empty_run, fastqs)
            self.assertEqual(len(cmds), 2)
            self.assertTrue(cmds[0][cmds[0].index("--sample-sheet") + 1].endswith("SampleSheet_I10.csv"))
            self.assertTrue(cmds[1][cmds[1].index("--sample-sheet") + 1].endswith("SampleSheet_U5N2_I10.csv"))

            # Existing fastqs_dir raises FileExistsError
            with self.assertRaises(FileExistsError):
                bclqc.demux_samples(empty_run, fastqs)

    def test_align_samples_cmds_and_missing_fastqs_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                bclqc.align_samples(Path(tmp) / "missing_fastqs", Path(tmp) / "bams", self.cfg)

            fastqs = Path(tmp) / "fastqs"
            bams = Path(tmp) / "bams"
            for grp, samples in (("I10", ["S_B", "S_A"]), ("U5N2_I10", ["H1"]), ("I8N2_N10", ["EX1"])):
                rep = fastqs / grp / "Reports"
                rep.mkdir(parents=True)
                (rep / "fastq_list.csv").write_text("RGSM,Read1File\n" + "".join(f"{s},r1.fq.gz\n" for s in samples))
            # Pre-existing dir is skipped (non-resume-safe semantics preserved)
            (bams / "S_A").mkdir(parents=True)

            cmds = []
            with mock.patch.object(bclqc, "run_cmd", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                bclqc.align_samples(fastqs, bams, self.cfg)

            self.assertEqual(len(cmds), 2)
            pcp_cmd, heme_cmd = cmds
            self.assertIn("S_B", pcp_cmd)
            self.assertIn("--enable-duplicate-marking", pcp_cmd)
            self.assertIn("H1", heme_cmd)
            self.assertIn("--umi-enable", heme_cmd)
            self.assertIn(self.cfg.align_groups["U5N2_I10"]["bed"], heme_cmd)

    def test_multiqc_removes_wgs_csvs(self):
        with tempfile.TemporaryDirectory() as tmp:
            bams, fastqs = Path(tmp) / "bams", Path(tmp) / "fastqs"
            sdir = bams / "S1"
            sdir.mkdir(parents=True)
            fastqs.mkdir()
            wgs = sdir / "S1.wgs_coverage_metrics.csv"
            keep = sdir / "S1.target_bed_coverage_metrics.csv"
            wgs.write_text("wgs\n")
            keep.write_text("keep\n")

            with mock.patch.object(bclqc, "run_cmd") as mock_run:
                bclqc.multiqc(fastqs, bams)
                mock_run.assert_called_once()
                self.assertEqual(mock_run.call_args[0][0][0], "multiqc")

            self.assertFalse(wgs.exists())
            self.assertTrue(keep.exists())

    def test_picard_jvm_flags_bounded_workers_and_log_streaming(self):
        self.assertEqual(bclqc.normalize_picard_xmx("4g"), "-Xmx4g")
        self.assertEqual(bclqc.normalize_picard_xmx("-Xmx8G"), "-Xmx8G")
        for bad in ("0g", "bad", "-Xmx"):
            with self.assertRaises(ValueError):
                bclqc.normalize_picard_xmx(bad)

        with mock.patch("os.cpu_count", return_value=64):
            self.assertEqual(bclqc.resolve_worker_count(20, None, self.cfg), 4)
            self.assertEqual(bclqc.resolve_worker_count(2, None, self.cfg), 2)
            self.assertEqual(bclqc.resolve_worker_count(20, 8, self.cfg), 8)

        captured = {}
        with mock.patch.object(bclqc, "run_cmd", side_effect=lambda cmd, **kw: captured.update(cmd=cmd, **kw)):
            bclqc.run_picard_hs_metrics(
                "in.cram", "/tmp/out/S1.hsm.txt", "b.ilist", "t.ilist", "/nvme/hg38.fa",
                picard_xmx="8g", log_file="/tmp/out/S1.picard.log",
            )
        self.assertEqual(captured["cmd"][:4], ["picard", "-Xmx8g", "-XX:ActiveProcessorCount=2", "CollectHsMetrics"])
        self.assertIn("R=/nvme/hg38.fa", captured["cmd"])
        self.assertEqual(captured["log_file"], "/tmp/out/S1.picard.log")

        # Streaming success & failure paths
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "sub.log"
            bclqc.run_cmd(["python3", "-c", "print('hello'); import sys; print('err', file=sys.stderr)"], log_file=log)
            self.assertIn("hello", log.read_text())
            self.assertIn("err", log.read_text())

            fail_log = Path(tmp) / "fail.log"
            with self.assertRaises(subprocess.CalledProcessError):
                bclqc.run_cmd(["python3", "-c", "print('boom'); import sys; sys.exit(3)"], log_file=fail_log)
            self.assertIn("boom", fail_log.read_text())

    def test_overlapped_align_and_qc_concurrent_picard_and_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            run_dir = tmp_path / "RUN_1"
            run_dir.mkdir()
            fastqs = tmp_path / "fastqs" / "RUN_1"
            bams = tmp_path / "bams" / "RUN_1"
            rep = fastqs / "I10" / "Reports"
            rep.mkdir(parents=True)
            (rep / "fastq_list.csv").write_text("RGSM,Read1File\nS1,r1.fq.gz\nS2,r2.fq.gz\n")

            # Duplicate sample ID in sampleinfo must fail BEFORE alignment starts
            dup_si = tmp_path / "dup.tsv"
            dup_si.write_text(f"Samples\tTumor\tBAM Path\nS1\tTumor\t{bams}/S1/S1.cram\nS1\tTumor\t{bams}/S1/S1.cram\n")
            with mock.patch.object(bclqc, "align_samples") as m_align:
                with self.assertRaises(ValueError):
                    bclqc.run_qc_stage(run_dir, fastqs, bams, dup_si, self.cfg, align_first=True)
                m_align.assert_not_called()

            # Valid overlapped run with historical normal + real Picard invocation concurrently with DRAGEN
            old_cram = tmp_path / "old" / "N_HIST.cram"
            old_cram.parent.mkdir()
            old_cram.write_bytes(b"cram")
            si = tmp_path / "pcp.tsv"
            si.write_text(
                "Samples\tTumor\tBAM Path\n"
                f"S2\tTumor\t{bams}/S2/S2.cram\n"
                f"N_HIST\tNormal\t{old_cram}\n"
                f"S1\tTumor\t{bams}/S1/S1.cram\n"
            )

            s1_picard_started = threading.Event()
            multiqc_called = False

            def fake_run_cmd(cmd, dry_run=False, log_file=None):
                nonlocal multiqc_called
                if cmd[0] == "dragen":
                    sid = cmd[cmd.index("--fastq-list-sample-id") + 1]
                    if sid == "S2":
                        s1_picard_started.wait(timeout=2.0)
                    out_dir = Path(cmd[cmd.index("--output-directory") + 1])
                    (out_dir / f"{sid}.cram").write_bytes(b"cram")
                elif cmd[0] == "picard":
                    out_hsm = Path(next(a.split("=", 1)[1] for a in cmd if a.startswith("O=")))
                    if out_hsm.name == "S1.hsm.txt":
                        s1_picard_started.set()
                    out_hsm.write_text(make_hsm_text())
                elif cmd[0] == "multiqc":
                    multiqc_called = True

            with (
                mock.patch.object(bclqc, "run_cmd", side_effect=fake_run_cmd),
                mock.patch.object(bclqc, "save_occ_pf_plot", return_value=None),
            ):
                summaries = bclqc.run_qc_stage(run_dir, fastqs, bams, si, self.cfg, align_first=True, workers=2)

            self.assertTrue(s1_picard_started.is_set(), "S1 Picard should start before S2 DRAGEN finishes")
            self.assertTrue(multiqc_called, "multiqc should run after overlapped align+qc")
            self.assertEqual([s.sample_id for s in summaries], ["S2", "N_HIST", "S1"])
            with (bams / "qcsum_mqc.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([(r["Sample"], r["Tumor_Normal"]) for r in rows], [
                ("S2", "Tumor"), ("N_HIST", "Normal"), ("S1", "Tumor"),
            ])

            # Empty sampleinfo writes header-only qcsum_mqc.csv
            empty_si = tmp_path / "empty.tsv"
            empty_si.write_text("Samples\tTumor\tBAM Path\n")
            bams_empty = tmp_path / "bams_empty"
            with (
                mock.patch.object(bclqc, "run_cmd"),
                mock.patch.object(bclqc, "save_occ_pf_plot", return_value=None),
            ):
                self.assertEqual(bclqc.run_qc_stage(run_dir, fastqs, bams_empty, empty_si, self.cfg), [])
            self.assertEqual((bams_empty / "qcsum_mqc.csv").read_text(), bclqc.QC_SUM_HEADER + "\n")

            # Mid-run DRAGEN failure propagates cleanly without hanging
            bams_fail = tmp_path / "bams_fail"
            si_fail = tmp_path / "fail.tsv"
            si_fail.write_text(f"Samples\tTumor\tBAM Path\nS1\tTumor\t{bams_fail}/S1/S1.cram\nS2\tTumor\t{bams_fail}/S2/S2.cram\n")

            def fail_on_s2(cmd, dry_run=False, log_file=None):
                if cmd[0] == "dragen":
                    sid = cmd[cmd.index("--fastq-list-sample-id") + 1]
                    if sid == "S2":
                        raise subprocess.CalledProcessError(1, cmd)
                    Path(cmd[cmd.index("--output-directory") + 1], f"{sid}.cram").write_bytes(b"cram")
                elif cmd[0] == "picard":
                    Path(next(a.split("=", 1)[1] for a in cmd if a.startswith("O="))).write_text(make_hsm_text())

            with mock.patch.object(bclqc, "run_cmd", side_effect=fail_on_s2):
                with self.assertRaises(subprocess.CalledProcessError):
                    bclqc.run_qc_stage(run_dir, fastqs, bams_fail, si_fail, self.cfg, align_first=True, workers=1)

    def test_cli_subcommands_validation_and_partial_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            out = p / "out"
            out.mkdir()
            crams = []
            for name in ("S_C", "S_A", "S_B"):
                (out / f"{name}.hsm.txt").write_text(make_hsm_text())
                crams.append(str(p / f"{name}.cram"))

            rc, stdout, _ = run_cli(["qcsum", *crams, "--panel", "pcp_normal", "-o", str(out), "-t", "2"])
            self.assertEqual(rc, 0)
            self.assertIn("S_C", stdout)
            with (out / "merged_qcsum.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([r["Sample"] for r in rows], ["S_C", "S_A", "S_B"])

            # Partial failure prints Failed summary and exits 1
            rc_part, _, stderr_part = run_cli([
                "qcsum", crams[0], str(p / "missing.cram"), "--panel", "pcp_normal", "-o", str(out),
            ])
            self.assertEqual(rc_part, 1)
            self.assertIn("Failed (1):", stderr_part)

            # qcsum --dry-run must return 1 if a task raises
            bad_cfg = p / "no_ref.yaml"
            bad_cfg.write_text("defaults: {}\nqc_panels:\n  p:\n    bait_intervals: b\n    target_intervals: t\n")
            rc_dry, _, _ = run_cli(["--config", str(bad_cfg), "qcsum", str(p / "x.cram"), "-p", "p", "--dry-run"])
            self.assertEqual(rc_dry, 1)

            # Missing --config exits 1 cleanly (no raw traceback)
            rc_cfg, _, _ = run_cli(["--config", str(p / "nope.yaml"), "panels"])
            self.assertEqual(rc_cfg, 1)

            # panels & merge subcommands
            rc_p, out_p, _ = run_cli(["panels"])
            self.assertEqual(rc_p, 0)
            self.assertIn("pcp_tumor", out_p)

            rc_m, _, _ = run_cli(["merge", str(out / "S_A.qcsum.txt"), str(out / "S_B.qcsum.txt"), "-o", str(out)])
            self.assertEqual(rc_m, 0)

            # Argparse validation errors
            for bad_argv in (
                ["qcsum", *crams, "--panel", "pcp_normal", "-s", "X"],
                ["run", "--run-dir", "/runs/R1", "--steps", "qc"],
                ["qcsum", crams[0], "--panel", "pcp_normal", "--picard-xmx", "bad"],
                ["qcsum", crams[0], "--panel", "pcp_normal", "-t", "0"],
            ):
                with self.assertRaises(SystemExit) as ctx:
                    run_cli(bad_argv)
                self.assertEqual(ctx.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
