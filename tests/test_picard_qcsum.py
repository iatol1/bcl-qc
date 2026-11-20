"""Unit and integration tests for picard_qcsum.py and its integration with bclqc.py."""

import csv
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import picard_qcsum


# Full Picard 3.3.0 CollectHsMetrics header columns (0-indexed to match Picard HsMetrics.java)
PICARD_HSM_HEADERS = [
    "BAIT_SET",                     # 0
    "BAIT_TERRITORY",               # 1
    "BAIT_DESIGN_EFFICIENCY",       # 2
    "ON_BAIT_BASES",                # 3
    "NEAR_BAIT_BASES",              # 4
    "OFF_BAIT_BASES",               # 5
    "PCT_SELECTED_BASES",           # 6
    "PCT_OFF_BAIT",                 # 7
    "ON_BAIT_VS_SELECTED",          # 8
    "MEAN_BAIT_COVERAGE",           # 9
    "PCT_USABLE_BASES_ON_BAIT",     # 10
    "PCT_USABLE_BASES_ON_TARGET",   # 11
    "FOLD_ENRICHMENT",              # 12
    "HS_LIBRARY_SIZE",              # 13
    "HS_PENALTY_10X",               # 14
    "HS_PENALTY_20X",               # 15
    "HS_PENALTY_30X",               # 16
    "HS_PENALTY_40X",               # 17
    "HS_PENALTY_50X",               # 18
    "HS_PENALTY_100X",              # 19
    "TARGET_TERRITORY",             # 20
    "GENOME_SIZE",                  # 21
    "TOTAL_READS",                  # 22
    "PF_READS",                     # 23
    "PF_BASES",                     # 24
    "PF_UNIQUE_READS",              # 25
    "PF_UQ_READS_ALIGNED",          # 26
    "PF_BASES_ALIGNED",             # 27
    "PF_UQ_BASES_ALIGNED",          # 28
    "ON_TARGET_BASES",              # 29
    "PCT_PF_READS",                 # 30
    "PCT_PF_UQ_READS",              # 31
    "PCT_PF_UQ_READS_ALIGNED",      # 32
    "MEAN_TARGET_COVERAGE",         # 33
    "MEDIAN_TARGET_COVERAGE",       # 34
    "MAX_TARGET_COVERAGE",          # 35
    "MIN_TARGET_COVERAGE",          # 36
    "ZERO_CVG_TARGETS_PCT",         # 37
    "PCT_EXC_DUPE",                 # 38
    "PCT_EXC_ADAPTER",              # 39
    "PCT_EXC_MAPQ",                 # 40
    "PCT_EXC_BASEQ",                # 41
    "PCT_EXC_OVERLAP",              # 42
    "PCT_EXC_OFF_TARGET",           # 43
    "FOLD_80_BASE_PENALTY",         # 44
    "PCT_TARGET_BASES_1X",          # 45
    "PCT_TARGET_BASES_2X",          # 46
    "PCT_TARGET_BASES_10X",         # 47
    "PCT_TARGET_BASES_20X",         # 48
    "PCT_TARGET_BASES_30X",         # 49
    "PCT_TARGET_BASES_40X",         # 50
    "PCT_TARGET_BASES_50X",         # 51
    "PCT_TARGET_BASES_100X",        # 52
    "PCT_TARGET_BASES_250X",        # 53
    "PCT_TARGET_BASES_500X",        # 54
    "PCT_TARGET_BASES_1000X",       # 55
    "AT_DROPOUT",                   # 56
    "GC_DROPOUT",                   # 57
    "HET_SNP_SENSITIVITY",          # 58
    "HET_SNP_Q",                    # 59
    "SAMPLE",                       # 60
    "LIBRARY",                      # 61
    "READ_GROUP",                   # 62
]


def make_hsm_text(**overrides: object) -> str:
    """Build a realistic Picard 3.3.0 CollectHsMetrics .hsm.txt file."""
    defaults: dict[str, object] = {
        "BAIT_SET": "test_baits",
        "BAIT_TERRITORY": "1000000",
        "BAIT_DESIGN_EFFICIENCY": "0.85",
        "ON_BAIT_BASES": "7200000000",       # 72.00% of PF_BASES_ALIGNED
        "NEAR_BAIT_BASES": "1600000000",     # 16.00% of PF_BASES_ALIGNED
        "OFF_BAIT_BASES": "1200000000",      # 12.00% of PF_BASES_ALIGNED
        "PCT_SELECTED_BASES": "0.880000",    # 88.00% on/near bait
        "PCT_OFF_BAIT": "0.120000",
        "ON_BAIT_VS_SELECTED": "0.818182",
        "MEAN_BAIT_COVERAGE": "1523.456789",
        # Deliberately distinct from ON_BAIT_BASES / PF_BASES_ALIGNED to catch wrong column use
        "PCT_USABLE_BASES_ON_BAIT": "0.310000",
        "PCT_USABLE_BASES_ON_TARGET": "0.240000",
        "FOLD_ENRICHMENT": "450.12",
        "HS_LIBRARY_SIZE": "12345678",
        "HS_PENALTY_10X": "1.1",
        "HS_PENALTY_20X": "1.2",
        "HS_PENALTY_30X": "1.3",
        "HS_PENALTY_40X": "1.4",
        "HS_PENALTY_50X": "1.5",
        "HS_PENALTY_100X": "1.6",
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
        "PCT_PF_UQ_READS_ALIGNED": "0.995000",  # 99.50%
        "MEAN_TARGET_COVERAGE": "1345.678901",
        "MEDIAN_TARGET_COVERAGE": "1310",
        "MAX_TARGET_COVERAGE": "4200",
        "MIN_TARGET_COVERAGE": "12",
        "ZERO_CVG_TARGETS_PCT": "0.0001",
        "PCT_EXC_DUPE": "0.35",
        "PCT_EXC_ADAPTER": "0.001",
        "PCT_EXC_MAPQ": "0.01",
        "PCT_EXC_BASEQ": "0.01",
        "PCT_EXC_OVERLAP": "0.05",
        "PCT_EXC_OFF_TARGET": "0.10",
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
        "SAMPLE": "",
        "LIBRARY": "",
        "READ_GROUP": "",
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
    reader = csv.DictReader(io.StringIO(f"{picard_qcsum.QC_SUM_HEADER}\n{csv_line}\n"))
    return next(reader)


class TestQCSumConfig(unittest.TestCase):
    def test_config_path_is_relative_to_repo(self):
        self.assertTrue(picard_qcsum.QCSUM_CONFIG_YAML.is_file())
        self.assertEqual(
            picard_qcsum.QCSUM_CONFIG_YAML,
            Path(picard_qcsum.__file__).parent.resolve() / "config" / "qcsum_config.yaml",
        )

    def test_available_panels(self):
        panels = picard_qcsum.QCSumConfig.available_panels()
        self.assertEqual(
            set(panels),
            {"pcp_tumor", "pcp_normal", "heme_comp", "heme_mpn", "heme_jak2", "heme_pblp"},
        )

    def test_heme_subpanel_inherits_without_mutating_parent(self):
        comp = picard_qcsum.QCSumConfig("heme_comp")
        mpn = picard_qcsum.QCSumConfig("heme_mpn")
        self.assertEqual(mpn.get("bait_intervals"), comp.get("bait_intervals"))
        self.assertEqual(mpn.get("pipeline_version"), comp.get("pipeline_version"))
        self.assertNotEqual(mpn.get("target_intervals"), comp.get("target_intervals"))
        self.assertIn("MPN.cds10bp.interval_list", mpn.get("target_intervals"))

    def test_unknown_panel_raises(self):
        with self.assertRaises(ValueError):
            picard_qcsum.QCSumConfig("nonexistent_panel")


class TestParseHsMetrics(unittest.TestCase):
    def test_parse_standard_picard_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            hsm = Path(tmp) / "S1.hsm.txt"
            hsm.write_text(make_hsm_text())
            metrics = picard_qcsum.parse_hs_metrics(hsm)
            self.assertEqual(metrics["TOTAL_READS"], "80000000")
            self.assertEqual(metrics["ON_BAIT_BASES"], "7200000000")
            self.assertEqual(metrics["PF_BASES_ALIGNED"], "10000000000")
            self.assertEqual(metrics["ON_TARGET_BASES"], "6150000000")
            self.assertEqual(metrics["SAMPLE"], "0")

    def test_parse_question_marks_as_strings(self):
        with tempfile.TemporaryDirectory() as tmp:
            hsm = Path(tmp) / "S1.hsm.txt"
            hsm.write_text(
                make_hsm_text(
                    TOTAL_READS="?",
                    FOLD_80_BASE_PENALTY="?",
                    MEDIAN_TARGET_COVERAGE="?",
                    MAX_TARGET_COVERAGE="?",
                )
            )
            metrics = picard_qcsum.parse_hs_metrics(hsm)
            self.assertEqual(metrics["TOTAL_READS"], "0")
            self.assertIsInstance(metrics["TOTAL_READS"], str)
            self.assertEqual(metrics["MEDIAN_TARGET_COVERAGE"], "0")
            # Must not raise TypeError in ",".join(fields)
            cfg = picard_qcsum.QCSumConfig("pcp_tumor")
            line = picard_qcsum.generate_qcsum_line("S1", metrics, cfg)
            row = parse_qcsum_row(line)
            self.assertEqual(row["Total_Reads"], "0")
            self.assertEqual(row["FOLD_80_BASE_PENALTY"], "0.00")
            self.assertEqual(row["MEDIAN_ROI_COVERAGE"], "0")
            self.assertEqual(row["MAX_ROI_COVERAGE"], "0")

    def test_truncated_hsm_raises_value_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            hsm = Path(tmp) / "truncated.hsm.txt"
            hsm.write_text("\t".join(PICARD_HSM_HEADERS) + "\n")
            with self.assertRaises(ValueError):
                picard_qcsum.parse_hs_metrics(hsm)

    def test_missing_header_raises_value_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            hsm = Path(tmp) / "invalid.hsm.txt"
            hsm.write_text("## empty metrics\n")
            with self.assertRaises(ValueError):
                picard_qcsum.parse_hs_metrics(hsm)


class TestEvaluateQCSum(unittest.TestCase):
    def _metrics(self, **overrides: object) -> dict[str, str]:
        with tempfile.TemporaryDirectory() as tmp:
            hsm = Path(tmp) / "S1.hsm.txt"
            hsm.write_text(make_hsm_text(**overrides))
            return picard_qcsum.parse_hs_metrics(hsm)

    def test_on_bait_and_on_target_formulas_and_formatting(self):
        metrics = self._metrics()
        cfg = picard_qcsum.QCSumConfig("pcp_tumor")
        summary = picard_qcsum.evaluate_qcsum("SAMPLE_T", metrics, cfg)
        row = parse_qcsum_row(summary.csv_line)

        self.assertEqual(row["Sample"], "SAMPLE_T")
        self.assertEqual(row["Tumor_Normal"], "Tumor")
        self.assertEqual(row["Sequencing_Platform"], "NovaSeq6000")
        self.assertEqual(row["Pipeline_version"], "PCPTNDx-v2.0")
        self.assertEqual(row["Alignment_QC"], "PASS")
        self.assertEqual(row["Coverage_QC"], "PASS")
        self.assertEqual(row["Total_Reads"], "80000000")
        self.assertEqual(row["%Reads_Aligned"], "99.50")
        # Matches main: capture only, not capture-version
        self.assertEqual(row["Capture"], "twist")
        self.assertEqual(row["Avg_Capture_Coverage"], "1523.46")
        self.assertEqual(row["%On/Near_Bait_Bases"], "88.00")
        # ON_BAIT_BASES / PF_BASES_ALIGNED * 100 = 72.00 (NOT PCT_USABLE_BASES_ON_BAIT = 31.00)
        self.assertEqual(row["%On_Bait_Bases"], "72.00")
        # ON_TARGET_BASES / PF_BASES_ALIGNED * 100 = 61.50 (NOT PCT_USABLE_BASES_ON_TARGET = 24.00)
        self.assertEqual(row["%On_Target_Bases"], "61.50")
        self.assertEqual(row["FOLD_80_BASE_PENALTY"], "1.35")
        self.assertEqual(row["AT_DROPOUT"], "1.23")
        self.assertEqual(row["GC_DROPOUT"], "2.35")
        self.assertEqual(row["Avg_ROI_Coverage"], "1345.68")
        self.assertEqual(row["MEDIAN_ROI_COVERAGE"], "1310")
        self.assertEqual(row["MAX_ROI_COVERAGE"], "4200")
        self.assertEqual(row["%ROI_1x"], "99.99")
        self.assertEqual(row["%ROI_20x"], "99.80")
        self.assertEqual(row["%ROI_100x"], "99.12")
        self.assertEqual(row["%ROI_250x"], "96.54")
        self.assertEqual(row["%ROI_500x"], "91.23")

    def test_sample_type_defaults_across_panels(self):
        metrics = self._metrics()
        expected = {
            "pcp_tumor": "Tumor",
            "pcp_normal": "Normal",
            "heme_comp": "",
            "heme_mpn": "",
            "heme_jak2": "",
            "heme_pblp": "",
        }
        for panel, exp_type in expected.items():
            cfg = picard_qcsum.QCSumConfig(panel)
            row = parse_qcsum_row(picard_qcsum.generate_qcsum_line("S1", metrics, cfg))
            self.assertEqual(row["Tumor_Normal"], exp_type, f"panel={panel}")

    def test_zero_pf_bases_aligned_safe(self):
        metrics = self._metrics(PF_BASES_ALIGNED="0", ON_BAIT_BASES="0", ON_TARGET_BASES="0")
        cfg = picard_qcsum.QCSumConfig("pcp_tumor")
        row = parse_qcsum_row(picard_qcsum.generate_qcsum_line("S1", metrics, cfg))
        self.assertEqual(row["%On_Bait_Bases"], "0.00")
        self.assertEqual(row["%On_Target_Bases"], "0.00")

    def test_alignment_and_coverage_threshold_boundaries(self):
        # pcp_tumor: fail_min_align_pct=97.0, fail_min_reads=66000000, covered=250, fail_min_roi_pct=90, fail_min_avgcov=400
        cfg = picard_qcsum.QCSumConfig("pcp_tumor")

        # Exactly at boundary -> PASS (strict <)
        at_bound = self._metrics(
            PCT_PF_UQ_READS_ALIGNED="0.970000",
            TOTAL_READS="66000000",
            PCT_TARGET_BASES_250X="0.900000",
            MEAN_TARGET_COVERAGE="400.0",
        )
        res = picard_qcsum.evaluate_qcsum("S1", at_bound, cfg)
        self.assertEqual(res.alignment_qc, "PASS")
        self.assertEqual(res.coverage_qc, "PASS")

        # Below align pct -> Alignment FAIL
        low_align = self._metrics(PCT_PF_UQ_READS_ALIGNED="0.969900")
        self.assertEqual(picard_qcsum.evaluate_qcsum("S1", low_align, cfg).alignment_qc, "FAIL")

        # Below total reads -> Alignment FAIL
        low_reads = self._metrics(TOTAL_READS="65999999")
        self.assertEqual(picard_qcsum.evaluate_qcsum("S1", low_reads, cfg).alignment_qc, "FAIL")

        # Below ROI pct -> Coverage FAIL
        low_roi = self._metrics(PCT_TARGET_BASES_250X="0.899900")
        self.assertEqual(picard_qcsum.evaluate_qcsum("S1", low_roi, cfg).coverage_qc, "FAIL")

        # Below avg cov -> Coverage FAIL
        low_cov = self._metrics(MEAN_TARGET_COVERAGE="399.99")
        self.assertEqual(picard_qcsum.evaluate_qcsum("S1", low_cov, cfg).coverage_qc, "FAIL")

    def test_missing_configured_coverage_column_raises(self):
        metrics = self._metrics()
        del metrics["PCT_TARGET_BASES_250X"]
        cfg = picard_qcsum.QCSumConfig("pcp_tumor")
        with self.assertRaises(ValueError):
            picard_qcsum.evaluate_qcsum("S1", metrics, cfg)


class TestRunAndMergeAndCLI(unittest.TestCase):
    def test_run_qcsum_reuses_existing_hsm_without_cram(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "out"
            out_dir.mkdir()
            hsm = out_dir / "S1.hsm.txt"
            hsm.write_text(make_hsm_text())

            with mock.patch.object(picard_qcsum, "run_picard_hs_metrics") as mock_picard:
                result = picard_qcsum.run_qcsum(
                    cram_file=str(Path(tmp) / "nonexistent" / "S1.cram"),
                    panel="pcp_tumor",
                    output_dir=str(out_dir),
                )
                mock_picard.assert_not_called()

            self.assertEqual(result["sample_id"], "S1")
            self.assertEqual(result["sample_type"], "Tumor")
            self.assertEqual(result["alignment_qc"], "PASS")
            self.assertEqual(result["coverage_qc"], "PASS")
            qcsum_lines = Path(result["qcsum_file"]).read_text().splitlines()
            self.assertEqual(len(qcsum_lines), 2)
            self.assertEqual(qcsum_lines[0], picard_qcsum.QC_SUM_HEADER)

    def test_run_qcsum_reruns_picard_when_hsm_is_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            cram = tmp_path / "S1.cram"
            cram.write_bytes(b"fake")
            hsm = tmp_path / "S1.hsm.txt"
            hsm.write_text("")  # 0-byte leftover

            def fake_picard(cram_file, output_file, **kwargs):
                Path(output_file).write_text(make_hsm_text())

            with mock.patch.object(picard_qcsum, "run_picard_hs_metrics", side_effect=fake_picard) as mock_picard:
                result = picard_qcsum.run_qcsum(cram_file=str(cram), panel="heme_comp")
                mock_picard.assert_called_once()
            self.assertEqual(result["alignment_qc"], "PASS")
            self.assertEqual(result["coverage_qc"], "PASS")

    def test_merge_qcsum_files_and_header_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            f1 = tmp_path / "S1.qcsum.txt"
            f2 = tmp_path / "S2.qcsum.txt"
            bad = tmp_path / "BAD.qcsum.txt"
            f1.write_text(f"{picard_qcsum.QC_SUM_HEADER}\nS1,Tumor,rest\n")
            f2.write_text(f"{picard_qcsum.QC_SUM_HEADER}\nS2,Normal,rest\n")
            bad.write_text("Sample,DifferentHeader\nS3,val\n")

            merged = Path(picard_qcsum.merge_qcsum_files([f1, f2], output_dir=tmp_path))
            self.assertEqual(
                merged.read_text().splitlines(),
                [picard_qcsum.QC_SUM_HEADER, "S1,Tumor,rest", "S2,Normal,rest"],
            )

            with self.assertRaises(ValueError):
                picard_qcsum.merge_qcsum_files([f1, bad], output_dir=tmp_path)

    def test_cli_batch_preserves_order_and_validates_args(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            out_dir = tmp_path / "out"
            out_dir.mkdir()
            crams = []
            for name in ("S_C", "S_A", "S_B"):
                (out_dir / f"{name}.hsm.txt").write_text(make_hsm_text())
                crams.append(str(tmp_path / f"{name}.cram"))

            rc = picard_qcsum.main([
                "--batch", *crams,
                "--panel", "pcp_normal",
                "--output-dir", str(out_dir),
                "--threads", "2",
            ])
            self.assertEqual(rc, 0)
            merged = out_dir / "merged_qcsum.csv"
            with merged.open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([r["Sample"] for r in rows], ["S_C", "S_A", "S_B"])
            self.assertEqual([r["Tumor_Normal"] for r in rows], ["Normal", "Normal", "Normal"])

            # Invalid --threads
            self.assertEqual(picard_qcsum.main(["--batch", *crams, "--panel", "pcp_normal", "--threads", "0"]), 1)
            # --sample-id with --batch rejected
            self.assertEqual(
                picard_qcsum.main(["--batch", *crams, "--panel", "pcp_normal", "--sample-id", "X"]),
                1,
            )
            # Duplicate basename collision into same --output-dir rejected
            dup_crams = [str(tmp_path / "d1" / "DUP.cram"), str(tmp_path / "d2" / "DUP.cram")]
            self.assertEqual(
                picard_qcsum.main(["--batch", *dup_crams, "--panel", "pcp_normal", "-o", str(out_dir)]),
                1,
            )

    def test_cli_list_panels_and_dry_run(self):
        self.assertEqual(picard_qcsum.main(["--list-panels"]), 0)
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "dry_out"
            rc = picard_qcsum.main([
                "--cram", str(Path(tmp) / "S1.cram"),
                "--panel", "pcp_tumor",
                "--output-dir", str(out_dir),
                "--dry-run",
            ])
            self.assertEqual(rc, 0)
            self.assertFalse(out_dir.exists())


class TestConfigFiles(unittest.TestCase):
    def test_multiqc_config_valid(self):
        import yaml

        cfg_path = picard_qcsum.SCRIPT_DIR / "config" / "multiqc_config.yaml"
        with cfg_path.open() as f:
            cfg = yaml.safe_load(f)
        self.assertIn("module_order", cfg)
        self.assertNotIn("general_stats", cfg["module_order"])
        self.assertEqual(cfg["exclude_modules"], ["snippy"])
        self.assertEqual(
            cfg["custom_content"]["order"],
            ["qcsum-module", "occ_pf_lane-module"],
        )


class TestBclqcIntegration(unittest.TestCase):
    @staticmethod
    def _import_bclqc():
        import importlib
        import sys
        import types

        stubs = {}
        for mod_name in ("matplotlib", "matplotlib.pyplot", "seaborn", "interop", "numpy"):
            if mod_name not in sys.modules:
                try:
                    importlib.import_module(mod_name)
                except ImportError:
                    mod = types.ModuleType(mod_name)
                    if mod_name == "interop":
                        mod.py_interop_run_metrics = None
                        mod.py_interop_run = None
                        mod.py_interop_table = None
                    elif mod_name == "numpy":
                        mod.zeros = None
                        mod.float32 = None
                    stubs[mod_name] = mod

        try:
            importlib.import_module("pandas")
        except ImportError:
            pd_mod = types.ModuleType("pandas")

            class _SimpleDF:
                def __init__(self, rows):
                    self._rows = rows

                def iterrows(self):
                    for idx, row in enumerate(self._rows):
                        yield idx, row

            def _read_csv(path, sep=","):
                with open(path) as f:
                    return _SimpleDF(list(csv.DictReader(f, delimiter=sep)))

            pd_mod.read_csv = _read_csv
            stubs["pandas"] = pd_mod

        sys.modules.update(stubs)
        return importlib.import_module("bclqc")

    def test_bclqc_qcsum_pcp_and_heme(self):
        bclqc = self._import_bclqc()
        with tempfile.TemporaryDirectory() as tmp:
            bams = Path(tmp) / "bams"
            bams.mkdir()
            for sample in ("T1", "N1", "H1"):
                sdir = bams / sample
                sdir.mkdir()
                (sdir / f"{sample}.hsm.txt").write_text(make_hsm_text())

            pcp_si = Path(tmp) / "pcp_sampleinfo.tsv"
            pcp_si.write_text(
                "Samples\tTumor\tBAM Path\n"
                f"T1\tTumor\t{bams}/T1/T1.cram\n"
                f"N1\tNormal\t{bams}/N1/N1.cram\n"
            )
            bclqc.qcsum(str(bams), str(pcp_si))
            with (bams / "qcsum_mqc.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(
                [(r["Sample"], r["Tumor_Normal"], r["Capture"], r["%On_Bait_Bases"], r["%On_Target_Bases"]) for r in rows],
                [
                    ("T1", "Tumor", "twist", "72.00", "61.50"),
                    ("N1", "Normal", "twist", "72.00", "61.50"),
                ],
            )

            heme_si = Path(tmp) / "heme_sampleinfo.tsv"
            heme_si.write_text(
                "Samples\tPanel\tBAM Path\n"
                f"H1\tComprehensive Heme Panel\t{bams}/H1/H1.cram\n"
            )
            bclqc.qcsum(str(bams), str(heme_si))
            with (bams / "qcsum_mqc.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(
                [(r["Sample"], r["Tumor_Normal"], r["Capture"], r["%On_Bait_Bases"], r["%On_Target_Bases"]) for r in rows],
                [("H1", "", "goal", "72.00", "61.50")],
            )


if __name__ == "__main__":
    unittest.main()
