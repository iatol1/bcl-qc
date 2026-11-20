"""Unit and integration tests for unified bclqc.py and config/panels.yaml."""

import csv
import io
import tempfile
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


class TestPanelsConfig(unittest.TestCase):
    def setUp(self):
        self.cfg = bclqc.PanelsConfig()

    def test_panels_and_aliases(self):
        self.assertEqual(
            set(self.cfg.qc_panels),
            {"pcp_tumor", "pcp_normal", "heme_comp", "heme_mpn", "heme_jak2", "heme_pblp"},
        )
        self.assertEqual(self.cfg.resolve_panel_name("Tumor"), "pcp_tumor")
        self.assertEqual(self.cfg.resolve_panel_name("Normal"), "pcp_normal")
        self.assertEqual(self.cfg.resolve_panel_name("Comprehensive Heme Panel"), "heme_comp")
        self.assertEqual(self.cfg.resolve_panel_name("MPN Screen Panel (JAK2, CALR, MPL)"), "heme_mpn")

    def test_inheritance_and_overrides(self):
        _, comp = self.cfg.get_qc_panel("heme_comp")
        _, mpn = self.cfg.get_qc_panel("heme_mpn")
        self.assertEqual(mpn["bait_intervals"], comp["bait_intervals"])
        self.assertEqual(mpn["pipeline_version"], comp["pipeline_version"])
        self.assertIn("MPN.cds10bp.interval_list", mpn["target_intervals"])

        _, tumor = self.cfg.get_qc_panel("pcp_tumor")
        _, normal = self.cfg.get_qc_panel("pcp_normal")
        self.assertEqual(normal["bait_intervals"], tumor["bait_intervals"])
        self.assertEqual(normal["sample_type"], "Normal")
        self.assertEqual(normal["covered"], "20")

    def test_circular_inheritance_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.yaml"
            bad.write_text("qc_panels:\n  a:\n    inherits: b\n  b:\n    inherits: a\n")
            with self.assertRaises(ValueError):
                bclqc.PanelsConfig(bad)

    def test_align_group_matching(self):
        name, grp = self.cfg.get_align_group("/reads/RUN/U5N2_I10/Reports/fastq_list.csv")
        self.assertEqual(name, "U5N2_I10")
        self.assertEqual(grp["panel"], "HEME")
        name_exome, grp_exome = self.cfg.get_align_group("/reads/RUN/I8N2_N10/Reports/fastq_list.csv")
        self.assertEqual(name_exome, "I8N2_N10")
        self.assertTrue(grp_exome["skip"])
        with self.assertRaises(ValueError):
            self.cfg.get_align_group("/reads/RUN/UNKNOWN/fastq_list.csv")

    def test_multiqc_config_valid(self):
        with bclqc.MULTIQC_CONFIG_YAML.open() as f:
            mc = yaml.safe_load(f)
        self.assertNotIn("general_stats", mc["module_order"])
        self.assertEqual(mc["exclude_modules"], ["snippy"])


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

    def test_formulas_formatting_and_thresholds(self):
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
        self.assertEqual(row["Alignment_QC"], "PASS")
        self.assertEqual(row["Coverage_QC"], "PASS")

        # Zero PF_BASES_ALIGNED
        zero_row = parse_qcsum_row(
            bclqc.evaluate_qcsum(
                "S0", "pcp_tumor",
                self._metrics(PF_BASES_ALIGNED="0", ON_BAIT_BASES="0", ON_TARGET_BASES="0"),
                pcfg,
            ).csv_line
        )
        self.assertEqual(zero_row["%On_Bait_Bases"], "0.00")
        self.assertEqual(zero_row["%On_Target_Bases"], "0.00")

        # Exact boundary vs below boundary
        at_bound = self._metrics(
            PCT_PF_UQ_READS_ALIGNED="0.970000", TOTAL_READS="66000000",
            PCT_TARGET_BASES_250X="0.900000", MEAN_TARGET_COVERAGE="400.0",
        )
        res = bclqc.evaluate_qcsum("S1", "pcp_tumor", at_bound, pcfg)
        self.assertEqual((res.alignment_qc, res.coverage_qc), ("PASS", "PASS"))

        low_cov = bclqc.evaluate_qcsum("S1", "pcp_tumor", self._metrics(MEAN_TARGET_COVERAGE="399.9"), pcfg)
        self.assertEqual(low_cov.coverage_qc, "FAIL")

    def test_merge_and_header_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            f1 = p / "S1.qcsum.txt"
            f2 = p / "S2.qcsum.txt"
            bad = p / "BAD.qcsum.txt"
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

    def test_samplesheets_strict_and_sorted(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = Path(tmp)
            (r / "SampleSheet_U5N2_I10.csv").write_text("x")
            (r / "SampleSheet_I10.csv").write_text("x")
            (r / "SampleSheet_I10.csv.bak").write_text("ignore")
            (r / "other_SampleSheet_X.csv").write_text("ignore")
            sheets = bclqc.get_samplesheets(r)
            self.assertEqual([s.name for s in sheets], ["SampleSheet_I10.csv", "SampleSheet_U5N2_I10.csv"])
            self.assertEqual([bclqc.get_index(s) for s in sheets], ["I10", "U5N2_I10"])

    def test_align_samples_builds_panel_specific_dragen_cmds_and_skips_exome(self):
        with tempfile.TemporaryDirectory() as tmp:
            fastqs = Path(tmp) / "fastqs"
            bams = Path(tmp) / "bams"
            for grp, samples in (("I10", ["S_B", "S_A"]), ("U5N2_I10", ["H1"]), ("I8N2_N10", ["EX1"])):
                rep = fastqs / grp / "Reports"
                rep.mkdir(parents=True)
                lines = ["RGSM,Read1File"] + [f"{s},r1.fq.gz" for s in samples]
                (rep / "fastq_list.csv").write_text("\n".join(lines) + "\n")
            # Pre-existing dir should be skipped (preserving non-resume-safe semantics)
            (bams / "S_A").mkdir(parents=True)

            cmds = []
            with mock.patch.object(bclqc, "run_cmd", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                bclqc.align_samples(fastqs, bams, self.cfg)

            self.assertEqual(len(cmds), 2)
            # S_B (PCP) has duplicate marking, H1 (HEME) has UMI flags with resolved BED
            pcp_cmd, heme_cmd = cmds
            self.assertIn("S_B", pcp_cmd)
            self.assertIn("--enable-duplicate-marking", pcp_cmd)
            self.assertIn("H1", heme_cmd)
            self.assertIn("--umi-enable", heme_cmd)
            self.assertIn(self.cfg.align_groups["U5N2_I10"]["bed"], heme_cmd)

    def test_multiqc_removes_wgs_csvs(self):
        with tempfile.TemporaryDirectory() as tmp:
            bams = Path(tmp) / "bams"
            fastqs = Path(tmp) / "fastqs"
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

    def test_cli_qcsum_batch_order_collisions_and_qc_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            out = p / "out"
            out.mkdir()
            crams = []
            for name in ("S_C", "S_A", "S_B"):
                (out / f"{name}.hsm.txt").write_text(make_hsm_text())
                crams.append(str(p / f"{name}.cram"))

            rc = bclqc.main(["qcsum", *crams, "--panel", "pcp_normal", "-o", str(out), "-t", "2"])
            self.assertEqual(rc, 0)
            with (out / "merged_qcsum.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([r["Sample"] for r in rows], ["S_C", "S_A", "S_B"])
            self.assertEqual([r["Tumor_Normal"] for r in rows], ["Normal", "Normal", "Normal"])

            # Validation checks
            with self.assertRaises(SystemExit) as ctx:
                bclqc.main(["qcsum", *crams, "--panel", "pcp_normal", "-s", "X"])
            self.assertEqual(ctx.exception.code, 2)

            with self.assertRaises(SystemExit) as ctx:
                bclqc.main(["run", "--run-dir", "/runs/R1", "--steps", "qc"])
            self.assertEqual(ctx.exception.code, 2)

            # Collision check
            rc_dup = bclqc.main([
                "qcsum", str(p / "d1" / "DUP.cram"), str(p / "d2" / "DUP.cram"),
                "--panel", "pcp_tumor", "-o", str(out),
            ])
            self.assertEqual(rc_dup, 1)

            # panels subcommand
            self.assertEqual(bclqc.main(["panels"]), 0)

    def test_qc_sampleinfo_pcp_and_heme(self):
        with tempfile.TemporaryDirectory() as tmp:
            bams = Path(tmp) / "bams"
            for sample in ("T1", "N1", "H1"):
                sdir = bams / sample
                sdir.mkdir(parents=True)
                (sdir / f"{sample}.hsm.txt").write_text(make_hsm_text())

            pcp_si = Path(tmp) / "pcp.tsv"
            pcp_si.write_text(f"Samples\tTumor\tBAM Path\nT1\tTumor\t{bams}/T1/T1.cram\nN1\tNormal\t{bams}/N1/N1.cram\n")
            summaries = bclqc.run_qcsum_for_sampleinfo(bams, pcp_si, self.cfg)
            self.assertEqual(len(summaries), 2)
            self.assertIn("SAMPLE", bclqc.format_summary_table(summaries).upper())
            with (bams / "qcsum_mqc.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([(r["Sample"], r["Tumor_Normal"], r["Capture"]) for r in rows], [
                ("T1", "Tumor", "twist"), ("N1", "Normal", "twist"),
            ])

            heme_si = Path(tmp) / "heme.tsv"
            heme_si.write_text(f"Samples\tPanel\tBAM Path\nH1\tComprehensive Heme Panel\t{bams}/H1/H1.cram\n")
            bclqc.run_qcsum_for_sampleinfo(bams, heme_si, self.cfg)
            with (bams / "qcsum_mqc.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([(r["Sample"], r["Tumor_Normal"], r["Capture"]) for r in rows], [("H1", "", "goal")])


if __name__ == "__main__":
    unittest.main()
