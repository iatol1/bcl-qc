"""Unit tests for Illumina SampleSheet filtering in samplesheet.py."""

import os
import tempfile
import unittest

from samplesheet import filter_samplesheet_lines, write_filtered_samplesheet


SAMPLE_SHEET_V2 = """\
[Header]
FileFormatVersion,2

[Reads]
Read1Cycles,151
Read2Cycles,151

[BCLConvert_Settings]
BarcodeMismatchesIndex1,1

[BCLConvert_Data]
Sample_ID,index,index2
25RR-001,CAGAGAGG,TATCCTCT
25RR-002,GCTACGCT,GTAAGGAG
25RR-Q01,AAGAGGCA,CTCTCTAT
"""

SAMPLE_SHEET_V1 = """\
[Header],,,
IEMFileVersion,5,,

[Data],,,
Lane,Sample_ID,index,index2
1,25RR-001,CAGAGAGG,TATCCTCT
1,25RR-002,GCTACGCT,GTAAGGAG
2,25RR-001,CAGAGAGG,TATCCTCT
"""


class TestSampleSheetFilter(unittest.TestCase):
    def test_v2_single_sample_preserves_headers_and_filters_rows(self):
        lines = SAMPLE_SHEET_V2.splitlines(keepends=True)
        filtered, matched = filter_samplesheet_lines(lines, ["25RR-002"])
        self.assertEqual(matched, {"25RR-002"})
        text = "".join(filtered)
        self.assertIn("[Header]", text)
        self.assertIn("[BCLConvert_Settings]", text)
        self.assertIn("Sample_ID,index,index2", text)
        self.assertIn("25RR-002,GCTACGCT,GTAAGGAG", text)
        self.assertNotIn("25RR-001", text)
        self.assertNotIn("25RR-Q01", text)

    def test_v1_non_first_column_sample_id_and_multi_lane(self):
        lines = SAMPLE_SHEET_V1.splitlines(keepends=True)
        filtered, matched = filter_samplesheet_lines(lines, ["25RR-001"])
        self.assertEqual(matched, {"25RR-001"})
        text = "".join(filtered)
        self.assertEqual(text.count("25RR-001"), 2)
        self.assertNotIn("25RR-002", text)

    def test_no_match_returns_empty_and_does_not_write_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            src = os.path.join(tmpdir, "SampleSheet_I10.csv")
            dst = os.path.join(tmpdir, "filtered", "SampleSheet_I10.csv")
            with open(src, "w") as f:
                f.write(SAMPLE_SHEET_V2)
            matched = write_filtered_samplesheet(src, dst, ["NONEXISTENT"])
            self.assertEqual(matched, set())
            self.assertFalse(os.path.exists(dst))

    def test_write_filtered_samplesheet_creates_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            src = os.path.join(tmpdir, "SampleSheet_U5N2_I10.csv")
            dst = os.path.join(tmpdir, "filtered", "SampleSheet_U5N2_I10.csv")
            with open(src, "w") as f:
                f.write(SAMPLE_SHEET_V2)
            matched = write_filtered_samplesheet(src, dst, ["25RR-001", "25RR-Q01"])
            self.assertEqual(matched, {"25RR-001", "25RR-Q01"})
            self.assertTrue(os.path.exists(dst))
            with open(dst) as f:
                content = f.read()
            self.assertIn("25RR-001", content)
            self.assertIn("25RR-Q01", content)
            self.assertNotIn("25RR-002", content)


if __name__ == "__main__":
    unittest.main()
