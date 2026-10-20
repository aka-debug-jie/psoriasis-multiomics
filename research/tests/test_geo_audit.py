import tempfile
import unittest
from pathlib import Path

from research.scripts.audit_geo import _record, _validate_processing, sample_blocks


class GeoAuditParserTests(unittest.TestCase):
    def test_repeated_characteristics_and_inclusion_rules(self):
        soft = """^SAMPLE = GSM1
!Sample_title = healthy
!Sample_characteristics_ch1 = subject id: donor D1
!Sample_characteristics_ch1 = subject status: healthy donor
!Sample_characteristics_ch1 = tissue type: Normal skin
!Sample_supplementary_file_1 = ftp://example/sample.csv.gz
^SAMPLE = GSM2
!Sample_title = peripheral
!Sample_characteristics_ch1 = subject id: donor D2
!Sample_characteristics_ch1 = subject status: psoriasis patient
!Sample_characteristics_ch1 = tissue type: Peripheral normal skin
"""
        with tempfile.TemporaryDirectory() as temporary_dir:
            soft_path = Path(temporary_dir) / "small.soft"
            soft_path.write_text(soft, encoding="utf-8")
            blocks = sample_blocks(soft_path.read_text(encoding="utf-8"))
        self.assertEqual(len(blocks), 2)
        healthy = _record(blocks[0], sc=True)
        peripheral = _record(blocks[1], sc=True)
        self.assertTrue(healthy["selected"])
        self.assertEqual(healthy["donor_id"], "D1")
        self.assertEqual(healthy["supplementary_url"], "ftp://example/sample.csv.gz")
        self.assertFalse(peripheral["selected"])

    def test_missing_donor_and_empty_supplement(self):
        soft = """^SAMPLE = GSM3
!Sample_title = no donor metadata
!Sample_characteristics_ch1 = tissue type: Normal skin
!Sample_supplementary_file_1 = NONE
"""
        with tempfile.TemporaryDirectory() as temporary_dir:
            soft_path = Path(temporary_dir) / "missing_donor.soft"
            soft_path.write_text(soft, encoding="utf-8")
            blocks = sample_blocks(soft_path.read_text(encoding="utf-8"))
        row = _record(blocks[0], sc=True)
        self.assertEqual(row["donor_id"], "")
        self.assertEqual(row["supplementary_url"], "")
        self.assertFalse(row["selected"])

    def test_processing_unit_mismatch_names_sample(self):
        blocks = sample_blocks("""^SAMPLE = GSM404
!Sample_data_processing = normalized RPKM matrix
""")
        with self.assertRaisesRegex(ValueError, "GSM404.*raw gene counts"):
            _validate_processing(blocks, unit="raw gene counts")


if __name__ == "__main__":
    unittest.main()
