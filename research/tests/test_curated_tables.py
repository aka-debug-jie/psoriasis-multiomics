import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "curate_release_tables.py"
SPEC = importlib.util.spec_from_file_location("curate_release_tables", SCRIPT)
curator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(curator)


SPECIAL = ["mr_mvp_cohorts.tsv", "mr_finngen_cohorts.tsv", "mr_stuart_cohorts.tsv"]
EXPECTED_HEADERS = ["rsid", "exposure_neg_log10_p", "mr_neg_log10_p", "F_exposure_source", "F_exposure_recomputed", "note"]


def make_fixture(root, unexpected_duplicate=False):
    data_root = root / "research" / "release" / "data"
    data_root.mkdir(parents=True)
    file_names = SPECIAL + [f"other_{i:02}.tsv" for i in range(27)]
    manifest = {"files": [{"output": name} for name in file_names]}
    (data_root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    special_bytes = b'rsid\tneg_log10_p\tneg_log10_p\tF_statistic\tF_statistic\tnote\nr1\tNA\t\t2.0\t3.0\t"line one\nline two"\n'
    for name in file_names:
        path = data_root / name
        if name in SPECIAL:
            path.write_bytes(special_bytes)
        elif unexpected_duplicate and name == file_names[3]:
            path.write_bytes(b"id\tid\nleft\tright\n")
        else:
            path.write_bytes(b"id\tvalue\nx\tNA\n")
    return data_root, file_names, special_bytes


class CuratedReleaseTablesTest(unittest.TestCase):
    def test_renames_only_approved_headers_and_preserves_cells(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data_root, file_names, special_bytes = make_fixture(root)
            result = subprocess.run([sys.executable, str(SCRIPT), "--root", str(root)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)

            output_root = root / "research" / "release" / "tables"
            for name in file_names:
                source = data_root / name
                output = output_root / name
                if name in SPECIAL:
                    self.assertNotEqual(output.read_bytes(), special_bytes)
                    with source.open(encoding="utf-8", newline="") as handle:
                        source_rows = list(csv.reader(handle, delimiter="\t"))
                    with output.open(encoding="utf-8", newline="") as handle:
                        output_rows = list(csv.reader(handle, delimiter="\t"))
                    self.assertEqual(output_rows[0], EXPECTED_HEADERS)
                    self.assertEqual(output_rows[1:], source_rows[1:])
                    self.assertEqual(output_rows[1][1:5], ["NA", "", "2.0", "3.0"])
                else:
                    self.assertEqual(output.read_bytes(), source.read_bytes())

            with (output_root / "field_mapping.tsv").open(encoding="utf-8", newline="") as handle:
                mappings = list(csv.DictReader(handle, delimiter="\t"))
            mvp_map = [r for r in mappings if r["table"] == SPECIAL[0]]
            self.assertEqual([r["new_name"] for r in mvp_map], EXPECTED_HEADERS)
            self.assertEqual([r["index"] for r in mvp_map], ["1", "2", "3", "4", "5", "6"])
            manifest = json.loads((output_root / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(len(manifest["files"]), 30)
            self.assertEqual(hashlib.sha256((output_root / "field_mapping.tsv").read_bytes()).hexdigest(), manifest["mapping"]["sha256"])

    def test_unapproved_duplicate_header_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            make_fixture(root, unexpected_duplicate=True)
            result = subprocess.run([sys.executable, str(SCRIPT), "--root", str(root)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Unexpected duplicate header", result.stderr)


if __name__ == "__main__":
    unittest.main()
