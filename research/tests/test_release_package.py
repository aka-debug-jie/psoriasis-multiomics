import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "package_release_data.py"
SPEC = importlib.util.spec_from_file_location("package_release_data", SCRIPT)
package_release_data = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(package_release_data)


class ReleasePackageTest(unittest.TestCase):
    def test_preserves_bytes_and_records_table_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            research = root / "research"
            sample = b'gene\tvalue\tnote\nG1\tNA\t\nG2\t2.5\t"multi\nline"\n'
            for index, source_rel in enumerate(package_release_data.FILES):
                source = research / source_rel
                source.parent.mkdir(parents=True, exist_ok=True)
                source.write_bytes(sample if index == 0 else b"id\tvalue\nX\t1\n")

            result = subprocess.run([sys.executable, str(SCRIPT), "--root", str(root)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)

            output = research / "release" / "data"
            first_target = output / next(iter(package_release_data.FILES.values()))
            self.assertEqual(first_target.read_bytes(), sample)
            for source_rel, output_rel in package_release_data.FILES.items():
                self.assertEqual((output / output_rel).read_bytes(), (research / source_rel).read_bytes(), output_rel)
            with (output / "catalog.tsv").open(encoding="utf-8", newline="") as handle:
                catalog = list(csv.DictReader(handle, delimiter="\t"))
            first = catalog[0]
            self.assertEqual(first["rows"], "2")
            self.assertEqual(first["columns"], "3")
            self.assertEqual(json.loads(first["column_names"]), ["gene", "value", "note"])
            digest = hashlib.sha256(sample).hexdigest()
            self.assertEqual(first["source_sha256"], digest)
            self.assertEqual(first["output_sha256"], digest)
            with (output / "dictionary.tsv").open(encoding="utf-8", newline="") as handle:
                dictionary = list(csv.DictReader(handle, delimiter="\t"))
            first_columns = [row for row in dictionary if row["output"] == first_target.name]
            self.assertEqual(first_columns[0]["column_position"], "1")
            self.assertEqual(first_columns[1]["na_count"], "1")
            self.assertEqual(first_columns[2]["empty_count"], "1")
            self.assertEqual(first_columns[1]["dtype_observed"], "numeric")


if __name__ == "__main__":
    unittest.main()
