import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from research.scripts.compare_sc_author import compare_tables


ROOT = Path(__file__).resolve().parents[2]


class ScAuthorCompareTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.contrast = self.root / "contrast.tsv"
        self.author = self.root / "author.tsv"
        self.features = self.root / "features.tsv"
        self.output = self.root / "out/comparison.tsv"
        contrast_rows = [
            self._contrast("gA", "A", "Fibroblasts", "psoriasis_vs_normal", 0.001, -0.5, 0.2),
            self._contrast("gR", "RIGI", "Keratinocytes 1", "psoriasis_vs_normal", 0.001, 0.4, 0.01),
            self._contrast("gD", "DDX58", "Fibroblasts", "psoriasis_vs_normal", 0.001, 0.5, 0.01),
            self._contrast("gR2", "RIGI", "Fibroblasts", "psoriasis_vs_normal", 0.001, 0.4, 0.01),
            self._contrast("gM1", "MULTI", "Fibroblasts", "psoriasis_vs_normal", 0.001, 0.5, 0.01),
            self._contrast("gM2", "MULTI", "Fibroblasts", "psoriasis_vs_normal", 0.001, 0.6, 0.01),
            self._contrast("gK", "KRT", "basal", "psoriasis_vs_normal", 0.001, 0.5, 0.01),
            self._contrast("gA", "A", "Fibroblasts", "another_comparison", 0.001, 0.5, 0.01),
        ]
        author_rows = [
            self._author("A", "Fibroblasts", 0.5, 2, 1, 0.01),
            self._author("DDX58", "Keratinocytes 1", 0.5, 2, 1, 0.01),
            self._author("DDX58", "Fibroblasts", 0.5, 2, 1, 0.01),
            self._author("MISSING", "Fibroblasts", 1.0, 4, 2, 0.2),
            self._author("MULTI", "Fibroblasts", 0.5, 2, 1, 0.01),
            self._author("KRT", "Keratinocytes 2", 0.5, 2, 1, 0.01),
        ]
        feature_rows = [
            {"gene_id": "gA", "symbol": "A", "mapping_status": "unique"},
            {"gene_id": "gR", "symbol": "RIGI", "mapping_status": "unique"},
            {"gene_id": "gD", "symbol": "DDX58", "mapping_status": "unique"},
            {"gene_id": "gR2", "symbol": "RIGI", "mapping_status": "unique"},
            {"gene_id": "gM1", "symbol": "MULTI", "mapping_status": "unique"},
            {"gene_id": "gM2", "symbol": "MULTI", "mapping_status": "unique"},
            {"gene_id": "gK", "symbol": "KRT", "mapping_status": "unique"},
        ]
        self._write_tsv(self.contrast, contrast_rows)
        self._write_tsv(self.author, author_rows)
        self._write_tsv(self.features, feature_rows)

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def _contrast(gene_id, symbol, cell_type, comparison, p_value, log2fc, fdr):
        return {
            "gene_id": gene_id,
            "symbol": symbol,
            "cell_type": cell_type,
            "comparison": comparison,
            "p_value": str(p_value),
            "log2FC": str(log2fc),
            "fdr_bh": str(fdr),
        }

    @staticmethod
    def _author(gene, cell_type, logfc, psoriasis_mean, normal_mean, adjusted_p):
        return {
            "Gene": gene,
            "cell type": cell_type,
            "logFC": str(logfc),
            "p.value": "0.01",
            "Psoriasis _Mean": str(psoriasis_mean),
            "Normal_Mean": str(normal_mean),
            "adjusted_pvalue": str(adjusted_p),
            "qval": str(adjusted_p),
        }

    @staticmethod
    def _write_tsv(path, rows):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_exact_alias_missing_multiple_and_exact_cell_labels(self):
        rows = compare_tables(self.contrast, self.author, self.features, self.output)
        self.assertEqual(len(rows), 6)
        self.assertEqual(rows[0]["comparison_status"], "unique_match")
        self.assertEqual(rows[0]["actual_gene_id"], "gA")
        self.assertEqual(rows[0]["actual_comparison"], "psoriasis_vs_normal")
        self.assertEqual(rows[0]["direction_agrees_with_author"], "False")
        self.assertEqual(rows[0]["actual_passes_reported_threshold"], "False")
        self.assertEqual(rows[0]["threshold_agrees_with_author"], "False")
        self.assertEqual(rows[0]["formula_diagnostic_status"], "formula_diagnostic")
        self.assertAlmostEqual(float(rows[0]["formula_log2FC"]), 1.0)
        self.assertEqual(rows[1]["match_method"], "explicit_alias")
        self.assertEqual(rows[1]["actual_gene_id"], "gR")
        self.assertEqual(rows[2]["match_method"], "exact_symbol")
        self.assertEqual(rows[2]["actual_gene_id"], "gD")
        self.assertEqual(rows[3]["comparison_status"], "missing_result")
        self.assertEqual(rows[3]["actual_passes_reported_threshold"], "NA")
        self.assertEqual(rows[3]["direction_agrees_with_author"], "NA")
        self.assertEqual(rows[4]["comparison_status"], "multiple_matches")
        self.assertEqual(rows[5]["comparison_status"], "no_exact_cell_label")
        self.assertEqual(rows[5]["actual_gene_id"], "")

        with self.output.open(encoding="utf-8", newline="") as handle:
            written = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(written[0]["author_cell type"], "Fibroblasts")
        self.assertEqual(written[5]["author_cell type"], "Keratinocytes 2")

    def test_cli_accepts_all_input_and_output_paths(self):
        script = ROOT / "research/scripts/compare_sc_author.py"
        process = subprocess.run(
            [sys.executable, str(script), "--contrast", str(self.contrast), "--author", str(self.author),
             "--features", str(self.features), "--output", str(self.output)],
            text=True, capture_output=True, check=False,
        )
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(json.loads(process.stdout)["author_rows"], 6)
        self.assertTrue(self.output.is_file())


if __name__ == "__main__":
    unittest.main()
