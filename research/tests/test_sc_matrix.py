import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/scripts"))
from prepare_sc_matrix import prepare  # noqa: E402


class PrepareScMatrixTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        qc_dir = self.root / "research/results/sc_qc"
        processed = self.root / "research/data/processed/sc_qc"
        metadata = self.root / "research/data/metadata"
        qc_dir.mkdir(parents=True)
        processed.mkdir(parents=True)
        metadata.mkdir(parents=True)

        sample_rows = [
            {"gsm": "GSM_B", "donor_id": "D_B", "tissue_type": "Normal skin", "retained_cells": "2"},
            {"gsm": "GSM_A", "donor_id": "D_A", "tissue_type": "Psoriatic skin", "retained_cells": "2"},
        ]
        self._write_tsv(qc_dir / "sample_qc.tsv", sample_rows)
        cell_rows = [
            self._cell("GSM_B", "D_B", "Normal skin", "b1", 2, 2, True),
            self._cell("GSM_B", "D_B", "Normal skin", "b_drop", 9, 100, False),
            self._cell("GSM_B", "D_B", "Normal skin", "b2", 1, 2, True),
            self._cell("GSM_A", "D_A", "Psoriatic skin", "a1", 1, 1, True),
            self._cell("GSM_A", "D_A", "Psoriatic skin", "a2", 3, 5, True),
        ]
        with gzip.open(qc_dir / "cell_qc.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=cell_rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(cell_rows)

        self.genes = ["ENSG1", "ENSG2", "ENSGZERO", "ENSGAMB"]
        matrices = {
            "GSM_B": np.array([[1, 0], [0, 2], [0, 0], [1, 0]], dtype=np.int32),
            "GSM_A": np.array([[0, 3], [1, 1], [0, 0], [0, 1]], dtype=np.int32),
        }
        barcodes = {"GSM_B": ["b1", "b2"], "GSM_A": ["a1", "a2"]}
        for gsm, matrix in matrices.items():
            sparse.save_npz(processed / f"{gsm}.npz", sparse.csr_matrix(matrix))
            (processed / f"{gsm}.genes.txt").write_text("\n".join(self.genes) + "\n", encoding="utf-8")
            (processed / f"{gsm}.cells.txt").write_text("\n".join(barcodes[gsm]) + "\n", encoding="utf-8")

        gene_info = metadata / "Homo_sapiens.gene_info.gz"
        with gzip.open(gene_info, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("#tax_id", "GeneID", "Symbol", "dbXrefs"))
            writer.writerow(("9606", "1", "SAME", "Ensembl:ENSG1"))
            writer.writerow(("9606", "2", "SAME", "Ensembl:ENSG2"))
            writer.writerow(("9606", "3", "ZERO", "Ensembl:ENSGZERO"))
            writer.writerow(("9606", "4", "ONE", "Ensembl:ENSGAMB"))
            writer.writerow(("9606", "5", "TWO", "Ensembl:ENSGAMB"))

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def _cell(gsm, donor, tissue, barcode, n_features, n_counts, keep):
        return {
            "gsm": gsm,
            "donor_id": donor,
            "tissue_type": tissue,
            "barcode": barcode,
            "n_features": str(n_features),
            "n_counts": str(n_counts),
            "percent_mito": "1.5",
            "keep": str(keep),
        }

    @staticmethod
    def _write_tsv(path, rows):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_merge_preserves_order_counts_metadata_and_unique_symbol_status(self):
        result = prepare(self.root)
        output = self.root / "research/data/processed/sc_analysis"
        self.assertEqual(result["matrix_shape"], [3, 4])
        self.assertEqual(result["global_zero_genes_removed"], 1)
        with gzip.open(output / "counts.mtx.gz", "rb") as handle:
            actual = mmread(handle).toarray()
        expected = np.array([[1, 0, 0, 3], [0, 2, 1, 1], [1, 0, 0, 1]])
        np.testing.assert_array_equal(actual, expected)

        with (output / "features.tsv").open(encoding="utf-8", newline="") as handle:
            features = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["gene_id"] for row in features], ["ENSG1", "ENSG2", "ENSGAMB"])
        self.assertEqual([row["symbol"] for row in features], ["SAME", "SAME", ""])
        self.assertEqual([row["mapping_status"] for row in features], ["unique", "unique", "ambiguous"])

        with (output / "cells.tsv").open(encoding="utf-8", newline="") as handle:
            cells = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["cell_id"] for row in cells], ["GSM_B_b1", "GSM_B_b2", "GSM_A_a1", "GSM_A_a2"])
        self.assertEqual([row["n_counts"] for row in cells], ["2", "2", "1", "5"])
        run = json.loads((output / "run.json").read_text(encoding="utf-8"))
        self.assertEqual(len(run["input_sha256"]), 9)
        self.assertEqual(set(run["output_sha256"]), {"counts.mtx.gz", "features.tsv", "cells.tsv"})

    def test_count_mismatch_fails_with_cell_context(self):
        qc_path = self.root / "research/results/sc_qc/cell_qc.tsv.gz"
        with gzip.open(qc_path, "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows[0]["n_counts"] = "999"
        with gzip.open(qc_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(ValueError, "GSM_B/b1.*n_counts"):
            prepare(self.root)

    def test_gene_order_mismatch_fails_before_hstack(self):
        gene_path = self.root / "research/data/processed/sc_qc/GSM_A.genes.txt"
        gene_path.write_text("\n".join(reversed(self.genes)) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "GSM_A: gene order differs"):
            prepare(self.root)


if __name__ == "__main__":
    unittest.main()
