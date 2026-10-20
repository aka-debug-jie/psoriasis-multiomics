import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri_exp

from research.scripts.prepare_sensitivity_inputs import prepare


ROOT = Path(__file__).resolve().parents[2]


class SensitivityInputPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / "research/data/processed/exposures").mkdir(parents=True)
        (self.root / "research/data/processed/smr_gwas").mkdir(parents=True)
        (self.root / "research/results/available_followup").mkdir(parents=True)
        (self.root / "research/results/mr_available").mkdir(parents=True)

        self._write_tsv(self.root / "research/results/available_followup/followup_proteins.tsv", [
            {"uniprot": "P1", "gene": "G1"},
            {"uniprot": "P2", "gene": "G2"},
        ])
        self._write_tsv(self.root / "research/results/mr_available/cohort_wald.tsv", [
            self._wald("Iceland", "A1", "G1", "P1", "rs1,rs11", 1, 100, "A", "G", 0.25, 0.2, 2.0, "101"),
            self._wald("Iceland", "A2", "G2", "P2", "rs999", 2, 200, "C", "T", 0.4, -0.3, 3.2, "202"),
            self._wald("UKB-PPP", "B1", "G1", "P1", "rs2", 1, 300, "G", "A", 0.15, 0.4, 4.0, "303"),
        ])
        self._write_tsv(self.root / "research/data/processed/exposures/ST19_somascan_pQTLs_smpnorm.tsv", [
            self._st19("P1", "A1", "G1", "rs1,rs11", "cis", "chr1", 100, "A", "G", 25, 0.2, 2.0),
            self._st19("P1", "A1", "G1", "chr1:101", "trans", "chr1", 101, "C", "A", 40, 0.01, 0.1),
            self._st19("P1", "OTHER", "G1", "rsNotSelected", "cis", "chr1", 102, "T", "C", 30, 2.0, 10.0),
            self._st19("P2", "A2", "G2", "rs999", "trans", "chr2", 200, "C", "T", 40, -0.3, 3.2),
            self._st19("P3", "A3", "G3", "rsUnfollowed", "cis", "chr3", 300, "A", "C", 20, 0.5, 5.0),
        ])
        self._write_tsv(self.root / "research/data/processed/exposures/ST21_olink_ukb_bi_pQTLs.tsv", [
            self._st21("P1", "B1", "G1", "rs2", "cis", "chr1", 300, "G", "A", 15, 0.4, 4.0),
            self._st21("P1", "OTHER", "G1", "rsNotSelected2", "trans", "chr1", 301, "T", "C", 20, 1.5, 8.0),
        ])
        self._write_tsv(self.root / "research/data/processed/smr_gwas/mvp_full.ma", [
            {"SNP": "rs1", "A1": "A", "A2": "G", "freq": "0.1", "b": "0.2", "se": "0.1", "p": "0.05", "n": "10"},
            {"SNP": "rs11", "A1": "A", "A2": "G", "freq": "0.2", "b": "0.3", "se": "0.1", "p": "0.04", "n": "10"},
            {"SNP": "rs1", "A1": "A", "A2": "G", "freq": "0.3", "b": "0.4", "se": "0.1", "p": "0.03", "n": "10"},
            {"SNP": "rs2", "A1": "G", "A2": "A", "freq": "0.4", "b": "0.5", "se": "0.1", "p": "0.02", "n": "10"},
            {"SNP": "rsUnrelated", "A1": "C", "A2": "T", "freq": "0.5", "b": "0.6", "se": "0.1", "p": "0.01", "n": "10"},
        ])

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def _wald(cohort, assay, gene, uniprot, variant, chrom, pos, ea, oa, eaf, beta, logp, row_id):
        return {
            "cohort": cohort, "assay": assay, "gene": gene, "uniprot": uniprot, "variant": variant,
            "chr.exposure": str(chrom), "pos.exposure": str(pos), "EA.exposure": ea, "OA.exposure": oa,
            "eaf.exposure": str(eaf), "beta.exposure": str(beta), "neg_log10_p": str(logp),
            "exposure_row_id": row_id,
        }

    @staticmethod
    def _st19(uniprot, assay, gene, variant, locus_type, chrom, pos, ea, oa, freq, beta, logp):
        return {
            "UniProt": uniprot, "SeqId": assay, "gene prot": gene, "variant": variant,
            "cis trans": locus_type, "chr var": chrom, "pos var": pos, "Amin": ea, "Amaj": oa,
            "MAF PC": freq, "beta unadj": beta, "mLog10pval unadj": logp,
        }

    @staticmethod
    def _st21(uniprot, assay, gene, variant, locus_type, chrom, pos, ea, oa, freq, beta, logp):
        return {
            "UniProt": uniprot, "OlinkID": assay, "Affected Protein Gene name": gene, "rsName": variant,
            "CisOrTrans": locus_type, "Chrom": chrom, "Pos": pos, "Amin": ea, "Amaj": oa,
            "MAF PC": freq, "Effect Amin unadj": beta, "Log10 pval gc cor unadj": logp,
        }

    @staticmethod
    def _write_tsv(path, rows):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_selected_assays_all_loci_and_streamed_mvp_aliases(self):
        result = prepare(self.root)
        out = self.root / "research/data/processed/sensitivity_inputs"
        loci = pd.read_csv(out / "all_loci.tsv", sep="\t")
        self.assertEqual(len(loci), 4)
        self.assertEqual(set(loci["cohort"]), {"Iceland", "UKB-PPP"})
        self.assertEqual(set(loci["locus_type"]), {"cis", "trans"})
        self.assertNotIn("OTHER", set(loci["assay"]))
        self.assertNotIn("rsUnfollowed", set(loci["variant"]))
        low_p_row = loci.loc[loci["variant"] == "chr1:101"].iloc[0]
        self.assertEqual(low_p_row["neg_log10_p"], 0.1)
        self.assertEqual(low_p_row["source_data_row"], 2)
        self.assertEqual(loci.loc[loci["variant"] == "rs1,rs11", "eaf"].iloc[0], 0.25)
        expected_se = abs(0.2) / -ndtri_exp(-2.0 * np.log(10) - np.log(2))
        self.assertAlmostEqual(loci.loc[loci["variant"] == "rs1,rs11", "se_reconstructed"].iloc[0], expected_se)
        self.assertIn("not_reported_SE", loci["se_source"].iloc[0])

        with (out / "mvp_subset.tsv").open(encoding="utf-8", newline="") as handle:
            mvp = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["SNP"] for row in mvp], ["rs1", "rs11", "rs1", "rs2"])
        self.assertEqual(result["mvp_subset"]["missing_rsids"], ["rs999"])
        self.assertEqual(result["mvp_subset"]["matched_rows"], 4)
        self.assertEqual(result["all_loci_rows"], 4)
        self.assertEqual(len(result["input_sha256"]), 5)
        self.assertEqual(set(result["output_sha256"]), {"all_loci.tsv", "mvp_subset.tsv"})
        run = json.loads((self.root / "research/records/advance_20261017/sensitivity_inputs_run.json").read_text(encoding="utf-8"))
        self.assertIn("no p/F/LD/direction filtering", run["scope"])


if __name__ == "__main__":
    unittest.main()
