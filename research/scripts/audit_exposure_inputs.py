"""Compare published exposure inputs with their original pQTL worksheets.

The inversion of two-sided p-values is checked against Table S2. It reconstructs
the table's numerical convention, not an independently reported standard error.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.special import ndtri_exp

ROOT = Path(__file__).resolve().parents[2]


def main():
    s2_path = ROOT / "research/references/author_table_S2.tsv"
    s2 = pd.read_csv(s2_path, sep="\t")
    output = ROOT / "research/results/exposure_audit"
    output.mkdir(parents=True, exist_ok=True)
    comparisons, inventories, inputs, normalized = [], [], [s2_path], []
    specifications = [
        ("ST19_somascan_pQTLs_smpnorm", "Iceland", "variant", "gene prot", "cis trans", "beta unadj",
         "mLog10pval unadj", "chr var", "pos var", "SeqId", "MAF PC"),
        ("ST21_olink_ukb_bi_pQTLs", "UKB-PPP", "rsName", "Affected Protein Gene name", "CisOrTrans",
         "Effect Amin unadj", "Log10 pval gc cor unadj", "Chrom", "Pos", "OlinkID", "MAF PC")]
    for sheet, cohort, rsid, gene, cis, beta, logp, chrom, pos, assay, maf in specifications:
        path = ROOT / "research/data/processed/exposures" / f"{sheet}.tsv"
        d = pd.read_csv(path, sep="\t", low_memory=False)
        inputs.append(path)
        c = d[d[cis] == "cis"].copy()
        strongest = c.sort_values(logp, ascending=False).drop_duplicates(gene).set_index(gene)
        inventories.append({"cohort": cohort, "all_pqtl_rows": len(d), "cis_rows": len(c),
                            "cis_gene_labels": c[gene].nunique(), "cis_assays": c[assay].nunique()})
        # Preserve every original cis row, including unresolved indels and markers.
        # No deduplication or selection is performed in this normalized input table.
        n = pd.DataFrame({"cohort": cohort, "assay": c[assay], "gene": c[gene], "uniprot": c["UniProt"],
                          "variant": c[rsid], "chr": c[chrom].astype(str).str.replace("chr", "", regex=False),
                          "pos": c[pos], "EA": c["Amin"], "OA": c["Amaj"], "eaf": c[maf] / 100,
                          "beta": c[beta], "neg_log10_p": c[logp], "build": "GRCh38"})
        n["se_reconstructed"] = np.abs(n["beta"]) / -ndtri_exp(-n["neg_log10_p"] * np.log(10) - np.log(2))
        n["se_source"] = "normal_inversion_of_original_LR_p_value;matches_author_S2_rounding"
        n["input_role"] = "unselected_cis_source_rows"
        normalized.append(n)
        for _, row in s2[s2.dat == cohort].iterrows():
            match = d[(d[rsid] == row.SNP) & (d[gene] == row.Protein)]
            if len(match) != 1:
                comparisons.append({"cohort": cohort, "protein": row.Protein, "SNP": row.SNP,
                                    "status": "nonunique_or_missing_source_row", "matches": len(match)})
                continue
            m = match.iloc[0]
            if (m.Amin, m.Amaj) == (row.EA, row.NEA):
                sign = 1
            elif (m.Amin, m.Amaj) == (row.NEA, row.EA):
                sign = -1
            else:
                raise ValueError(f"Allele mismatch in author source audit: {cohort} {row.Protein}")
            se = abs(m[beta]) / -ndtri_exp(-m[logp] * np.log(10) - np.log(2))
            comparisons.append({"cohort": cohort, "protein": row.Protein, "SNP": row.SNP,
                                "status": "matched_original_row", "source_beta_aligned": sign * m[beta],
                                "author_beta": row["beta.exposure"], "beta_error": sign * m[beta] - row["beta.exposure"],
                                "reconstructed_se": se, "author_se": row["se.exposure"],
                                "se_error": se - row["se.exposure"],
                                "same_as_strongest_cis_by_gene": strongest.loc[row.Protein, rsid] == row.SNP,
                                "interpretation": "SE numerical agreement only; index selection still requires outcome availability"})
    normalized_path = ROOT / "research/data/processed/exposures/eldjarn_all_cis_unselected.tsv"
    pd.concat(normalized, ignore_index=True).to_csv(normalized_path, sep="\t", index=False)
    pd.DataFrame(comparisons).to_csv(output / "author_exposure_source_comparison.tsv", sep="\t", index=False)
    result = {"status": "source_inputs_checked_not_discovery", "cohorts": inventories,
              "author_rows_compared": len(comparisons), "all_rows_matched": all(x["status"] == "matched_original_row" for x in comparisons),
              "max_absolute_beta_error": max(abs(x["beta_error"]) for x in comparisons if "beta_error" in x),
              "max_absolute_se_error": max(abs(x["se_error"]) for x in comparisons if "se_error" in x),
              "missing_cohort": "Fenland", "independent_source_of_se": False,
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "normalized_exposure_sha256": hashlib.sha256(normalized_path.read_bytes()).hexdigest(),
              "input_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (output / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
