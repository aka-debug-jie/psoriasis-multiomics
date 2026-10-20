"""Identify numerical agreement of S2 outcomes with acquired GWAS data.

Agreement is evidence about the table, not a declaration of author intent.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main():
    source = ROOT / "research/references/author_table_S2.tsv"
    reported = pd.read_csv(source, sep="\t")
    inputs = [source]
    comparisons = []
    for outcome in ["mvp", "finngen_r12", "stuart"]:
        path = ROOT / "research/data/processed/outcomes" / f"{outcome}_normalized.tsv"
        inputs.append(path)
        actual = pd.read_csv(path, sep="\t")
        for i, r in reported.iterrows():
            matches = []
            for _, a in actual[actual.rsid == r.SNP].iterrows():
                if (r.EA, r.NEA) == (a.EA, a.OA):
                    matches.append((a.beta, a.se))
                elif (r.EA, r.NEA) == (a.OA, a.EA):
                    matches.append((-a.beta, a.se))
            matches = list(set(matches))
            row = {"s2_row": i + 1, "protein": r.Protein, "exposure_cohort": r.dat, "rsid": r.SNP,
                   "outcome": outcome, "matching_effects": len(matches),
                   "reported_beta": r["beta.outcome"], "reported_se": r["se.outcome"]}
            if len(matches) == 1:
                b, se = matches[0]
                row.update(original_beta_aligned=b, original_se=se,
                           beta_error=b-r["beta.outcome"], se_error=se-r["se.outcome"],
                           agrees_at_four_decimal_places=bool(abs(b-r["beta.outcome"]) <= .00005001
                                                               and abs(se-r["se.outcome"]) <= .00005001))
            comparisons.append(row)
    data = pd.DataFrame(comparisons)
    output = ROOT / "research/results/supplement_outcome_audit"
    output.mkdir(parents=True, exist_ok=True)
    target = output / "comparison.tsv"
    data.to_csv(target, sep="\t", index=False)
    summary = []
    for outcome, g in data.groupby("outcome"):
        summary.append({"outcome": outcome, "s2_rows": len(g),
                        "uniquely_matched": int((g.matching_effects == 1).sum()),
                        "beta_and_se_agree_at_4dp": int(g.agrees_at_four_decimal_places.eq(True).sum())})
    result = {"status": "cross_dataset_numeric_audit", "outcomes": summary,
              "interpretation": "Numerical agreement only; S2 does not name its outcome dataset",
              "input_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
              "output_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (output / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
