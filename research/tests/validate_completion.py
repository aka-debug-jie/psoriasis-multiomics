from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main():
    checks = []
    out = ROOT / "research/results/bulk_ncbi"
    run = json.loads((out / "run.json").read_text())
    assert run["result"]["status"] == "same_study_reprocessed_bulk_completed"
    for p, expected in run["inputs"].items(): assert sha(Path(p)) == expected, p
    for p, expected in run["outputs"].items(): assert sha(out / p) == expected, p
    tab = pd.read_csv(out / "all_gene_differential_expression.tsv", sep="\t")
    np.testing.assert_allclose(false_discovery_control(tab["P.Value"].to_numpy()), tab["adj.P.Val"], rtol=1e-8, atol=1e-10)
    assert len(tab) == 37065 and int(tab.significant.sum()) == 4266
    targets = pd.read_csv(out / "author_target_comparison.tsv", sep="\t")
    assert targets.GeneID.is_unique and len(targets) == 29
    reported = pd.read_csv(ROOT / "research/references/author_main_table1.tsv", sep="\t")
    expected = set(reported.loc[reported.bulk == "↑", "target"])
    found = set(targets.loc[targets.significant, "paper_target"])
    assert found == expected and len(found) == 13
    assert (targets.loc[targets.significant, "logFC"] > 0).all()
    checks.append("Bulk: integer-count source, 171 samples, 37065-gene BH, all 13 reported target directions and hashes checked")
    a = ROOT / "research/results/public_alternatives"
    r = json.loads((a / "run.json").read_text())
    assert r["status"] == "public_summary_lookup_completed"
    assert r["individual_PGS_or_Cox_performed"] is False
    assert r["replacement_for_original_PGS_score_slots"] is False
    assert sha(ROOT / "research/scripts/alternative_gwas_lookup.py") == r["script_sha256"]
    for p, expected in r["input_hashes"].items(): assert sha(ROOT / p) == expected, p
    for p, expected in r["output_hashes"].items(): assert sha(ROOT / p) == expected, p
    coords = json.loads((ROOT / "research/records/coordinate_validation.json").read_text())
    assert coords["checked"] == 3 and all(x["matches_Ensembl_GRCh37"] for x in coords["checks"])
    for name in ["japan_psoriasis", "european_psa"]:
        data = pd.read_csv(a / f"{name}_variant_comparison.tsv", sep="\t")
        assert data.rsid.is_unique
        assert len(data) == r[name]["matched_variants"]
        assert ((data.se_alt > 0) & np.isfinite(data.beta_alt_aligned_to_mvp)).all()
    checks.append("Alternative GWAS: unique effects, hashes, 3 independent coordinate probes; PGS/Cox slots not populated")
    scores = ROOT / "research/results/author_score_audit"
    s = json.loads((scores / "summary.json").read_text())
    assert s["matched_rows"] == 29 and s["mismatched_rows"] == 0
    assert sha(scores / "author_score_slots.tsv") == s["output_sha256"]
    assert sha(ROOT / "research/scripts/audit_author_scores.R") == s["script_sha256"]
    checks.append("Author scores: all 29 match nine-slot arithmetic; kept separate from reconstructed evidence")
    (ROOT / "research/records/completion_validation.json").write_text(json.dumps({"status": "passed", "checks": checks}, indent=2) + "\n")
    print("\n".join("PASS " + x for x in checks))


if __name__ == "__main__":
    main()
