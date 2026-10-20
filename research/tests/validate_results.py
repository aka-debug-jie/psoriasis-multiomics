"""Independent numerical checks and byte-level provenance verification."""
from pathlib import Path
import hashlib
import json
import math
import csv
from scipy.stats import norm
from scipy import sparse
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main():
    checks = []
    out = ROOT / "research/results/supplement_check"
    rows = list(csv.DictReader((out / "wald_by_cohort.tsv").open(), delimiter="\t"))
    for r in rows:
        b = float(r["beta.outcome"]) / float(r["beta.exposure"])
        se = float(r["se.outcome"]) / abs(float(r["beta.exposure"]))
        for actual, expected in [(float(r["beta"]), b), (float(r["se"]), se),
                                 (float(r["p"]), 2 * norm.sf(abs(b / se)))]:
            assert math.isclose(actual, expected, rel_tol=1e-8, abs_tol=1e-10)
    combined = list(csv.DictReader((out / "fixed_effect_by_protein.tsv").open(), delimiter="\t"))
    for r in combined:
        group = [x for x in rows if x["Protein"] == r["protein"]]
        weights = [1 / float(x["se"]) ** 2 for x in group]
        b = sum(w * float(x["beta"]) for w, x in zip(weights, group)) / sum(weights)
        assert math.isclose(float(r["beta"]), b, rel_tol=1e-8, abs_tol=1e-10)
        assert math.isclose(float(r["se"]), math.sqrt(1 / sum(weights)), rel_tol=1e-8, abs_tol=1e-10)
    run = json.loads((out / "run.json").read_text())
    for path, expected in run["inputs"].items():
        assert sha256(Path(path)) == expected, path
    for name, expected in run["outputs"].items():
        assert sha256(out / name) == expected, name
    checks.append(f"MR: {len(rows)} Wald rows and {len(combined)} fixed-effect rows independently checked; run hashes verified")
    geo = json.loads((ROOT / "research/results/geo_audit/summary.json").read_text())
    for name, expected in geo["input_sha256"].items():
        assert sha256(ROOT / "research/data/metadata" / name) == expected
    checks.append("GEO: original SOFT byte hashes verified")
    mr_path = ROOT / "research/results/mr_available"
    if (mr_path / "run.json").exists():
        mr_run = json.loads((mr_path / "run.json").read_text())
        assert mr_run["result"]["status"] == "partial_available_data_calculated"
        assert mr_run["result"]["formal_candidate_list_produced"] is False
        for name, expected in mr_run["inputs"].items():
            assert sha256(Path(name)) == expected, name
        for name, expected in mr_run["outputs"].items():
            assert sha256(mr_path / name) == expected, name
        instruments = pd.read_csv(mr_path / "cohort_wald.tsv", sep="\t")
        source_keys = set(zip(instruments.cohort, instruments.uniprot, instruments.rsid))
        for prefix, rows_name in [("", "cohort_wald.tsv"), ("finngen", "finngen_same_instrument_wald.tsv"),
                                  ("stuart", "stuart_same_instrument_wald.tsv")]:
            tab = pd.read_csv(mr_path / rows_name, sep="\t")
            assert not tab.duplicated(["cohort", "uniprot"]).any()
            assert set(zip(tab.cohort, tab.uniprot, tab.rsid)) <= source_keys
            np.testing.assert_allclose(tab.beta_mr, tab.beta_outcome_aligned / tab["beta.exposure"], rtol=1e-8, atol=1e-10)
            np.testing.assert_allclose(tab.se_mr, tab.se / np.abs(tab["beta.exposure"]), rtol=1e-8, atol=1e-10)
            meta_name = "available_meta.tsv" if not prefix else f"{prefix}_available_meta.tsv"
            meta = pd.read_csv(mr_path / meta_name, sep="\t").set_index("uniprot")
            for protein, group in tab.groupby("uniprot"):
                w = 1 / group.se_mr ** 2
                assert math.isclose(meta.loc[protein, "beta"], float((w * group.beta_mr).sum() / w.sum()), rel_tol=1e-8, abs_tol=1e-10)
                assert math.isclose(meta.loc[protein, "se"], math.sqrt(1 / w.sum()), rel_tol=1e-8, abs_tol=1e-10)
        checks.append("Available MR: three outcome calculations, fixed IVs, unique effects, meta weights and hashes independently checked")
    sc_path = ROOT / "research/results/sc_qc/run.json"
    if sc_path.exists():
        sc = json.loads(sc_path.read_text())
        assert sc["status"] == "qc_completed_not_single_cell_reproduction"
        assert sha256(ROOT / "research/scripts/sc_qc.py") == sc["script_sha256"]
        for item in sc["input_hashes"]:
            assert sha256(ROOT / item["path"]) == item["sha256"]
        for name, expected in sc["metadata_hashes"].items():
            assert sha256(ROOT / name) == expected
        for name, expected in sc["output_hashes"].items():
            assert sha256(ROOT / name) == expected, name
        qc = pd.read_csv(ROOT / "research/results/sc_qc/cell_qc.tsv.gz", sep="\t")
        for gsm, group in qc.groupby("gsm", sort=False):
            kept = group[group["keep"]]
            matrix = sparse.load_npz(ROOT / "research/data/processed/sc_qc" / f"{gsm}.npz")
            np.testing.assert_array_equal(np.asarray(matrix.sum(axis=0)).ravel(), kept["n_counts"].to_numpy())
            np.testing.assert_array_equal(matrix.getnnz(axis=0), kept["n_features"].to_numpy())
        assert len(qc) == sc["input_cells"] and int(qc["keep"].sum()) == sc["retained_cells"]
        checks.append("SC QC: code/output hashes, sparse matrix cell counts and feature counts independently checked")
    print("\n".join("PASS " + c for c in checks))
    (ROOT / "research/records/validation.json").write_text(json.dumps({"status": "passed", "checks": checks}, indent=2) + "\n")


if __name__ == "__main__":
    main()
