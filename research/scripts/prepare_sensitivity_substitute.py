#!/usr/bin/env python3
"""Prepare published loci for the separate Fenland-substitute sensitivity run."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "research/data/processed/sensitivity_substitute"
FOLLOWUP = ROOT / "research/results/mr_fenland_substitute/followup_proteins.tsv"
WALD = ROOT / "research/results/mr_fenland_substitute/mvp_cohort_wald.tsv"
INDICES = ROOT / "research/data/processed/Fenland_Olink/source_indices.tsv"
FENLAND_RAW = ROOT / "research/data/raw/Fenland_Olink/Summary.statistics.cis.Olink.proteins.Fenland.20220107.txt.gz"
MVP = ROOT / "research/data/processed/smr_gwas/mvp_full.ma"
SE_SOURCES = Path(__file__).resolve().parents[1] / "scripts/prepare_sensitivity_inputs.py"


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_old_helpers():
    spec = importlib.util.spec_from_file_location("prepare_sensitivity_inputs", SE_SOURCES)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def numeric_match(value, expected, field, context, rel=1e-10, abs_=1e-10):
    try:
        a, b = float(value), float(expected)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context}: nonnumeric {field}: {value!r}, {expected!r}") from exc
    if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a, b, rel_tol=rel, abs_tol=abs_):
        raise ValueError(f"{context}: {field} mismatch: {a!r} vs {b!r}")


def adapt_selected_wald(row):
    """Map merged MVP wald fields to the legacy helper's unsuffixed schema."""
    return {
        "cohort": str(row["cohort"]), "assay": str(row["assay"]), "gene": str(row["gene"]),
        "uniprot": str(row["uniprot"]), "variant": str(row["variant"]),
        "chr.exposure": str(row["chr"]), "pos.exposure": str(row["pos"]),
        "EA.exposure": str(row["EA_exposure"]), "OA.exposure": str(row["OA_exposure"]),
        "eaf.exposure": str(row["eaf_exposure"]), "beta.exposure": str(row["beta_exposure"]),
        "neg_log10_p": str(row["neg_log10_p"]), "exposure_row_id": str(row["source_row"]),
    }


def selected_inputs():
    follow = pd.read_csv(FOLLOWUP, sep="\t", dtype={"uniprot": str})
    if "uniprot" not in follow or follow["uniprot"].isna().any() or follow["uniprot"].astype(str).duplicated().any():
        raise ValueError("followup_proteins.tsv must contain 17 unique nonmissing UniProt IDs")
    followups = set(follow["uniprot"].astype(str))
    if len(followups) != 17:
        raise ValueError(f"expected 17 follow-up UniProt IDs, found {len(followups)}")
    wald = pd.read_csv(WALD, sep="\t", dtype=str, keep_default_na=False)
    required = {"cohort", "assay", "gene", "uniprot", "variant", "chr", "pos", "EA_exposure", "OA_exposure",
                "eaf_exposure", "beta_exposure", "se_exposure", "neg_log10_p", "source_row", "sdY_estimated",
                "raw_beta", "raw_se"}
    missing = sorted(required - set(wald.columns))
    if missing:
        raise ValueError(f"mvp_cohort_wald.tsv lacks required fields: {missing}")
    selected = wald[wald["uniprot"].isin(followups)].copy()
    if selected[["cohort", "uniprot", "assay"]].eq("").any().any():
        raise ValueError("selected mvp_cohort_wald row has empty cohort/uniprot/assay")
    counts = selected.groupby(["cohort", "uniprot"], dropna=False).size()
    bad = counts[counts != 1]
    if len(bad):
        raise ValueError(f"expected one frozen assay row per cohort/UniProt: {bad.to_dict()}")
    allowed = {"Iceland", "UKB-PPP", "Fenland_Olink_substitute_485"}
    unexpected = set(selected["cohort"]) - allowed
    if unexpected:
        raise ValueError(f"unexpected cohorts in selected wald rows: {sorted(unexpected)}")
    return followups, selected


def fenland_records(fen_rows: pd.DataFrame, indices: pd.DataFrame, raw_path: Path):
    # Each selected assay must be the exact frozen source index in the main wald table.
    idx = indices.set_index("id", drop=False)
    chosen = {}
    for _, row in fen_rows.iterrows():
        aid = str(row["assay"])
        if aid not in idx.index:
            raise ValueError(f"selected Fenland assay missing from source_indices.tsv: {aid}")
        source = idx.loc[aid]
        if isinstance(source, pd.DataFrame):
            raise ValueError(f"nonunique source_indices id: {aid}")
        context = f"Fenland/{row['uniprot']}/{aid}/{row['variant']}"
        if str(source["Uniprot_ID"]) != str(row["uniprot"]):
            raise ValueError(f"{context}: frozen index UniProt differs")
        for field, actual, expected in [
            ("variant", source["rsid"], row["variant"]), ("chr", source["chr"], row["chr"]),
            ("pos", source["pos"], row["pos"]), ("EA", source["Allele1"], row["EA_exposure"]),
            ("OA", source["Allele2"], row["OA_exposure"]),
        ]:
            if str(actual) != str(expected):
                raise ValueError(f"{context}: {field} mismatch: index={actual!r}, wald={expected!r}")
        for field, actual, expected in [
            ("eaf", source["Freq1"], row["eaf_exposure"]),
            ("beta", source["beta_estimated_SD"], row["beta_exposure"]),
            ("se", source["se_estimated_SD"], row["se_exposure"]),
            ("raw_beta", source["Effect"], row["raw_beta"]),
            ("raw_se", source["StdErr"], row["raw_se"]),
            ("sdY_estimated", source["sdY_estimated"], row["sdY_estimated"]),
        ]:
            numeric_match(actual, expected, field, context)
        chosen[aid] = {"index": source.to_dict(), "wald": row.to_dict()}
    if not chosen:
        return pd.DataFrame()

    selected_ids = set(chosen)
    required = {"id", "chr", "rsid", "pos", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "log10p"}
    out = []
    offset = 0
    for chunk in pd.read_csv(raw_path, sep="\t", chunksize=150000, low_memory=False, dtype={"id": str, "rsid": str}):
        miss = sorted(required - set(chunk.columns))
        if miss:
            raise ValueError(f"Fenland raw source missing fields: {miss}")
        selected = chunk[chunk["id"].isin(selected_ids)].copy()
        if len(selected):
            selected["_source_data_row"] = np.arange(offset + 1, offset + len(chunk) + 1, dtype=np.int64)[chunk["id"].isin(selected_ids).to_numpy()]
            for _, raw in selected.iterrows():
                aid = str(raw["id"])
                target = chosen[aid]
                idxrow, wald = target["index"], target["wald"]
                try:
                    beta_raw, se_raw = float(raw["Effect"]), float(raw["StdErr"])
                    freq, sd = float(raw["Freq1"]), float(idxrow["sdY_estimated"])
                    logp = float(raw["log10p"])
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Fenland {aid} raw row {raw['_source_data_row']}: nonnumeric scale fields") from exc
                if not all(math.isfinite(x) for x in (beta_raw, se_raw, freq, sd, logp)) or se_raw <= 0 or sd <= 0:
                    raise ValueError(f"Fenland {aid} raw row {raw['_source_data_row']}: invalid beta/SE/freq/sdY/log10p")
                # Preserve every source pair for the selected assay; no P or outcome filter.
                out.append({
                    "cohort": "Fenland_Olink_substitute_485", "assay": aid,
                    "gene": wald["gene"], "uniprot": wald["uniprot"], "locus_type": "cis",
                    "variant": str(raw["rsid"]), "chr": str(raw["chr"]), "pos": raw["pos"],
                    "EA": raw["Allele1"], "OA": raw["Allele2"], "eaf": freq,
                    "beta": beta_raw / sd, "neg_log10_p": logp, "se_reconstructed": "NA",
                    "se_used": se_raw / sd,
                    "se_source": "frozen_source_index_sdY_scaled_raw_StdErr_not_reconstructed",
                    "source_table": raw_path.name, "source_data_row": int(raw["_source_data_row"]),
                    "cohort_wald_exposure_row_id": wald["source_row"],
                    "input_role": "followup_protein_fixed_assay_all_source_cis_rows_no_p_filter",
                    "raw_beta": beta_raw, "raw_se": se_raw, "sdY_estimated": sd,
                    "unit": "per_estimated_SD_not_observed_SD",
                })
        offset += len(chunk)
    counts = pd.Series([x["assay"] for x in out]).value_counts().to_dict()
    absent = selected_ids - set(counts)
    if absent:
        raise ValueError(f"selected Fenland assay has no source rows: {sorted(absent)}")
    for aid in selected_ids:
        index_row = chosen[aid]["index"]
        if str(index_row["rsid"]) not in {x["variant"] for x in out if x["assay"] == aid}:
            raise ValueError(f"Fenland selected index variant absent from full raw assay region: {aid}/{index_row['rsid']}")
    return pd.DataFrame(out)


def prepare():
    helpers = load_old_helpers()
    followups, selected = selected_inputs()
    wald_rows, selected_pairs, pair_rows = {}, {}, {}
    for _, row in selected.iterrows():
        mapped = adapt_selected_wald(row)
        key = (mapped["cohort"], mapped["uniprot"])
        wald_rows[key] = mapped
        selected_pairs[key] = mapped["assay"]
        pair_rows[key] = mapped

    st_tables, st_source_counts = [], {}
    for source_spec in helpers.SOURCE_SPECS:
        frame, path, total = helpers._source_frame(ROOT, source_spec, followups, selected_pairs, pair_rows)
        st_tables.append(frame)
        st_source_counts[source_spec["cohort"]] = {"source_rows": int(total), "selected_loci": len(frame)}

    indices = pd.read_csv(INDICES, sep="\t", dtype={"id": str, "Uniprot_ID": str}, low_memory=False)
    fen_rows = selected[selected["cohort"] == "Fenland_Olink_substitute_485"]
    fen = fenland_records(fen_rows, indices, FENLAND_RAW)
    if fen.empty:
        raise ValueError("No selected Fenland substitute assay rows")
    all_loci = pd.concat(st_tables + [fen], ignore_index=True, sort=False)
    if all_loci.empty or not set(all_loci["locus_type"].dropna().astype(str)) <= {"cis", "trans"}:
        raise ValueError("all_loci empty or contains unexpected locus_type")
    all_loci["se_used"] = all_loci["se_used"].where(all_loci["se_used"].notna(), all_loci["se_reconstructed"])

    OUT.mkdir(parents=True, exist_ok=True)
    loci_path, mvp_path = OUT / "all_loci.tsv", OUT / "mvp_subset.tsv"
    all_loci.to_csv(loci_path, sep="\t", index=False, na_rep="NA")
    requested = helpers.explicit_rsids(all_loci["variant"])
    mvp = helpers.stream_mvp_subset(MVP, mvp_path, requested)
    mvp_summary = {
        "sha256": mvp["sha256"], "data_rows_scanned": mvp["data_rows_scanned"],
        "matched_rows": mvp["matched_rows"], "matched_unique_rsids": len(mvp["matched_unique_rsids"]),
        "missing_rsids": len(mvp["missing_rsids"]), "missing_rsid_examples": mvp["missing_rsids"][:20],
    }
    inputs = [FOLLOWUP, WALD, INDICES, FENLAND_RAW, MVP,
              ROOT / "research/data/processed/exposures/ST19_somascan_pQTLs_smpnorm.tsv",
              ROOT / "research/data/processed/exposures/ST21_olink_ukb_bi_pQTLs.tsv",
              SE_SOURCES]
    row_counts = {str(k): int(v) for k, v in all_loci["cohort"].value_counts().items()}
    type_counts = {str(k): int(v) for k, v in all_loci.groupby(["cohort", "locus_type"], dropna=False).size().items()}
    record = {
        "status": "prepared_sensitivity_substitute_inputs",
        "scope": "17 follow-up UniProt IDs; fixed mvp_cohort_wald assays; all published ST19/ST21 cis/trans rows and all source Fenland selected-assay cis rows; no outcome/P filter",
        "followup_uniprots": len(followups), "selected_wald_rows": len(selected),
        "selected_assays_by_cohort": {str(k): int(v) for k, v in selected.groupby("cohort").size().items()},
        "source_rows_by_cohort": {**st_source_counts, "Fenland_Olink_substitute_485": {"selected_assays": len(fen_rows), "selected_cis_loci": len(fen)}},
        "all_loci_rows": len(all_loci), "all_loci_rows_by_cohort": row_counts,
        "all_loci_rows_by_cohort_and_type": type_counts,
        "fenland_scope": "selected source assays only; cis region pairs from raw summary; no trans source available and none represented",
        "fenland_index_checks": {"selected_rows": len(fen_rows), "source_indices_beta_se_allele_variant_matches": True,
                                 "beta_scale": "Effect/sdY_estimated", "se_scale": "StdErr/sdY_estimated; actual source SE, not reverse-inferred"},
        "mvp_subset": mvp_summary,
        "input_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
        "script_sha256": sha256(Path(__file__)),
        "output_sha256": {loci_path.name: sha256(loci_path), mvp_path.name: sha256(mvp_path)},
        "outputs": [str(loci_path.relative_to(ROOT)), str(mvp_path.relative_to(ROOT))],
    }
    (OUT / "run.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(prepare(), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
