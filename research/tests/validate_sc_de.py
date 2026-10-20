#!/usr/bin/env python3
"""Mechanical consistency checks for existing sc_de and pseudobulk outputs."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXPECTED_GENES = 27678
EXPECTED_CELLS = 46225
RTOL = 1e-8
ATOL = 1e-10


class Audit:
    def __init__(self):
        self.checks = []

    def add(self, name, passed, detail):
        self.checks.append({"name": name, "passed": bool(passed), "detail": str(detail)})

    @property
    def failures(self):
        return [check for check in self.checks if not check["passed"]]


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_manifest_path(root, raw_path, output_base):
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return (output_base if output_base is not None else root) / path


def verify_manifest_hashes(root, run_path, run_data, audit, hash_cache):
    reports = []
    for section in ("inputs", "outputs"):
        entries = run_data.get(section, {})
        if not isinstance(entries, dict):
            audit.add(f"{run_path.name}.{section}_manifest_shape", False, "manifest section is not a path-to-SHA256 object")
            continue
        for raw_path, expected in entries.items():
            path = resolve_manifest_path(root, raw_path, run_path.parent if section == "outputs" else None)
            if not path.is_file():
                actual = None
                passed = False
                detail = "missing file"
                size = None
            else:
                resolved = path.resolve()
                key = str(resolved)
                if key not in hash_cache:
                    hash_cache[key] = sha256_file(resolved)
                actual = hash_cache[key]
                size = resolved.stat().st_size
                passed = actual == expected
                detail = "match" if passed else f"expected {expected}; found {actual}"
            audit.add(f"{run_path.name}.{section}:{raw_path}", passed, detail)
            reports.append({
                "section": section,
                "manifest_path": raw_path,
                "resolved_path": str(path),
                "expected_sha256": expected,
                "actual_sha256": actual,
                "size_bytes": size,
                "passed": passed,
            })
    return reports


def read_tsv(path, **kwargs):
    return pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=True, **kwargs)


def parse_numeric(frame, column, audit, context):
    converted = pd.to_numeric(frame[column], errors="coerce")
    bad_text = frame[column].notna() & converted.isna()
    audit.add(f"{context}.{column}_numeric_parse", not bool(bad_text.any()),
              f"unparseable values={int(bad_text.sum())}")
    return converted.to_numpy(dtype=float)


def bh_adjust(p_values, family_size=None):
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(p_values.shape, np.nan, dtype=float)
    indices = np.flatnonzero(np.isfinite(p_values))
    if not len(indices):
        return adjusted
    order = indices[np.argsort(p_values[indices], kind="mergesort")]
    m = int(family_size if family_size is not None else len(indices))
    ranks = np.arange(1, len(order) + 1, dtype=float)
    raw_adjusted = p_values[order] * m / ranks
    monotone = np.minimum.accumulate(raw_adjusted[::-1])[::-1]
    adjusted[order] = np.minimum(monotone, 1.0)
    return adjusted


def validate_p_domain(values, audit, name):
    finite = np.isfinite(values)
    passed = bool(np.all((values[finite] >= 0) & (values[finite] <= 1)))
    audit.add(name, passed, f"finite={int(finite.sum())}; missing={int((~finite).sum())}")


def validate_de_table(path, comparison, annotated, expected_cell_types, audit):
    frame = read_tsv(path)
    required = {
        "gene_id", "cell_type", "comparison", "p_value", "log2FC", "fdr_bh",
        "n_group1", "n_group0", "mean_normalized_group1", "mean_normalized_group0",
        "pct_group1", "pct_group0", "test_status", "BH_family_size",
    }
    missing_columns = sorted(required - set(frame.columns))
    audit.add(f"{path.name}.required_columns", not missing_columns, f"missing={missing_columns}")
    if missing_columns:
        return {"rows": len(frame), "cell_types": {}, "groups": 0}

    audit.add(f"{path.name}.comparison_value", bool((frame["comparison"] == comparison).all()), comparison)
    p_all = parse_numeric(frame, "p_value", audit, path.name)
    fdr_all = parse_numeric(frame, "fdr_bh", audit, path.name)
    logfc_all = parse_numeric(frame, "log2FC", audit, path.name)
    validate_p_domain(p_all, audit, f"{path.name}.p_value_domain")
    validate_p_domain(fdr_all, audit, f"{path.name}.fdr_bh_domain")
    same_missing = np.array_equal(np.isnan(p_all), np.isnan(fdr_all))
    audit.add(f"{path.name}.p_fdr_missing_mask", same_missing,
              f"p_missing={int(np.isnan(p_all).sum())}; fdr_missing={int(np.isnan(fdr_all).sum())}")

    cell_types = list(pd.unique(frame["cell_type"]))
    audit.add(f"{path.name}.cell_type_set", set(cell_types) == set(expected_cell_types),
              f"observed={cell_types}")
    audit.add(f"{path.name}.row_count", len(frame) == len(expected_cell_types) * EXPECTED_GENES,
              f"observed={len(frame)}; expected={len(expected_cell_types) * EXPECTED_GENES}")

    group_details = {}
    max_bh_error = 0.0
    max_fc_error = 0.0
    bad_bh = 0
    bad_fc = 0
    gene_sets = {}
    all_zero_count = 0
    all_zero_bad = 0
    for cell_type, group in frame.groupby("cell_type", sort=False, dropna=False):
        key = str(cell_type)
        unique_genes = group["gene_id"].nunique(dropna=False)
        exact_genes = len(group) == EXPECTED_GENES and unique_genes == EXPECTED_GENES
        audit.add(f"{path.name}.{key}.complete_unique_genes", exact_genes,
                  f"rows={len(group)}; unique_gene_id={unique_genes}; expected={EXPECTED_GENES}")
        gene_sets[key] = set(group["gene_id"].astype(str))

        raw_p = pd.to_numeric(group["p_value"], errors="coerce").to_numpy(dtype=float)
        observed_fdr = pd.to_numeric(group["fdr_bh"], errors="coerce").to_numpy(dtype=float)
        if comparison == "type_vs_rest":
            family_size = EXPECTED_GENES
        else:
            family_size = EXPECTED_GENES
        expected_fdr = bh_adjust(raw_p, family_size=family_size)
        valid = np.isfinite(expected_fdr) & np.isfinite(observed_fdr)
        equal = np.isclose(expected_fdr[valid], observed_fdr[valid], rtol=RTOL, atol=ATOL)
        local_bad = int((~equal).sum()) + int((np.isfinite(expected_fdr) != np.isfinite(observed_fdr)).sum())
        bad_bh += local_bad
        if valid.any():
            max_bh_error = max(max_bh_error, float(np.max(np.abs(expected_fdr[valid] - observed_fdr[valid]))))

        n1 = pd.to_numeric(group["n_group1"], errors="coerce").to_numpy(dtype=float)
        n0 = pd.to_numeric(group["n_group0"], errors="coerce").to_numpy(dtype=float)
        mean1 = pd.to_numeric(group["mean_normalized_group1"], errors="coerce").to_numpy(dtype=float)
        mean0 = pd.to_numeric(group["mean_normalized_group0"], errors="coerce").to_numpy(dtype=float)
        pct1 = pd.to_numeric(group["pct_group1"], errors="coerce").to_numpy(dtype=float)
        pct0 = pd.to_numeric(group["pct_group0"], errors="coerce").to_numpy(dtype=float)
        if np.any(n1 <= 0) or np.any(n0 <= 0) or np.any(mean1 < 0) or np.any(mean0 < 0):
            bad_fc += len(group)
            expected_fc = np.full(len(group), np.nan)
        else:
            expected_fc = np.log2(mean1 + 1.0 / n1) - np.log2(mean0 + 1.0 / n0)
            actual_fc = pd.to_numeric(group["log2FC"], errors="coerce").to_numpy(dtype=float)
            fc_equal = np.isclose(expected_fc, actual_fc, rtol=RTOL, atol=ATOL)
            bad_fc += int((~fc_equal).sum())
            max_fc_error = max(max_fc_error, float(np.max(np.abs(expected_fc - actual_fc))))
        all_zero = (mean1 == 0) & (mean0 == 0) & (pct1 == 0) & (pct0 == 0)
        all_zero_count += int(all_zero.sum())
        all_zero_bad += int((all_zero & (np.isfinite(raw_p) | np.isfinite(observed_fdr))).sum())

        n1_values = pd.unique(n1)
        n0_values = pd.unique(n0)
        if comparison == "type_vs_rest":
            expected_n1 = int((annotated["cell_type"] == cell_type).sum())
            expected_n0 = EXPECTED_CELLS - expected_n1
            sample_size_ok = len(n1_values) == 1 and len(n0_values) == 1 and int(n1_values[0]) == expected_n1 and int(n0_values[0]) == expected_n0 and expected_n1 + expected_n0 == EXPECTED_CELLS
            sample_detail = f"n1={n1_values.tolist()}; n0={n0_values.tolist()}; annotated type n={expected_n1}; total={expected_n1 + expected_n0}"
        else:
            sub = annotated[annotated["cell_type"] == cell_type]
            expected_n1 = int((sub["condition"] == "psoriasis").sum())
            expected_n0 = int((sub["condition"] == "normal").sum())
            sample_size_ok = len(n1_values) == 1 and len(n0_values) == 1 and int(n1_values[0]) == expected_n1 and int(n0_values[0]) == expected_n0
            sample_detail = f"n1={n1_values.tolist()}; n0={n0_values.tolist()}; annotated psoriasis={expected_n1}; normal={expected_n0}"
        audit.add(f"{path.name}.{key}.sample_sizes", sample_size_ok, sample_detail)
        group_details[key] = {
            "rows": len(group),
            "unique_gene_ids": int(unique_genes),
            "n_group1": n1_values.tolist(),
            "n_group0": n0_values.tolist(),
        }

    reference = next(iter(gene_sets.values()), set())
    same_gene_set = all(values == reference for values in gene_sets.values())
    audit.add(f"{path.name}.gene_set_same_across_cell_types", same_gene_set,
              f"cell_type_count={len(gene_sets)}")
    audit.add(f"{path.name}.bh_recalculation", bad_bh == 0,
              f"N={EXPECTED_GENES}; mismatches={bad_bh}; max_abs_error={max_bh_error:.3g}; rtol={RTOL}; atol={ATOL}")
    audit.add(f"{path.name}.log2fc_recalculation", bad_fc == 0,
              f"mismatches={bad_fc}; max_abs_error={max_fc_error:.3g}; formula=log2(mean_group1+1/n1)-log2(mean_group0+1/n0)")
    audit.add(f"{path.name}.all_zero_p_fdr_missing", all_zero_bad == 0,
              f"all_zero_rows={all_zero_count}; violations={all_zero_bad}")
    return {"rows": len(frame), "cell_types": group_details, "groups": len(gene_sets), "all_zero_rows": all_zero_count}


def validate_author_s5(path, audit):
    frame = read_tsv(path, dtype=str)
    required = {
        "author_Gene", "author_cell type", "comparison_status", "actual_gene_id",
        "actual_p_value", "actual_log2FC", "actual_fdr_bh",
        "author_passes_reported_threshold", "actual_passes_reported_threshold",
        "direction_agrees_with_author", "threshold_agrees_with_author",
    }
    missing = sorted(required - set(frame.columns))
    audit.add("author_S5_comparison.required_columns", not missing, f"missing={missing}")
    audit.add("author_S5_comparison.row_count", len(frame) == 26, f"observed={len(frame)}; expected=26")
    status_counts = frame["comparison_status"].value_counts(dropna=False).to_dict()
    non_unique = frame["comparison_status"] != "unique_match"
    fields_to_be_missing = [
        "actual_gene_id", "actual_p_value", "actual_log2FC", "actual_fdr_bh",
        "author_passes_reported_threshold", "actual_passes_reported_threshold",
        "direction_agrees_with_author", "threshold_agrees_with_author",
    ]
    missing_semantics = True
    violations = []
    for column in fields_to_be_missing:
        invalid = non_unique & frame[column].notna() & (frame[column] != "") & (frame[column] != "NA")
        if invalid.any():
            missing_semantics = False
            violations.append(f"{column}:{int(invalid.sum())}")
    flags = ["author_passes_reported_threshold", "actual_passes_reported_threshold", "direction_agrees_with_author", "threshold_agrees_with_author"]
    false_missing = sum(int((non_unique & (frame[column].fillna("").astype(str).str.upper() == "FALSE")).sum()) for column in flags)
    audit.add("author_S5_comparison.nonunique_missing_is_not_false", missing_semantics and false_missing == 0,
              f"status_counts={status_counts}; invalid_fields={violations}; false_flags_on_nonunique={false_missing}")
    unique_rows = frame[~non_unique]
    unique_has_actual = unique_rows["actual_gene_id"].notna() & (unique_rows["actual_gene_id"] != "")
    audit.add("author_S5_comparison.unique_has_actual_gene", bool(unique_has_actual.all()),
              f"unique_rows={len(unique_rows)}; missing_actual={int((~unique_has_actual).sum())}")
    return {"rows": len(frame), "status_counts": {str(key): int(value) for key, value in status_counts.items()}}


def validate_pseudobulk(root, annotated, de_run, pb_run, audit):
    counts_path = root / "research/results/sc_de/pseudobulk_raw_counts.tsv.gz"
    samples_path = root / "research/results/sc_de/pseudobulk_samples.tsv"
    inclusion_path = root / "research/results/sc_pseudobulk/sample_inclusion.tsv"
    counts = read_tsv(counts_path, dtype={"gene_id": str})
    samples = read_tsv(samples_path, dtype={"sample_id": str, "donor_id": str, "cell_type": str, "condition": str})
    inclusion = read_tsv(inclusion_path, dtype={"sample_id": str, "donor_id": str, "cell_type": str, "condition": str})
    gene_ids = counts["gene_id"].astype(str)
    matrix = counts.drop(columns=["gene_id"]).apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    sample_ids = list(counts.columns[1:])
    all_integer = bool(np.isfinite(matrix).all() and np.all(matrix >= 0) and np.all(matrix == np.floor(matrix)))
    audit.add("pseudobulk.counts_nonnegative_integer", all_integer,
              f"shape={list(matrix.shape)}; finite={bool(np.isfinite(matrix).all())}; min={float(np.nanmin(matrix))}; max={float(np.nanmax(matrix))}")
    audit.add("pseudobulk.gene_id_unique", gene_ids.nunique() == len(gene_ids),
              f"rows={len(gene_ids)}; unique={gene_ids.nunique()}")
    audit.add("pseudobulk.sample_id_unique", samples["sample_id"].nunique() == len(samples),
              f"rows={len(samples)}; unique={samples['sample_id'].nunique()}")
    audit.add("pseudobulk.matrix_columns_equal_samples_order", sample_ids == samples["sample_id"].tolist(),
              f"matrix_columns={len(sample_ids)}; sample_rows={len(samples)}")
    audit.add("pseudobulk.sample_column_is_donor_type", not samples.duplicated(["donor_id", "cell_type"]).any(),
              f"duplicate donor/type rows={int(samples.duplicated(['donor_id', 'cell_type']).sum())}")

    group_columns = ["donor_id", "cell_type", "condition"]
    metadata_groups = annotated.groupby(group_columns, dropna=False).agg(
        metadata_cells=("cell_id", "size"), metadata_counts=("n_counts", "sum")
    )
    cell_size_errors = []
    count_errors = []
    sample_keys = []
    for index, row in samples.iterrows():
        key = (str(row["donor_id"]), str(row["cell_type"]), str(row["condition"]))
        sample_keys.append(key)
        if key not in metadata_groups.index:
            cell_size_errors.append(f"{row['sample_id']}: absent metadata group {key}")
            continue
        expected_cells = int(metadata_groups.loc[key, "metadata_cells"])
        expected_counts = int(metadata_groups.loc[key, "metadata_counts"])
        if int(row["n_cells"]) != expected_cells:
            cell_size_errors.append(f"{row['sample_id']}: table={row['n_cells']} metadata={expected_cells}")
        observed_counts = int(matrix[:, index].sum(dtype=np.int64))
        if observed_counts != expected_counts:
            count_errors.append(f"{row['sample_id']}: matrix={observed_counts} metadata={expected_counts}")
    audit.add("pseudobulk.n_cells_match_annotated_metadata", not cell_size_errors,
              f"samples={len(samples)}; mismatches={cell_size_errors[:5]}")
    audit.add("pseudobulk.column_counts_match_annotated_metadata", not count_errors,
              f"samples={len(samples)}; mismatches={count_errors[:5]}")
    annotated_total = int(pd.to_numeric(annotated["n_counts"], errors="raise").sum())
    matrix_total = int(matrix.sum(dtype=np.int64))
    audit.add("pseudobulk.total_counts_include_all_annotated_cells", matrix_total == annotated_total,
              f"matrix_sum={matrix_total}; annotated_n_counts_sum={annotated_total}; Unassigned_cells={int((annotated['cell_type'] == 'Unassigned').sum())}")
    audit.add("pseudobulk.includes_unassigned_group", bool((samples["cell_type"] == "Unassigned").any()),
              f"Unassigned_sample_groups={int((samples['cell_type'] == 'Unassigned').sum())}")

    inclusion_columns = ["sample_id", "donor_id", "cell_type", "condition", "n_cells", "included"]
    missing_inclusion_columns = sorted(set(inclusion_columns) - set(inclusion.columns))
    audit.add("pseudobulk.sample_inclusion_columns", not missing_inclusion_columns, f"missing={missing_inclusion_columns}")
    if not missing_inclusion_columns:
        joined = samples[inclusion_columns[:-1]].merge(
            inclusion[inclusion_columns], on="sample_id", how="outer", suffixes=("_sample", "_inclusion"), indicator=True
        )
        merge_ok = bool((joined["_merge"] == "both").all())
        fields_ok = True
        for column in ("donor_id", "cell_type", "condition", "n_cells"):
            fields_ok &= bool((joined[f"{column}_sample"] == joined[f"{column}_inclusion"]).all())
        audit.add("pseudobulk.inclusion_matches_sample_manifest", merge_ok and fields_ok,
                  f"sample_inclusion_rows={len(inclusion)}; mismatched_rows={int((joined['_merge'] != 'both').sum())}")

    donor_path = root / "research/results/sc_pseudobulk/donor_level_differential_expression.tsv"
    donor = read_tsv(donor_path)
    donor_run_result = pb_run.get("result", {})
    tested_reported = int(donor_run_result.get("tested_gene_contrasts", -1))
    significant_reported = int(donor_run_result.get("significant_gene_contrasts", -1))
    donor_sig = donor["significant"].astype(str).str.upper() == "TRUE"
    audit.add("pseudobulk.donor_result_row_count", len(donor) == tested_reported,
              f"table={len(donor)}; run_manifest={tested_reported}")
    audit.add("pseudobulk.donor_result_significant_count", int(donor_sig.sum()) == significant_reported,
              f"table={int(donor_sig.sum())}; run_manifest={significant_reported}")
    donor_p = parse_numeric(donor, "PValue", audit, donor_path.name)
    donor_fdr = parse_numeric(donor, "FDR", audit, donor_path.name)
    validate_p_domain(donor_p, audit, f"{donor_path.name}.PValue_domain")
    validate_p_domain(donor_fdr, audit, f"{donor_path.name}.FDR_domain")
    donor_bh_bad = 0
    donor_bh_max_error = 0.0
    donor_type_counts = {}
    for (cell_type, comparison), group in donor.groupby(["cell_type", "comparison"], sort=False):
        p = pd.to_numeric(group["PValue"], errors="coerce").to_numpy(dtype=float)
        q = pd.to_numeric(group["FDR"], errors="coerce").to_numpy(dtype=float)
        family_size = int(np.isfinite(p).sum())
        expected = bh_adjust(p, family_size=family_size)
        observed = q
        valid = np.isfinite(expected) & np.isfinite(observed)
        donor_bh_bad += int((~np.isclose(expected[valid], observed[valid], rtol=RTOL, atol=ATOL)).sum())
        donor_bh_bad += int((np.isfinite(expected) != np.isfinite(observed)).sum())
        if valid.any():
            donor_bh_max_error = max(donor_bh_max_error, float(np.max(np.abs(expected[valid] - observed[valid]))))
        donor_type_counts[str(cell_type)] = {
            "tested_rows": len(group), "tested_p_values": family_size,
            "significant_rows": int((group["significant"].astype(str).str.upper() == "TRUE").sum()),
        }
    audit.add("pseudobulk.donor_bh_recalculation", donor_bh_bad == 0,
              f"actual tested families; mismatches={donor_bh_bad}; max_abs_error={donor_bh_max_error:.3g}")

    comparison_path = root / "research/results/sc_pseudobulk/cell_vs_donor_comparison.tsv"
    comparison = read_tsv(comparison_path, dtype={"donor_result_missing_reason": str})
    missing_reason = comparison["donor_result_missing_reason"].fillna("").astype(str).str.strip() != ""
    donor_missing_fields = ["logFC", "PValue", "FDR", "significant_donor", "direction_agrees"]
    missing_donor_false = 0
    missing_donor_non_na = 0
    for column in donor_missing_fields:
        values = comparison.loc[missing_reason, column]
        missing_donor_non_na += int(values.notna().sum())
        if column in ("significant_donor", "direction_agrees"):
            missing_donor_false += int(values.astype(str).str.upper().eq("FALSE").sum())
    audit.add("pseudobulk.missing_donor_not_filled_false", missing_donor_false == 0 and missing_donor_non_na == 0,
              f"missing_donor_rows={int(missing_reason.sum())}; non_NA_values={missing_donor_non_na}; False_flags={missing_donor_false}")
    present_reason = ~missing_reason
    donor_missing_without_reason = comparison.loc[present_reason, ["PValue", "FDR"]].isna().any(axis=1).sum()
    audit.add("pseudobulk.donor_missing_reason_consistency", int(donor_missing_without_reason) == 0,
              f"nonmissing_reason_rows={int(present_reason.sum())}; rows_missing_P_or_FDR={int(donor_missing_without_reason)}")

    return {
        "matrix_shape": list(matrix.shape),
        "sample_columns": len(sample_ids),
        "sample_inclusion_rows": len(inclusion),
        "n_counts_total": matrix_total,
        "annotated_n_counts_total": annotated_total,
        "donor_result_rows": len(donor),
        "donor_tested_by_cell_type": donor_type_counts,
        "missing_donor_rows": int(missing_reason.sum()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    root = args.root.resolve()
    audit = Audit()
    hash_cache = {}
    sc_de_dir = root / "research/results/sc_de"
    pb_dir = root / "research/results/sc_pseudobulk"
    sc_de_run_path = sc_de_dir / "run.json"
    pb_run_path = pb_dir / "run.json"
    sc_de_run = json.loads(sc_de_run_path.read_text(encoding="utf-8"))
    pb_run = json.loads(pb_run_path.read_text(encoding="utf-8"))
    sc_hashes = verify_manifest_hashes(root, sc_de_run_path, sc_de_run, audit, hash_cache)
    pb_hashes = verify_manifest_hashes(root, pb_run_path, pb_run, audit, hash_cache)

    de_result = sc_de_run.get("result", {})
    gene_count = int(de_result.get("genes", -1))
    cell_count = int(de_result.get("cells", -1))
    audit.add("run_manifest.fixed_gene_count", gene_count == EXPECTED_GENES, f"run.json genes={gene_count}")
    audit.add("run_manifest.fixed_cell_count", cell_count == EXPECTED_CELLS, f"run.json cells={cell_count}")
    audit.add("run_manifest_rds_not_deserialized", True,
              "RDS input SHA256 is streamed from disk; no RDS deserialization occurs")

    annotated_path = sc_de_dir / "annotated_cells.tsv"
    annotated = read_tsv(annotated_path, dtype={"gsm": str, "barcode": str, "cell_id": str, "donor_id": str, "cell_type": str, "condition": str})
    annotated_cells = len(annotated)
    annotated_types = annotated["cell_type"].value_counts(dropna=False).to_dict()
    annotated_unassigned = int((annotated["cell_type"] == "Unassigned").sum())
    audit.add("annotated_cells.unique_cell_id", annotated["cell_id"].nunique() == annotated_cells,
              f"rows={annotated_cells}; unique_cell_id={annotated['cell_id'].nunique()}")
    audit.add("annotated_cells.run_cell_count", annotated_cells == cell_count,
              f"table={annotated_cells}; run={cell_count}")

    type_path = sc_de_dir / "all_gene_type_vs_rest.tsv"
    disease_path = sc_de_dir / "all_gene_psoriasis_vs_normal.tsv"
    expected_cell_types = list(de_result.get("tested_cell_types", []))
    type_summary = validate_de_table(type_path, "type_vs_rest", annotated, expected_cell_types, audit)
    disease_summary = validate_de_table(disease_path, "psoriasis_vs_normal", annotated, expected_cell_types, audit)
    sample_sizes_path = sc_de_dir / "contrast_sample_sizes.tsv"
    sample_sizes = read_tsv(sample_sizes_path)
    sample_sizes_ok = len(sample_sizes) == len(expected_cell_types)
    if sample_sizes_ok:
        for row in sample_sizes.itertuples(index=False):
            group = annotated[annotated["cell_type"] == row.cell_type]
            sample_sizes_ok &= (
                len(group) == int(row.cells)
                and int((group["condition"] == "normal").sum()) == int(row.normal)
                and int((group["condition"] == "psoriasis").sum()) == int(row.psoriasis)
            )
    audit.add("sc_de.contrast_sample_sizes_match_annotated_cells", sample_sizes_ok,
              f"types={len(sample_sizes)}; unassigned_cells={annotated_unassigned}; annotated_total={annotated_cells}")

    author_path = sc_de_dir / "author_S5_comparison.tsv"
    author_summary = validate_author_s5(author_path, audit)
    pseudobulk_summary = validate_pseudobulk(root, annotated, sc_de_run, pb_run, audit)

    script_path = Path(__file__).resolve()
    result = {
        "status": "passed_mechanical_checks" if not audit.failures else "checks_failed",
        "scope": "Mechanical file/hash/numeric consistency only; no statistical rerun, RDS loading, or scientific interpretation.",
        "validated_script_sha256": sha256_file(script_path),
        "RDS_policy": "The 3.5GB checkpoint was not deserialized; any manifest hash was computed by streaming chunks.",
        "sc_de_manifest_hash_checks": sc_hashes,
        "pseudobulk_manifest_hash_checks": pb_hashes,
        "sc_de_tables": {"type_vs_rest": type_summary, "psoriasis_vs_normal": disease_summary},
        "annotated_cells": {"rows": annotated_cells, "unassigned": annotated_unassigned,
                            "cell_type_counts": {str(key): int(value) for key, value in annotated_types.items()}},
        "author_S5_comparison": author_summary,
        "pseudobulk": pseudobulk_summary,
        "checks_passed": sum(check["passed"] for check in audit.checks),
        "checks_failed": len(audit.failures),
        "checks": audit.checks,
        "failures": audit.failures,
    }
    output = root / "research/records/sc_de_validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "checks_passed", "checks_failed", "author_S5_comparison", "pseudobulk"
    )}, ensure_ascii=False, indent=2))
    if audit.failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
