#!/usr/bin/env python3
"""Compare author-reported S5 rows with a supplied psoriasis-vs-normal contrast."""

import argparse
import csv
import json
import math
from pathlib import Path


AUTHOR_COLUMNS = (
    "Gene", "cell type", "logFC", "p.value", "Psoriasis _Mean",
    "Normal_Mean", "adjusted_pvalue", "qval",
)
CONTRAST_COLUMNS = ("gene_id", "symbol", "cell_type", "comparison", "p_value", "log2FC", "fdr_bh")
FEATURE_COLUMNS = ("gene_id", "symbol", "mapping_status")
ALIASES = {"DDX58": "RIGI"}
LOG2FC_THRESHOLD = 0.25
FDR_THRESHOLD = 0.05
EVIDENCE_ORIGIN = "author_table_vs_supplied_sc_contrast_table_only"


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty TSV: {path}")
        return reader.fieldnames, list(reader)


def _require_columns(path, fieldnames, required):
    missing = [name for name in required if name not in fieldnames]
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(missing)}")


def _number(value):
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _format_number(value):
    return "" if value is None else format(value, ".15g")


def _formula_diagnostic(author):
    psoriasis = _number(author["Psoriasis _Mean"])
    normal = _number(author["Normal_Mean"])
    author_logfc = _number(author["logFC"])
    if psoriasis is None or normal is None or psoriasis <= 0 or normal <= 0:
        return "invalid_or_nonpositive_mean", "", ""
    formula = math.log2(psoriasis / normal)
    difference = None if author_logfc is None else author_logfc - formula
    return "formula_diagnostic", _format_number(formula), _format_number(difference)


def _feature_labels(row, feature_rows_by_id):
    labels = set()
    feature_rows = feature_rows_by_id.get(row["gene_id"], [])
    if row["symbol"]:
        labels.add(row["symbol"])
    if len(feature_rows) == 1 and feature_rows[0]["mapping_status"] == "unique" and feature_rows[0]["symbol"]:
        labels.add(feature_rows[0]["symbol"])
    return labels


def _find_match(author, contrast_rows, feature_rows_by_id):
    label = author["Gene"]
    author_cell = author["cell type"]
    cell_rows = [row for row in contrast_rows if row["cell_type"] == author_cell]
    if not cell_rows:
        return "no_exact_cell_label", "", []

    exact = [row for row in cell_rows if label in _feature_labels(row, feature_rows_by_id)]
    method = "exact_symbol"
    candidates = exact
    if not candidates and label in ALIASES:
        alias = ALIASES[label]
        candidates = [row for row in cell_rows if alias in _feature_labels(row, feature_rows_by_id)]
        method = "explicit_alias"
    if len(candidates) > 1:
        return "multiple_matches", method, []
    if len(candidates) == 1:
        return "unique_match", method, candidates

    unresolved_feature = [
        row for row in cell_rows
        if label in {row["symbol"], ALIASES.get(label, "")}
        and (row["gene_id"] not in feature_rows_by_id or len(feature_rows_by_id[row["gene_id"]]) != 1)
    ]
    if unresolved_feature:
        return "missing_or_multiple_feature_mapping", "", []
    return "missing_result", "", []


def _compare_row(author, match_status, match_method, actual_rows):
    row = {f"author_{key}": author[key] for key in AUTHOR_COLUMNS}
    formula_status, formula_logfc, formula_difference = _formula_diagnostic(author)
    row.update({
        "comparison_status": match_status,
        "match_method": match_method,
        "evidence_origin": EVIDENCE_ORIGIN,
        "actual_gene_id": "",
        "actual_symbol": "",
        "actual_cell_type": "",
        "actual_comparison": "",
        "actual_p_value": "",
        "actual_log2FC": "",
        "actual_fdr_bh": "",
        "author_passes_reported_threshold": "NA",
        "actual_passes_reported_threshold": "NA",
        "direction_agrees_with_author": "NA",
        "threshold_agrees_with_author": "NA",
        "formula_diagnostic_status": formula_status,
        "formula_log2FC": formula_logfc,
        "author_logFC_minus_formula_log2FC": formula_difference,
    })
    if match_status != "unique_match":
        return row

    actual = actual_rows[0]
    author_logfc = _number(author["logFC"])
    author_fdr = _number(author["adjusted_pvalue"])
    actual_logfc = _number(actual["log2FC"])
    actual_fdr = _number(actual["fdr_bh"])
    author_pass = None if author_logfc is None or author_fdr is None else abs(author_logfc) >= LOG2FC_THRESHOLD and author_fdr < FDR_THRESHOLD
    actual_pass = None if actual_logfc is None or actual_fdr is None else abs(actual_logfc) >= LOG2FC_THRESHOLD and actual_fdr < FDR_THRESHOLD
    if author_logfc is None or actual_logfc is None or author_logfc == 0 or actual_logfc == 0:
        direction_agrees = "indeterminate_zero_or_missing_effect"
    else:
        direction_agrees = str((author_logfc > 0) == (actual_logfc > 0))
    threshold_agrees = "NA" if author_pass is None or actual_pass is None else str(author_pass == actual_pass)
    row.update({
        "actual_gene_id": actual["gene_id"],
        "actual_symbol": actual["symbol"],
        "actual_cell_type": actual["cell_type"],
        "actual_comparison": actual["comparison"],
        "actual_p_value": actual["p_value"],
        "actual_log2FC": actual["log2FC"],
        "actual_fdr_bh": actual["fdr_bh"],
        "author_passes_reported_threshold": "NA" if author_pass is None else str(author_pass),
        "actual_passes_reported_threshold": "NA" if actual_pass is None else str(actual_pass),
        "direction_agrees_with_author": direction_agrees,
        "threshold_agrees_with_author": threshold_agrees,
    })
    return row


def compare_tables(contrast_path, author_path, features_path, output_path):
    contrast_header, contrast = read_tsv(contrast_path)
    author_header, authors = read_tsv(author_path)
    feature_header, features = read_tsv(features_path)
    _require_columns(contrast_path, contrast_header, CONTRAST_COLUMNS)
    _require_columns(author_path, author_header, AUTHOR_COLUMNS)
    _require_columns(features_path, feature_header, FEATURE_COLUMNS)

    feature_rows_by_id = {}
    for row in features:
        feature_rows_by_id.setdefault(row["gene_id"], []).append(row)
    psoriasis_contrast = [row for row in contrast if row["comparison"] == "psoriasis_vs_normal"]
    output_rows = []
    for author in authors:
        status, method, matches = _find_match(author, psoriasis_contrast, feature_rows_by_id)
        output_rows.append(_compare_row(author, status, method, matches))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [f"author_{key}" for key in AUTHOR_COLUMNS] + [
        "comparison_status", "match_method", "evidence_origin", "actual_gene_id", "actual_symbol",
        "actual_cell_type", "actual_comparison", "actual_p_value", "actual_log2FC",
        "actual_fdr_bh", "author_passes_reported_threshold",
        "actual_passes_reported_threshold", "direction_agrees_with_author",
        "threshold_agrees_with_author", "formula_diagnostic_status", "formula_log2FC",
        "author_logFC_minus_formula_log2FC",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    return output_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contrast", type=Path, required=True, help="Full gene-level contrast TSV")
    parser.add_argument("--author", type=Path, required=True, help="Author table_S5.tsv")
    parser.add_argument("--features", type=Path, required=True, help="features.tsv from sc matrix preparation")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = compare_tables(args.contrast, args.author, args.features, args.output)
    counts = {}
    for row in rows:
        status = row["comparison_status"]
        counts[status] = counts.get(status, 0) + 1
    print(json.dumps({
        "output": str(args.output),
        "author_rows": len(rows),
        "comparison_status_counts": counts,
        "evidence_origin": "table_comparison_only; formula_diagnostic is not an expression reanalysis",
        "threshold": {"absolute_log2FC_minimum": LOG2FC_THRESHOLD, "fdr_bh_exclusive_maximum": FDR_THRESHOLD},
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
