#!/usr/bin/env python3
"""Audit whether rounded MVP OR/CI/P fields can agree under a Wald CI model."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from decimal import Decimal, InvalidOperation
from pathlib import Path
from statistics import NormalDist


SOURCE_REL = "data/processed/outcomes/mvp_original_rows.tsv"
CURRENT_REL = "results/mr_fenland_substitute/mvp_cohort_wald.tsv"
OUTPUT_REL = "results/mvp_precision_audit"
REQUIRED_SOURCE = ("rsid", "odds_ratio", "ci_lower", "ci_upper", "p_value")
Q975 = NormalDist().inv_cdf(0.975)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def raw_decimal(text: str) -> Decimal | None:
    if text.strip().upper() in {"", "NA", "NAN", "N/A"}:
        return None
    try:
        value = Decimal(text)
    except InvalidOperation:
        return None
    return value if value.is_finite() else None


def half_ulp(text: str) -> Decimal | None:
    value = raw_decimal(text)
    if value is None:
        return None
    return Decimal("0.5") * (Decimal(10) ** value.as_tuple().exponent)


def fmt(value: float | None) -> str:
    return "" if value is None else format(value, ".17g")


def decimal_text(value: Decimal | None) -> str:
    return "" if value is None else str(value)


def audit_row(row: dict[str, str], record_number: int, current_ids: set[str]) -> dict[str, str]:
    original = {name: row.get(name, "") for name in REQUIRED_SOURCE}
    out: dict[str, str] = {
        "source_record": str(record_number),
        **{f"original_{name}": value for name, value in original.items()},
        "currently_used_rsid": "true" if original["rsid"] in current_ids else "false",
        "status": "unresolved",
        "reason": "",
        "beta_abs_min": "",
        "beta_abs_max": "",
        "se_min": "",
        "se_max": "",
        "p_ci_point_estimate": "",
        "p_ci_feasible_lower": "",
        "p_ci_feasible_upper": "",
        "or_half_ulp": decimal_text(half_ulp(original["odds_ratio"])),
        "ci_lower_half_ulp": decimal_text(half_ulp(original["ci_lower"])),
        "ci_upper_half_ulp": decimal_text(half_ulp(original["ci_upper"])),
        "p_half_ulp": decimal_text(half_ulp(original["p_value"])),
        "p_rounding_lower": "",
        "p_rounding_upper": "",
        "p_rounding_intersects_ci_feasible": "unresolved",
    }

    values = {name: raw_decimal(original[name]) for name in ("odds_ratio", "ci_lower", "ci_upper", "p_value")}
    if any(values[name] is None for name in ("odds_ratio", "ci_lower", "ci_upper")):
        out["reason"] = "missing_or_nonfinite_OR_or_CI"
        return out
    odds_ratio, lower, upper = (values[name] for name in ("odds_ratio", "ci_lower", "ci_upper"))
    if odds_ratio <= 0 or lower <= 0 or upper <= lower:
        out["reason"] = "OR_or_CI_not_positive_ordered"
        return out

    or_half, lower_half, upper_half = (half_ulp(original[name]) for name in ("odds_ratio", "ci_lower", "ci_upper"))
    assert or_half is not None and lower_half is not None and upper_half is not None
    or_min, or_max = odds_ratio - or_half, odds_ratio + or_half
    lower_min, lower_max = lower - lower_half, lower + lower_half
    upper_min, upper_max = upper - upper_half, upper + upper_half
    if min(or_min, lower_min, upper_min) <= 0:
        out["reason"] = "rounding_interval_reaches_nonpositive_value"
        return out
    if upper_min <= lower_max or upper_max <= lower_min:
        out["reason"] = "rounding_intervals_do_not_define_positive_width_CI"
        return out

    beta_abs_min = 0.0 if or_min <= 1 <= or_max else min(abs(math.log(float(or_min))), abs(math.log(float(or_max))))
    beta_abs_max = max(abs(math.log(float(or_min))), abs(math.log(float(or_max))))
    se_min = (math.log(float(upper_min)) - math.log(float(lower_max))) / (2 * Q975)
    se_max = (math.log(float(upper_max)) - math.log(float(lower_min))) / (2 * Q975)
    if not all(math.isfinite(x) for x in (beta_abs_min, beta_abs_max, se_min, se_max)) or se_min <= 0 or se_max < se_min:
        out["reason"] = "invalid_SE_or_beta_bounds"
        return out

    p_ci_lower = math.erfc(beta_abs_max / se_min / math.sqrt(2))
    p_ci_upper = math.erfc(beta_abs_min / se_max / math.sqrt(2))
    beta_point = abs(math.log(float(odds_ratio)))
    se_point = (math.log(float(upper)) - math.log(float(lower))) / (2 * Q975)
    p_ci_point = math.erfc(beta_point / se_point / math.sqrt(2))
    out.update({
        "beta_abs_min": fmt(beta_abs_min),
        "beta_abs_max": fmt(beta_abs_max),
        "se_min": fmt(se_min),
        "se_max": fmt(se_max),
        "p_ci_point_estimate": fmt(p_ci_point),
        "p_ci_feasible_lower": fmt(p_ci_lower),
        "p_ci_feasible_upper": fmt(p_ci_upper),
    })
    p = values["p_value"]
    p_u = half_ulp(original["p_value"])
    if p is None or p_u is None or p < 0 or p > 1:
        out["reason"] = "missing_or_out_of_range_P"
        return out
    p_low = max(Decimal(0), p - p_u)
    p_high = min(Decimal(1), p + p_u)
    out["p_rounding_lower"] = str(p_low)
    out["p_rounding_upper"] = str(p_high)
    intersects = float(p_low) <= p_ci_upper and float(p_high) >= p_ci_lower
    out["p_rounding_intersects_ci_feasible"] = "true" if intersects else "false"
    out["status"] = "consistent" if intersects else "inconsistent"
    out["reason"] = "rounding_intervals_overlap" if intersects else "no_overlap_under_Wald_95CI"
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2], help="project root (defaults to this checkout)")
    args = parser.parse_args()
    project_root = args.root.resolve()
    research = project_root / "research"
    source = research / SOURCE_REL
    current = research / CURRENT_REL
    if not source.is_file() or not current.is_file():
        missing = [str(path) for path in (source, current) if not path.is_file()]
        raise FileNotFoundError("Required audit input(s) missing:\n" + "\n".join(missing))

    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        absent = [field for field in REQUIRED_SOURCE if field not in (reader.fieldnames or [])]
        if absent:
            raise ValueError(f"Original MVP rows missing required columns: {absent}")
        source_rows = list(reader)
    with current.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "rsid" not in (reader.fieldnames or []):
            raise ValueError("Current MVP cohort table lacks required rsid column")
        current_ids = {row["rsid"] for row in reader if row.get("rsid", "")}

    source_records = len(source_rows)
    covered = {r['rsid'] for r in source_rows} & current_ids
    missing = current_ids-covered
    full_source = research/'data/raw/gwas/GCST90476190.tsv.gz'
    extra=[]
    if missing:
        with gzip.open(full_source,'rt') as handle:
            header=handle.readline().rstrip('\r\n').split('\t');rs_col=header.index('rsid')
            for number,line in enumerate(handle,start=1):
                fields=line.rstrip('\r\n').split('\t')
                if fields[rs_col] in missing:extra.append((number,dict(zip(header,fields))))
    records = [dict(audit_row(row,i,current_ids),source_file=SOURCE_REL) for i,row in enumerate(source_rows,start=1)]
    records.extend(dict(audit_row(row,i,current_ids),source_file='data/raw/gwas/GCST90476190.tsv.gz') for i,row in extra)
    destination = research / OUTPUT_REL
    destination.mkdir(parents=True, exist_ok=True)
    audit_path = destination / "mvp_precision_audit.tsv"
    anomaly_path = destination / "anomalies.tsv"
    fields = list(records[0]) if records else ["source_record", *[f"original_{name}" for name in REQUIRED_SOURCE], "currently_used_rsid", "status", "reason", "beta_abs_min", "beta_abs_max", "se_min", "se_max", "p_ci_point_estimate", "p_ci_feasible_lower", "p_ci_feasible_upper", "or_half_ulp", "ci_lower_half_ulp", "ci_upper_half_ulp", "p_half_ulp", "p_rounding_lower", "p_rounding_upper", "p_rounding_intersects_ci_feasible"]
    for path, selected in ((audit_path, records), (anomaly_path, [r for r in records if r["status"] != "consistent"])):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(selected)

    scopes = {
        "original_retained_records": records[:source_records],
        "all_audited_records": records,
        "current_mvp_table_rsid_matches": [r for r in records if r["currently_used_rsid"] == "true"],
    }
    summary = {
        "audit": "Rounded source OR/CI/P fields checked against the Wald 95% CI relation under displayed-value rounding intervals.",
        "source_record_numbering": "1-based TSV data-record ordinal; header excluded and quoted embedded newlines parsed as one record.",
        "interpretation_limit": "A mismatch can reflect different models or other processing; this audit does not identify a root cause or establish source error or scientific validity.",
        "scope_limit": "Exposure-selected original rows plus current instrument rsIDs extracted from the complete GWAS where absent; not a whole-GWAS error-rate audit.",
        "current_unique_rsids": len(current_ids),
        "initially_covered_current_rsids": len(covered),
        "supplemented_records_from_full_GWAS":len(extra),
        "current_rsids_still_missing":sorted(current_ids-{r['original_rsid'] for r in records}),
        "inputs": {"source": SOURCE_REL, "source_sha256": digest(source), "current_table": CURRENT_REL, "current_table_sha256": digest(current)},
        "scopes": {},
        "outputs": {},
    }
    if extra:summary['inputs']['complete_GWAS_sha256']=digest(full_source)
    for label, subset in scopes.items():
        summary["scopes"][label] = {
            "retained_records": len(subset),
            "unique_nonempty_rsids": len({r["original_rsid"] for r in subset if r["original_rsid"]}),
            "consistent": sum(r["status"] == "consistent" for r in subset),
            "inconsistent": sum(r["status"] == "inconsistent" for r in subset),
            "unresolved": sum(r["status"] == "unresolved" for r in subset),
            "current_rsid_matches": sum(r["currently_used_rsid"] == "true" for r in subset),
            "inconsistent_rsid_examples": [r["original_rsid"] for r in subset if r["status"] == "inconsistent"][:20],
        }
    for path, count in ((audit_path, len(records)), (anomaly_path, len([r for r in records if r["status"] != "consistent"]))):
        summary["outputs"][path.name] = {"sha256": digest(path), "data_records": count}
    summary_path = destination / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Audited {len(records)} source records; {summary['scopes']['current_mvp_table_rsid_matches']['retained_records']} match rsIDs in the current MVP table.")
    print(f"Output: {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
