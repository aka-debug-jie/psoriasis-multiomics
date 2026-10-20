#!/usr/bin/env python3
"""Normalize the selected follow-up protein pQTL loci and stream matching MVP rows."""

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri_exp


ROOT = Path(__file__).resolve().parents[2]
SOURCE_SPECS = (
    {
        "cohort": "Iceland",
        "filename": "ST19_somascan_pQTLs_smpnorm.tsv",
        "assay": "SeqId",
        "gene": "gene prot",
        "uniprot": "UniProt",
        "locus_type": "cis trans",
        "variant": "variant",
        "chr": "chr var",
        "pos": "pos var",
        "ea": "Amin",
        "oa": "Amaj",
        "freq_pct": "MAF PC",
        "beta": "beta unadj",
        "neg_log10_p": "mLog10pval unadj",
    },
    {
        "cohort": "UKB-PPP",
        "filename": "ST21_olink_ukb_bi_pQTLs.tsv",
        "assay": "OlinkID",
        "gene": "Affected Protein Gene name",
        "uniprot": "UniProt",
        "locus_type": "CisOrTrans",
        "variant": "rsName",
        "chr": "Chrom",
        "pos": "Pos",
        "ea": "Amin",
        "oa": "Amaj",
        "freq_pct": "MAF PC",
        "beta": "Effect Amin unadj",
        "neg_log10_p": "Log10 pval gc cor unadj",
    },
)
FOLLOWUP_PATH = "research/results/available_followup/followup_proteins.tsv"
COHORT_WALD_PATH = "research/results/mr_available/cohort_wald.tsv"
MVP_PATH = "research/data/processed/smr_gwas/mvp_full.ma"
SE_SOURCE = "normal_inversion_of_original_LR_p_value;matches_author_S2_rounding;not_reported_SE"


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _numeric_equal(source_value, target_value, label, context):
    try:
        left = float(source_value)
        right = float(target_value)
    except (TypeError, ValueError):
        raise ValueError(f"{context}: nonnumeric {label}: source={source_value!r}, cohort_wald={target_value!r}")
    if not math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError(f"{context}: {label} differs: source={left}, cohort_wald={right}")


def _validate_selected_wald_row(source, source_spec, wald_row):
    context = f"{wald_row['cohort']}/{wald_row['uniprot']}/{wald_row['variant']}"
    _numeric_equal(source[source_spec["beta"]], wald_row["beta.exposure"], "beta", context)
    _numeric_equal(source[source_spec["neg_log10_p"]], wald_row["neg_log10_p"], "neg_log10_p", context)
    _numeric_equal(float(source[source_spec["freq_pct"]]) / 100.0, wald_row["eaf.exposure"], "eaf", context)
    _numeric_equal(source[source_spec["pos"]], wald_row["pos.exposure"], "position", context)
    source_chr = str(source[source_spec["chr"]]).replace("chr", "", 1)
    wald_chr = str(wald_row["chr.exposure"]).replace("chr", "", 1)
    if source_chr != wald_chr:
        raise ValueError(f"{context}: chromosome differs: source={source_chr!r}, cohort_wald={wald_chr!r}")
    for field, column in ((source_spec["ea"], "EA.exposure"), (source_spec["oa"], "OA.exposure")):
        if str(source[field]) != str(wald_row[column]):
            raise ValueError(f"{context}: {field} allele differs: source={source[field]!r}, cohort_wald={wald_row[column]!r}")
    if str(source[source_spec["gene"]]) != str(wald_row["gene"]):
        raise ValueError(f"{context}: gene label differs between source and selected cohort_wald row")


def _source_frame(root, spec, followup_proteins, selected_pairs, pair_rows):
    path = root / "research/data/processed/exposures" / spec["filename"]
    frame = pd.read_csv(path, sep="\t", low_memory=False, dtype={spec["uniprot"]: str})
    required = {spec[field] for field in (
        "assay", "gene", "uniprot", "locus_type", "variant", "chr", "pos", "ea", "oa", "freq_pct", "beta", "neg_log10_p"
    )}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{path}: missing mapped source fields {missing}")
    frame["source_data_row"] = np.arange(1, len(frame) + 1, dtype=np.int64)
    frame[spec["assay"]] = frame[spec["assay"]].astype(str)
    frame[spec["uniprot"]] = frame[spec["uniprot"]].astype(str)
    selected_keys = {(uniprot, assay) for (cohort, uniprot), assay in selected_pairs.items() if cohort == spec["cohort"]}
    mask = [(uniprot, assay) in selected_keys for uniprot, assay in zip(frame[spec["uniprot"]], frame[spec["assay"]])]
    selected = frame.loc[mask].copy()
    selected_locus_types = set(selected[spec["locus_type"]].dropna().astype(str))
    unexpected_types = selected_locus_types - {"cis", "trans"}
    if unexpected_types:
        raise ValueError(f"{spec['cohort']}: selected assay rows have unexpected cis/trans labels: {sorted(unexpected_types)}")

    missing_source_pairs = []
    for key in selected_keys:
        if not (((selected[spec["uniprot"]] == key[0]) & (selected[spec["assay"]] == key[1])).any()):
            missing_source_pairs.append(key)
    if missing_source_pairs:
        raise ValueError(f"{spec['cohort']}: selected cohort_wald assays missing from source table: {missing_source_pairs}")

    # Confirm the cohort_wald-selected instrument maps back to exactly one original source row.
    for key, wald_row in pair_rows.items():
        if key[0] != spec["cohort"]:
            continue
        matches = selected[
            (selected[spec["uniprot"]] == key[1])
            & (selected[spec["assay"]] == str(wald_row["assay"]))
            & (selected[spec["variant"]].astype(str) == str(wald_row["variant"]))
        ]
        if len(matches) != 1:
            raise ValueError(f"{key[0]}/{key[1]}: cohort_wald variant maps to {len(matches)} original source rows")
        _validate_selected_wald_row(matches.iloc[0], spec, wald_row)

    beta = pd.to_numeric(selected[spec["beta"]], errors="coerce").to_numpy(dtype=float)
    logp = pd.to_numeric(selected[spec["neg_log10_p"]], errors="coerce").to_numpy(dtype=float)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        z_score = -ndtri_exp(-logp * np.log(10.0) - np.log(2.0))
        se = np.abs(beta) / z_score
    normalized_chr = selected[spec["chr"]].astype(str).str.replace("chr", "", regex=False)
    normalized = pd.DataFrame({
        "cohort": spec["cohort"],
        "assay": selected[spec["assay"]].astype(str),
        "gene": selected[spec["gene"]],
        "uniprot": selected[spec["uniprot"]].astype(str),
        "locus_type": selected[spec["locus_type"]].astype(str),
        "variant": selected[spec["variant"]].astype(str),
        "chr": normalized_chr,
        "pos": selected[spec["pos"]],
        "EA": selected[spec["ea"]],
        "OA": selected[spec["oa"]],
        "eaf": pd.to_numeric(selected[spec["freq_pct"]], errors="coerce") / 100.0,
        "beta": beta,
        "neg_log10_p": logp,
        "se_reconstructed": se,
        "se_source": SE_SOURCE,
        "source_table": spec["filename"],
        "source_data_row": selected["source_data_row"].to_numpy(dtype=np.int64),
        "cohort_wald_exposure_row_id": [
            pair_rows[(spec["cohort"], uniprot)]["exposure_row_id"]
            for uniprot in selected[spec["uniprot"]]
        ],
        "input_role": "followup_protein_selected_assay_all_cis_trans_loci_no_p_filter",
    })
    return normalized, path, len(frame)


def explicit_rsids(variants):
    ids = set()
    for value in variants:
        if pd.isna(value):
            continue
        for token in str(value).split(","):
            candidate = token.strip()
            if re.fullmatch(r"rs[0-9]+", candidate):
                ids.add(candidate)
    return ids


def stream_mvp_subset(source_path, output_path, requested_ids):
    digest = hashlib.sha256()
    matched_ids = set()
    matched_rows = 0
    data_rows = 0
    with source_path.open("rb") as source, output_path.open("wb") as output:
        header = source.readline()
        if not header:
            raise ValueError(f"Empty MVP summary file: {source_path}")
        digest.update(header)
        header_fields = header.rstrip(b"\r\n").decode("utf-8").split("\t")
        if "SNP" not in header_fields:
            raise ValueError(f"MVP header lacks explicit SNP column: {header_fields[:12]}")
        snp_index = header_fields.index("SNP")
        output.write(header)
        for line_number, raw_line in enumerate(source, start=2):
            digest.update(raw_line)
            data_rows += 1
            fields = raw_line.rstrip(b"\r\n").split(b"\t")
            if len(fields) <= snp_index:
                raise ValueError(f"Malformed MVP row at line {line_number}: expected SNP field index {snp_index}")
            snp = fields[snp_index].decode("utf-8")
            if snp in requested_ids:
                output.write(raw_line)
                matched_ids.add(snp)
                matched_rows += 1
    return {
        "sha256": digest.hexdigest(),
        "data_rows_scanned": data_rows,
        "matched_rows": matched_rows,
        "matched_unique_rsids": sorted(matched_ids),
        "missing_rsids": sorted(requested_ids - matched_ids),
    }


def prepare(root):
    root = Path(root).resolve()
    followup_path = root / FOLLOWUP_PATH
    cohort_wald_path = root / COHORT_WALD_PATH
    followup = pd.read_csv(followup_path, sep="\t", dtype={"uniprot": str})
    if "uniprot" not in followup.columns:
        raise ValueError(f"{followup_path}: missing uniprot column")
    followup_uniprots = set(followup["uniprot"].dropna().astype(str))
    if not followup_uniprots:
        raise ValueError("followup_proteins.tsv has no UniProt identifiers")
    if followup["uniprot"].dropna().astype(str).duplicated().any():
        raise ValueError("followup_proteins.tsv contains duplicate UniProt identifiers")

    wald = pd.read_csv(cohort_wald_path, sep="\t", low_memory=False, dtype={
        "cohort": str, "assay": str, "gene": str, "uniprot": str, "variant": str, "exposure_row_id": str,
    })
    required_wald = {
        "cohort", "assay", "gene", "uniprot", "variant", "chr.exposure", "pos.exposure",
        "EA.exposure", "OA.exposure", "eaf.exposure", "beta.exposure", "neg_log10_p", "exposure_row_id",
    }
    missing_wald = sorted(required_wald - set(wald.columns))
    if missing_wald:
        raise ValueError(f"{cohort_wald_path}: missing columns {missing_wald}")
    selected = wald[wald["uniprot"].isin(followup_uniprots)].copy()
    pair_counts = selected.groupby(["cohort", "uniprot"], dropna=False).size()
    duplicate_pairs = pair_counts[pair_counts != 1]
    if len(duplicate_pairs):
        raise ValueError(f"cohort_wald must contain one fixed selected row per (cohort, UniProt): {duplicate_pairs.to_dict()}")
    selected_pairs = {}
    pair_rows = {}
    for _, row in selected.iterrows():
        key = (str(row["cohort"]), str(row["uniprot"]))
        selected_pairs[key] = str(row["assay"])
        pair_rows[key] = row.to_dict()

    normalized_tables = []
    source_hashes = {str(followup_path.relative_to(root)): sha256_file(followup_path),
                     str(cohort_wald_path.relative_to(root)): sha256_file(cohort_wald_path)}
    source_rows = {}
    for spec in SOURCE_SPECS:
        normalized, path, total_rows = _source_frame(root, spec, followup_uniprots, selected_pairs, pair_rows)
        normalized_tables.append(normalized)
        source_hashes[str(path.relative_to(root))] = sha256_file(path)
        source_rows[spec["cohort"]] = {"source_rows": total_rows, "selected_loci": len(normalized)}

    all_loci = pd.concat(normalized_tables, ignore_index=True)
    if all_loci.empty:
        raise ValueError("No ST19/ST21 loci matched followup UniProt and cohort_wald-selected assay")
    if not set(all_loci["locus_type"]) <= {"cis", "trans"}:
        raise ValueError("Normalized loci contain labels outside cis/trans")

    output_dir = root / "research/data/processed/sensitivity_inputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    all_loci_path = output_dir / "all_loci.tsv"
    mvp_subset_path = output_dir / "mvp_subset.tsv"
    all_loci.to_csv(all_loci_path, sep="\t", index=False, na_rep="NA")

    requested_ids = explicit_rsids(all_loci["variant"])
    mvp_path = root / MVP_PATH
    mvp_report = stream_mvp_subset(mvp_path, mvp_subset_path, requested_ids)
    source_hashes[str(mvp_path.relative_to(root))] = mvp_report.pop("sha256")

    script_path = Path(__file__).resolve()
    result = {
        "status": "prepared_unselected_sensitivity_inputs",
        "scope": "Published ST19/ST21 cis/trans locus rows for followup UniProt and fixed cohort_wald assays only; not a full GWAS and no p/F/LD/direction filtering.",
        "se_reconstruction": "abs(beta) / -ndtri_exp(-neg_log10_p*ln(10)-ln(2)); a numerical inversion, not an original reported SE",
        "followup_proteins": len(followup_uniprots),
        "cohort_wald_selected_assay_pairs": len(selected_pairs),
        "source_row_counts": source_rows,
        "all_loci_rows": len(all_loci),
        "all_loci_rows_by_cohort": {str(k): int(v) for k, v in all_loci["cohort"].value_counts().items()},
        "all_loci_rows_by_type": {str(k): int(v) for k, v in all_loci["locus_type"].value_counts().items()},
        "explicit_variant_rsids": len(requested_ids),
        "mvp_subset": mvp_report,
        "input_sha256": source_hashes,
        "script_sha256": sha256_file(script_path),
        "output_sha256": {
            all_loci_path.name: sha256_file(all_loci_path),
            mvp_subset_path.name: sha256_file(mvp_subset_path),
        },
        "outputs": [str(all_loci_path.relative_to(root)), str(mvp_subset_path.relative_to(root))],
    }
    records_dir = root / "research/records/advance_20261017"
    records_dir.mkdir(parents=True, exist_ok=True)
    (records_dir / "sensitivity_inputs_run.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    print(json.dumps(prepare(args.root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
