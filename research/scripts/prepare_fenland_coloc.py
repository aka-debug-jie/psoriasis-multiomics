#!/usr/bin/env python3
"""Prepare frozen Fenland Olink cis inputs for the downstream coloc module."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "research/data/processed/coloc_protein_fenland"
REGIONS = ROOT / "research/data/processed/Fenland_Olink/reference_and_followup_regions.tsv.gz"
INDICES = ROOT / "research/data/processed/Fenland_Olink/source_indices.tsv"
TARGETS = ROOT / "research/data/processed/Fenland_Olink/target_registry.tsv"
GENES = ROOT / "research/data/processed/coloc_all_genes/gene_registry.tsv"
MVP = ROOT / "research/data/processed/smr_gwas/mvp_full.ma"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def write_tsv(path: Path, header, records):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=header, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def is_autosome(value: str) -> bool:
    try:
        return 1 <= int(value) <= 22
    except ValueError:
        return False


def build(input_paths: dict[str, Path], out: Path):
    genes = list(rows(input_paths["genes"]))
    gene_by_id = {r["gene_id"]: r for r in genes}
    targets = list(rows(input_paths["targets"]))
    by_id = {r["id"]: r for r in targets}
    frozen_indices = {r["id"]: r for r in rows(input_paths["indices"])}

    eligible = []
    exclusions = defaultdict(list)
    for target in targets:
        gid = target["ensembl_gene_id"]
        if gid not in gene_by_id:
            continue
        if target["assay_multi_protein"].lower() == "true":
            exclusions[gid].append((target, "composite_assay_excluded"))
        elif not is_autosome(target["chr"]):
            exclusions[gid].append((target, "non_autosomal_target"))
        else:
            eligible.append(target)

    max_logp = {}
    source_pairs = defaultdict(list)
    for row in rows(input_paths["regions"]):
        if row["id"] not in {t["id"] for t in eligible}:
            continue
        source_pairs[row["id"]].append(row)
        value = float(row["log10p"])
        max_logp[row["id"]] = max(max_logp.get(row["id"], -math.inf), value)

    # The source-region maximum (not the frozen index row) selects one assay per UniProt.
    chosen = {}
    for target in eligible:
        if target["id"] not in max_logp:
            continue
        u = target["Uniprot_ID"]
        rank = (-max_logp[target["id"]], target["id"])
        if u not in chosen or rank < chosen[u][0]:
            chosen[u] = (rank, target)

    registry = []
    selected_targets = {item[1]["id"]: item[1] for item in chosen.values()}
    for gene in genes:
        gid = gene["gene_id"]
        matched = [t for t in targets if t["ensembl_gene_id"] == gid]
        gene_targets = [t for t in eligible if t["ensembl_gene_id"] == gid]
        if not matched:
            status = "no_target_mapping"
            target = None
        elif not gene_targets:
            status = "composite_or_non_autosomal_only"
            target = None
        else:
            uni = {t["Uniprot_ID"] for t in gene_targets}
            chosen_for_gene = [selected_targets[t["id"]] for t in gene_targets if t["id"] in selected_targets]
            if not chosen_for_gene:
                status, target = "no_source_region", None
            elif len(uni) != 1:
                # Do not silently collapse distinct proteins under one gene identifier.
                status, target = "multiple_uniprot_mappings", None
            else:
                target = next((t for t in chosen_for_gene if t["Uniprot_ID"] == next(iter(uni))), None)
                status = "selected" if target else "assay_not_selected_for_uniprot"
        rec = {"gene_id": gid, "symbol": gene["symbol"], "uniprot": "", "assay_id": "", "assay": "",
               "source_max_log10p": "", "frozen_index_rsid": "", "frozen_index_log10p": "",
               "sdY_estimated": "", "status": status, "excluded_assays": ""}
        if target:
            idx = frozen_indices.get(target["id"])
            rec.update(uniprot=target["Uniprot_ID"], assay_id=target["id"], assay=target["Assay"],
                       source_max_log10p=max_logp[target["id"]],
                       frozen_index_rsid=idx["rsid"] if idx else "",
                       frozen_index_log10p=idx["log10p"] if idx else "",
                       sdY_estimated=idx["sdY_estimated"] if idx else "",
                       status="selected" if idx else "selected_missing_frozen_index")
        excluded = [t["id"] + ":" + reason for t, reason in exclusions[gid]]
        rec["excluded_assays"] = ";".join(excluded)
        registry.append(rec)

    chosen_ids = set(selected_targets)
    raw_kept = []
    all_qc_rows = 0
    for identifier in chosen_ids:
        target = by_id[identifier]
        idx = frozen_indices.get(identifier)
        if not idx:
            continue
        sd = float(idx["sdY_estimated"])
        if not math.isfinite(sd) or sd <= 0:
            raise ValueError(f"invalid sdY_estimated for frozen index {identifier}: {sd}")
        for row in source_pairs.get(identifier, []):
            se, freq, info = float(row["StdErr"]), float(row["Freq1"]), float(row["info"])
            if not (se > 0 and info >= 0.4 and 0.001 <= freq <= 0.999):
                continue
            all_qc_rows += 1
            raw_kept.append({"SNP": row["rsid"], "chr": row["chr"], "pos": row["pos"],
                             "A1": row["Allele1"], "A2": row["Allele2"],
                             "beta": float(row["Effect"]) / sd, "se": se / sd,
                             "freq": freq, "assay": target["Assay"], "gene_id": target["ensembl_gene_id"],
                             "uniprot": target["Uniprot_ID"], "raw_beta": float(row["Effect"]),
                             "raw_se": se, "assay_id": identifier, "sdY_estimated": sd})

    # A variant ID is retained only when it has one unambiguous coordinate/allele identity.
    identities = defaultdict(set)
    for row in raw_kept:
        identities[row["SNP"]].add((row["chr"], row["pos"], row["A1"], row["A2"]))
    conflict_rsids = {rsid for rsid, keys in identities.items() if len(keys) > 1}
    final = [r for r in raw_kept if r["SNP"] not in conflict_rsids]
    out.mkdir(parents=True, exist_ok=True)
    registry_path = out / "index_registry.tsv"
    write_tsv(registry_path, list(registry[0]) if registry else [], registry)
    cis_path = out / "full_cis.tsv.gz"
    cis_fields = ["SNP", "chr", "pos", "A1", "A2", "beta", "se", "freq", "assay", "gene_id", "uniprot", "raw_beta", "raw_se", "assay_id", "sdY_estimated"]
    with gzip.open(cis_path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=cis_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(final)

    wanted = {r["SNP"] for r in final}
    mvp_path = out / "mvp_matched.ma.gz"
    mvp_rows = 0
    with input_paths["mvp"].open("r", encoding="utf-8", newline="") as src, gzip.open(mvp_path, "wt", encoding="utf-8", newline="") as dst:
        header = src.readline()
        columns = header.rstrip("\r\n").split("\t")
        snp_col = columns.index("SNP")
        dst.write(header)
        for line in src:
            fields = line.rstrip("\r\n").split("\t")
            if fields[snp_col] in wanted:
                dst.write(line); mvp_rows += 1
    observed = set()
    with gzip.open(mvp_path, "rt", encoding="utf-8", newline="") as stream:
        next(stream)
        for line in stream:
            observed.add(line.split("\t", 1)[0])
    missing_mvp = sorted(wanted - observed)
    missing_path = out / "mvp_missing_ids.tsv"
    write_tsv(missing_path, ["SNP"], ({"SNP": value} for value in missing_mvp))
    record = {
        "status": "prepared_input_adapter_only", "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_rule": "maximum full-source log10p per UniProt; ties by assay_id lexical order",
        "source_index_rule": "frozen source_indices.tsv exact assay_id only; no replacement",
        "pair_filter": "StdErr>0; info>=0.4; Freq1 in [0.001,0.999]; conflicting rsID groups excluded",
        "scale": "beta and se divided by frozen sdY_estimated; raw_beta/raw_se retained",
        "gene_registry_rows": len(genes), "selected_assays": len(chosen_ids),
        "selected_assays_missing_frozen_index": sum(r["status"] == "selected_missing_frozen_index" for r in registry),
        "quality_pass_before_conflict": all_qc_rows, "conflict_rsids_excluded": len(conflict_rsids),
        "full_cis_rows": len(final), "matched_mvp_rows": mvp_rows,
        "matched_mvp_unique_rsids": len(observed), "requested_ids_missing_from_mvp": len(missing_mvp),
        "missing_id_examples": missing_mvp[:20],
        "input_sha256": {p.name: sha256(p) for p in input_paths.values()},
        "script_sha256": sha256(Path(__file__)),
        "output_sha256": {p.name: sha256(p) for p in (registry_path, cis_path, mvp_path, missing_path)},
    }
    (out / "run.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", type=Path, default=GENES)
    parser.add_argument("--targets", type=Path, default=TARGETS)
    parser.add_argument("--indices", type=Path, default=INDICES)
    parser.add_argument("--regions", type=Path, default=REGIONS)
    parser.add_argument("--mvp", type=Path, default=MVP)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    record = build({k: getattr(args, k) for k in ("genes", "targets", "indices", "regions", "mvp")}, args.out)
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
