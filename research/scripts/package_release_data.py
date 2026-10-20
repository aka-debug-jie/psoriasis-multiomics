#!/usr/bin/env python3
"""Copy selected result tables into a byte-preserving release data package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import tempfile
from pathlib import Path


FILES = {
    "results/mr_fenland_substitute/all_proteins_comparison.tsv": "mr_all.tsv",
    "results/mr_fenland_substitute/followup_proteins.tsv": "targets.tsv",
    "results/mr_fenland_substitute/mvp_cohort_wald.tsv": "mr_mvp_cohorts.tsv",
    "results/mr_fenland_substitute/finngen_cohort_wald.tsv": "mr_finngen_cohorts.tsv",
    "results/mr_fenland_substitute/stuart_cohort_wald.tsv": "mr_stuart_cohorts.tsv",
    "results/clinical_mr_substitute/all_meta_results.tsv": "clinical_mr.tsv",
    "results/coloc_mrna_substitute/coloc_summary.tsv": "coloc_mrna.tsv",
    "results/coloc_protein_substitute/coloc_summary.tsv": "coloc_protein.tsv",
    "results/multi_instrument_substitute/method_results.tsv": "sensitivity.tsv",
    "results/steiger_substitute/candidate_robustness.tsv": "steiger.tsv",
    "results/evidence_substitute/followup_evidence.tsv": "evidence.tsv",
    "results/evidence_substitute/expression_gene_identity.tsv": "gene_identity.tsv",
    "results/evidence_substitute/cell_Wilcoxon_followup_contrasts.tsv": "expression_cell.tsv",
    "results/evidence_substitute/donor_pseudobulk_followup_contrasts.tsv": "expression_donor.tsv",
    "results/bulk_ncbi/all_gene_differential_expression.tsv": "bulk_all.tsv",
    "results/drug_network_substitute/all_drug_target_relations.tsv": "drug_relations.tsv",
    "results/drug_network_substitute/candidate_network_drug_associations.tsv": "drug_network_links.tsv",
    "results/drug_network_substitute/STRING_network.tsv": "string_edges.tsv",
    "results/drug_network_substitute/STRING_mapping.tsv": "string_nodes.tsv",
    "results/functional_substitute/GO_biological_process.tsv": "functional/GO_biological_process.tsv",
    "results/functional_substitute/GO_cellular_component.tsv": "functional/GO_cellular_component.tsv",
    "results/functional_substitute/GO_molecular_function.tsv": "functional/GO_molecular_function.tsv",
    "results/functional_substitute/KEGG.tsv": "functional/KEGG.tsv",
    "data/processed/clinical_finngen/phenotype_manifest.tsv": "clinical_endpoints.tsv",
    "results/shared_outcome_audit/followup_comparison.tsv": "audit_single_iv_targets.tsv",
    "results/shared_outcome_audit/all_proteins_comparison.tsv": "audit_single_iv_all.tsv",
    "results/shared_outcome_audit/clinical_comparison.tsv": "audit_single_iv_clinical.tsv",
    "results/shared_outcome_audit/single_IV_frozen.tsv": "audit_single_iv_instruments.tsv",
    "results/mvp_precision_audit/mvp_precision_audit.tsv": "audit_mvp_ci.tsv",
    "results/mvp_precision_audit/anomalies.tsv": "audit_mvp_ci_anomalies.tsv",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_tsv(path: Path) -> tuple[int, int, list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"TSV is empty and has no header: {path}")
        counts = [0] * len(header)
        kinds = [set() for _ in header]
        rows = 0
        for row in reader:
            rows += 1
            if len(row) != len(header):
                raise ValueError(f"Malformed TSV row in {path}: expected {len(header)} columns, got {len(row)}")
            for i, value in enumerate(row):
                counts[i] += 1
                if value == "" or value.upper() == "NA":
                    continue
                if value.lower() in {"true", "false"}:
                    kinds[i].add("boolean")
                elif re.fullmatch(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?", value):
                    kinds[i].add("numeric")
                else:
                    kinds[i].add("string")
    columns = []
    for i, name in enumerate(header):
        observed = sorted(kinds[i])
        columns.append({
            "column_position": str(i + 1),
            "column": name,
            "empty_count": "0",
            "dtype_observed": observed[0] if len(observed) == 1 else ("empty" if not observed else "mixed:" + ",".join(observed)),
        })
    # Empty strings and literal NA are reported separately for transparent missingness counts.
    # Re-read the table to count empty cells without coercing values.
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        empty_counts = [0] * len(header)
        na_counts = [0] * len(header)
        for row in reader:
            for i, value in enumerate(row):
                empty_counts[i] += value == ""
                na_counts[i] += value.upper() == "NA"
    for i, column in enumerate(columns):
        column["empty_count"] = str(empty_counts[i])
        column["na_count"] = str(na_counts[i])
    return rows, len(header), header, columns


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2], help="project root (defaults to this checkout)")
    args = parser.parse_args()
    project_root = args.root.resolve()
    research_root = project_root / "research"
    output_root = research_root / "release" / "data"

    sources = {relative: research_root / relative for relative in FILES}
    missing = [str(path) for path in sources.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required release input(s) missing:\n" + "\n".join(missing))

    expected_outputs = {name for name in FILES.values()} | {"catalog.tsv", "dictionary.tsv", "manifest.json"}
    manifest_path = output_root / "manifest.json"
    owned = set()
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        owned = set(old.get("outputs", []))
    collisions = [name for name in expected_outputs if (output_root / name).exists() and name not in owned]
    if collisions:
        raise FileExistsError("Refusing to overwrite outputs not owned by the existing package manifest: " + ", ".join(sorted(collisions)))

    catalog_rows = []
    dictionary_rows = []
    for source_rel, output_rel in FILES.items():
        source = sources[source_rel]
        rows, cols, header, columns = inspect_tsv(source)
        catalog_rows.append({
            "output": output_rel,
            "source": source_rel,
            "source_sha256": sha256(source),
            "output_sha256": sha256(source),
            "rows": rows,
            "columns": cols,
            "column_names": json.dumps(header, ensure_ascii=False),
        })
        for column in columns:
            dictionary_rows.append({"output": output_rel, **column})

    output_root.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".release-data-", dir=output_root.parent))
    try:
        for source_rel, output_rel in FILES.items():
            destination = stage / output_rel
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(sources[source_rel], destination)
        for filename, rows, fields in [
            ("catalog.tsv", catalog_rows, ["output", "source", "source_sha256", "output_sha256", "rows", "columns", "column_names"]),
            ("dictionary.tsv", dictionary_rows, ["output", "column_position", "column", "empty_count", "na_count", "dtype_observed"]),
        ]:
            with (stage / filename).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        manifest = {
            "format": "release-data-package-v1",
            "outputs": sorted(expected_outputs),
            "files": catalog_rows,
        }
        (stage / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        # Replace only paths named by this package; staged files are complete before any destination changes.
        for name in sorted(expected_outputs):
            destination = output_root / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(stage / name, destination)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
    print(f"Packaged {len(FILES)} tables into {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
