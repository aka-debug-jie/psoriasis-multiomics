#!/usr/bin/env python3
"""Merge existing QC-filtered single-cell matrices without changing cell selection."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_cell_qc(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_unique_ensembl_symbols(path):
    symbols_by_ensembl = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("#tax_id") != "9606":
                continue
            symbol = row.get("Symbol", "").strip()
            if not symbol or symbol == "-":
                continue
            for xref in row.get("dbXrefs", "").split("|"):
                if xref.startswith("Ensembl:"):
                    ensembl = xref.removeprefix("Ensembl:").strip()
                    if ensembl:
                        symbols_by_ensembl.setdefault(ensembl, set()).add(symbol)
    return {
        ensembl: (next(iter(symbols)), "unique") if len(symbols) == 1 else ("", "ambiguous")
        for ensembl, symbols in symbols_by_ensembl.items()
    }


def _validate_integer_counts(matrix, gsm):
    values = matrix.data
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError(f"{gsm}: matrix contains non-finite or negative counts")
    if not np.all(values == np.floor(values)):
        raise ValueError(f"{gsm}: matrix contains non-integer counts")
    matrix = matrix.astype(np.int64, copy=True).tocsr()
    matrix.eliminate_zeros()
    return matrix


def _load_sample(processed_dir, sample, expected_cells):
    gsm = sample["gsm"]
    prefix = processed_dir / gsm
    matrix_path = prefix.with_suffix(".npz")
    genes_path = processed_dir / f"{gsm}.genes.txt"
    cells_path = processed_dir / f"{gsm}.cells.txt"
    for path in (matrix_path, genes_path, cells_path):
        if not path.is_file():
            raise FileNotFoundError(f"{gsm}: required QC artifact missing: {path}")

    matrix = _validate_integer_counts(sparse.load_npz(matrix_path), gsm)
    genes = genes_path.read_text(encoding="utf-8").splitlines()
    cells = cells_path.read_text(encoding="utf-8").splitlines()
    if matrix.shape != (len(genes), len(cells)):
        raise ValueError(f"{gsm}: matrix shape {matrix.shape} does not match {len(genes)} genes and {len(cells)} cells")
    if len(set(genes)) != len(genes):
        raise ValueError(f"{gsm}: duplicate gene identifiers")
    if len(set(cells)) != len(cells):
        raise ValueError(f"{gsm}: duplicate cell barcodes")
    if cells != [row["barcode"] for row in expected_cells]:
        raise ValueError(f"{gsm}: cells.txt order does not match cell_qc.tsv.gz keep==True rows")
    if int(sample["retained_cells"]) != len(cells):
        raise ValueError(f"{gsm}: sample_qc retained_cells={sample['retained_cells']} but found {len(cells)} cells")

    counts = np.asarray(matrix.sum(axis=0), dtype=np.int64).ravel()
    features = np.asarray(matrix.getnnz(axis=0), dtype=np.int64).ravel()
    for index, row in enumerate(expected_cells):
        if counts[index] != int(row["n_counts"]):
            raise ValueError(f"{gsm}/{cells[index]}: matrix count {counts[index]} != cell QC n_counts {row['n_counts']}")
        if features[index] != int(row["n_features"]):
            raise ValueError(f"{gsm}/{cells[index]}: matrix features {features[index]} != cell QC n_features {row['n_features']}")
    input_paths = (matrix_path, genes_path, cells_path)
    return genes, cells, matrix, [
        {"path": str(path.relative_to(processed_dir.parents[3])), "sha256": sha256(path)}
        for path in input_paths
    ]


def prepare(root):
    root = Path(root).resolve()
    qc_dir = root / "research/results/sc_qc"
    processed_dir = root / "research/data/processed/sc_qc"
    gene_info_path = root / "research/data/metadata/Homo_sapiens.gene_info.gz"
    sample_qc_path = qc_dir / "sample_qc.tsv"
    cell_qc_path = qc_dir / "cell_qc.tsv.gz"
    samples = read_tsv(sample_qc_path)
    all_cells = read_cell_qc(cell_qc_path)
    if not samples:
        raise ValueError("sample_qc.tsv contains no samples")
    sample_order = [row["gsm"] for row in samples]
    if len(sample_order) != len(set(sample_order)):
        raise ValueError("sample_qc.tsv contains duplicate GSM identifiers")
    if any(row["gsm"] not in sample_order for row in all_cells):
        raise ValueError("cell_qc.tsv.gz contains a GSM absent from sample_qc.tsv")
    qc_sample_order = list(dict.fromkeys(row["gsm"] for row in all_cells))
    if qc_sample_order != sample_order:
        raise ValueError(f"cell_qc sample order {qc_sample_order} != sample_qc order {sample_order}")
    if any(row.get("keep") not in ("True", "False") for row in all_cells):
        raise ValueError("cell_qc.tsv.gz keep values must be True or False")

    keep_by_gsm = {
        gsm: [row for row in all_cells if row["gsm"] == gsm and row["keep"] == "True"]
        for gsm in sample_order
    }
    matrices = []
    input_hashes = []
    gene_order = None
    cell_rows = []
    for sample in samples:
        gsm = sample["gsm"]
        expected_cells = keep_by_gsm[gsm]
        for row in expected_cells:
            if row["donor_id"] != sample["donor_id"] or row["tissue_type"] != sample["tissue_type"]:
                raise ValueError(f"{gsm}/{row['barcode']}: donor/tissue metadata disagrees with sample_qc.tsv")
        genes, barcodes, matrix, sample_hashes = _load_sample(processed_dir, sample, expected_cells)
        if gene_order is None:
            gene_order = genes
        elif genes != gene_order:
            raise ValueError(f"{gsm}: gene order differs from the first sample; refusing to hstack")
        matrices.append(matrix)
        input_hashes.extend(sample_hashes)
        for row in expected_cells:
            cell_rows.append({
                "cell_id": f"{gsm}_{row['barcode']}",
                "gsm": gsm,
                "barcode": row["barcode"],
                "donor_id": row["donor_id"],
                "tissue_type": row["tissue_type"],
                "n_features": row["n_features"],
                "n_counts": row["n_counts"],
                "percent_mito": row["percent_mito"],
            })
    cell_ids = [row["cell_id"] for row in cell_rows]
    if len(cell_ids) != len(set(cell_ids)):
        raise ValueError("Constructed cell IDs are not unique")

    merged = sparse.hstack(matrices, format="csr", dtype=np.int64)
    nonzero_rows = np.asarray(merged.getnnz(axis=1)).ravel() > 0
    retained_genes = [gene for gene, keep in zip(gene_order, nonzero_rows) if keep]
    merged = merged[nonzero_rows, :].tocsr()
    symbols = read_unique_ensembl_symbols(gene_info_path)

    output_dir = root / "research/data/processed/sc_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = output_dir / "counts.mtx.gz"
    feature_path = output_dir / "features.tsv"
    cell_path = output_dir / "cells.tsv"
    with gzip.open(matrix_path, "wb") as handle:
        mmwrite(handle, merged)
    with feature_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_id", "symbol", "mapping_status"))
        writer.writerows((gene, *symbols.get(gene, ("", "unmapped"))) for gene in retained_genes)
    cell_fields = ("cell_id", "gsm", "barcode", "donor_id", "tissue_type", "n_features", "n_counts", "percent_mito")
    with cell_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cell_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(cell_rows)

    metadata_inputs = (sample_qc_path, cell_qc_path, gene_info_path)
    input_hashes.extend({"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in metadata_inputs)
    script_path = Path(__file__).resolve()
    result = {
        "status": "merged_qc_selected_cells",
        "samples": sample_order,
        "cells": len(cell_rows),
        "genes_before_global_zero_filter": len(gene_order),
        "genes_after_global_zero_filter": len(retained_genes),
        "global_zero_genes_removed": len(gene_order) - len(retained_genes),
        "matrix_shape": list(merged.shape),
        "matrix_format": "gzip-compressed Matrix Market coordinate",
        "selection_source": "research/results/sc_qc/cell_qc.tsv.gz rows with keep == True, in recorded order",
        "gene_symbol_mapping": "Homo sapiens gene_info dbXrefs Ensembl mapping; only a single distinct Symbol per ENSG is labeled; ENSG IDs are retained and rows are never merged",
        "input_sha256": input_hashes,
        "script_sha256": sha256(script_path),
        "output_sha256": {
            path.name: sha256(path) for path in (matrix_path, feature_path, cell_path)
        },
        "outputs": [str(path.relative_to(root)) for path in (matrix_path, feature_path, cell_path)],
    }
    (output_dir / "run.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    print(json.dumps(prepare(args.root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
