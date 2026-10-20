"""Apply the paper's explicit cell QC thresholds to the selected GEO matrices.

Only QC and lossless sparse conversion are performed. No clustering, annotation,
differential-expression inference or target selection is performed here.
"""
from pathlib import Path
import csv
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import scipy
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def qc_mask(features, counts, percent_mito, settings):
    return ((features >= settings["min_genes"]) & (features <= settings["max_genes"])
            & (counts >= settings["min_counts"]) & (counts <= settings["max_counts"])
            & (percent_mito <= settings["max_percent_mito"]))


def main(config_path):
    started = datetime.now(timezone.utc).isoformat()
    out = ROOT / "research/results/sc_qc"
    processed = ROOT / "research/data/processed/sc_qc"
    out.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    code_hash = digest(Path(__file__))
    (out / "run.json").write_text(json.dumps({"status": "running", "started_utc": started,
                                             "script_sha256": code_hash}, indent=2) + "\n")
    config = json.loads(config_path.read_text())
    mapping_path = ROOT / "research/data/metadata/ensembl_mt_genes.json"
    mapping = json.loads(mapping_path.read_text())
    mt_genes = {r["gene_id"] if "gene_id" in r else r["id"] for r in mapping
                if r.get("seq_region_name") == "MT" and r.get("assembly_name") == "GRCh38"}
    if not mt_genes:
        raise ValueError("No confirmed GRCh38 mitochondrial genes in downloaded mapping")
    manifest_path = ROOT / "research/results/geo_audit/sc_selected_samples.tsv"
    samples = list(csv.DictReader(manifest_path.open(), delimiter="\t"))
    summaries, inputs, all_qc = [], [], []
    generated = []
    global_features = set()
    for sample in samples:
        gsm = sample["gsm"]
        source = ROOT / "research/data/raw/GSE173706" / sample["supplementary_url"].rsplit("/", 1)[-1]
        totals = features = mito = None
        gene_ids, matrices, barcodes = [], [], None
        for block in pd.read_csv(source, index_col=0, chunksize=256):
            ids = block.index.astype(str).tolist()
            values = block.to_numpy()
            if np.any(values < 0) or not np.all(np.isfinite(values)) or not np.all(values == np.floor(values)):
                raise ValueError(f"Non-count values in {gsm}")
            if totals is None:
                barcodes = block.columns.astype(str).tolist()
                totals = np.zeros(len(barcodes), dtype=np.int64)
                features = np.zeros_like(totals)
                mito = np.zeros_like(totals)
            totals += values.sum(axis=0).astype(np.int64)
            features += np.count_nonzero(values, axis=0)
            is_mt = [g.split(".")[0] in mt_genes for g in ids]
            mito += values[is_mt].sum(axis=0).astype(np.int64)
            gene_ids.extend(ids)
            matrices.append(sparse.csr_matrix(values, dtype=np.int32))
        if len(set(gene_ids)) != len(gene_ids) or len(set(barcodes)) != len(barcodes):
            raise ValueError(f"Duplicate feature/cell identifiers in {gsm}")
        matched_mt = len({g.split(".")[0] for g in gene_ids} & mt_genes)
        if not matched_mt:
            raise ValueError(f"No mitochondrial feature matches in {gsm}")
        pct = np.divide(100.0 * mito, totals, out=np.full(len(totals), np.inf), where=totals > 0)
        keep = qc_mask(features, totals, pct, config["sc_qc"])
        qc = pd.DataFrame({"gsm": gsm, "donor_id": sample["donor_id"],
                           "tissue_type": sample["tissue_type"], "barcode": barcodes,
                           "n_features": features, "n_counts": totals,
                           "percent_mito": pct, "keep": keep})
        all_qc.append(qc)
        matrix = sparse.vstack(matrices, format="csr")[:, keep]
        # Save the original feature rows; no gene filtering beyond removal of zero
        # rows across all samples is assumed by this preprocessing checkpoint.
        sparse.save_npz(processed / f"{gsm}.npz", matrix)
        (processed / f"{gsm}.genes.txt").write_text("\n".join(gene_ids) + "\n")
        (processed / f"{gsm}.cells.txt").write_text("\n".join(np.asarray(barcodes)[keep]) + "\n")
        generated.extend(processed / f"{gsm}{suffix}" for suffix in (".npz", ".genes.txt", ".cells.txt"))
        nonzero = np.asarray(matrix.getnnz(axis=1)).ravel() > 0
        global_features.update(np.asarray(gene_ids)[nonzero])
        summaries.append({"gsm": gsm, "donor_id": sample["donor_id"], "tissue_type": sample["tissue_type"],
                          "input_features": len(gene_ids), "input_cells": len(barcodes),
                          "retained_cells": int(keep.sum()), "matched_mito_genes": matched_mt,
                          "expressed_features_after_qc": int(nonzero.sum())})
        inputs.append({"path": str(source.relative_to(ROOT)), "sha256": digest(source)})
        print(f"{gsm}: {len(barcodes)} -> {int(keep.sum())} cells; {matched_mt} mitochondrial genes", flush=True)
    pd.concat(all_qc, ignore_index=True).to_csv(out / "cell_qc.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(out / "sample_qc.tsv", sep="\t", index=False)
    generated.extend([out / "cell_qc.tsv.gz", out / "sample_qc.tsv"])
    result = {"status": "qc_completed_not_single_cell_reproduction", "samples": len(samples),
              "input_cells": sum(s["input_cells"] for s in summaries),
              "retained_cells": sum(s["retained_cells"] for s in summaries),
              "expressed_feature_union": len(global_features),
              "paper_reported_cells": 58867, "paper_reported_features": 23128,
              "mitochondrial_annotation": "Current Ensembl GRCh38 MT gene IDs; historical CellRanger annotation not provided",
              "unresolved": ["historical annotation release", "whether paper totals are before or after QC",
                             "clustering parameters and annotation mapping", "DE testing and correction families"],
              "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
              "settings": config["sc_qc"], "input_hashes": inputs,
              "metadata_hashes": {str(p.relative_to(ROOT)): digest(p) for p in [manifest_path, config_path, mapping_path]},
              "script_sha256": code_hash,
              "output_hashes": {str(p.relative_to(ROOT)): digest(p) for p in generated},
              "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__}
    (out / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k:result[k] for k in ["status","samples","input_cells","retained_cells","expressed_feature_union"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=ROOT / "research/config/analysis.json")
    arguments = parser.parse_args()
    try:
        main(arguments.config)
    except Exception as exc:
        record = ROOT / "research/results/sc_qc/run.json"
        if record.exists():
            failure = json.loads(record.read_text())
            failure.update(status="failed", error=str(exc), finished_utc=datetime.now(timezone.utc).isoformat())
            record.write_text(json.dumps(failure, indent=2) + "\n")
        raise
