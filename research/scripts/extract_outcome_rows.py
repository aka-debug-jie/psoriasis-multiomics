"""Stream complete GWAS files and retain explicit rsIDs present in cis inputs.

Coordinates are retained for inspection but are not used to join different builds.
"""
from pathlib import Path
import csv
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main(selected):
    exposure = ROOT / "research/data/processed/exposures/eldjarn_all_cis_unselected.tsv"
    x = pd.read_csv(exposure, sep="\t")
    rsids = {token for v in x.variant.dropna().astype(str) for token in v.split(",") if token.startswith("rs")}
    out = ROOT / "research/data/processed/outcomes"
    out.mkdir(parents=True, exist_ok=True)
    record_path = ROOT / "research/records/outcome_extraction.json"
    previous = json.loads(record_path.read_text()) if record_path.exists() else {}
    inventories = [x for x in previous.get("outcomes", []) if x["outcome"] not in selected]
    for item in inventories:
        item.setdefault("script_sha256", previous.get("script_sha256"))
    specs = [("mvp", "GCST90476190.tsv.gz", "rsid"),
             ("finngen_r12", "finngen_R12_L12_PSORIASIS.gz", "rsids"),
             ("stuart", "34927100-GCST90019016-EFO_0000676.h.tsv.gz", "hm_rsid")]
    for name, filename, rs_col in specs:
        if name not in selected:
            continue
        source = ROOT / "research/data/raw/gwas" / filename
        if not source.exists():
            raise FileNotFoundError(f"Complete, checksum-verified input required: {source}")
        parts, count = [], 0
        for chunk in pd.read_csv(source, sep="\t", chunksize=100000, dtype=str, na_values=["#NA"], keep_default_na=True):
            count += len(chunk)
            # rsID aliases are explicit in source files; no position-based proxy inference.
            ids = chunk[rs_col].fillna("")
            mask = ids.isin(rsids)
            aliased = ids.str.contains("[;,]", regex=True)
            if aliased.any():
                mask |= ids.where(aliased, "").map(lambda s: any(t in rsids for t in s.replace(";", ",").split(",")))
            if mask.any():
                parts.append(chunk.loc[mask].copy())
        retained = pd.concat(parts, ignore_index=True)
        retained.to_csv(out / f"{name}_original_rows.tsv", sep="\t", index=False)
        if name == "mvp":
            number = lambda key: pd.to_numeric(retained[key], errors="coerce")
            se = number("standard_error")
            from_ci = se.isna()
            se_ci = (np.log(number("ci_upper")) - np.log(number("ci_lower"))) / (2 * 1.959963984540054)
            se = se.where(~from_ci, se_ci)
            n = pd.DataFrame({"rsid": retained.rsid, "EA": retained.effect_allele, "OA": retained.other_allele,
                              "beta": np.log(number("odds_ratio")), "se": se, "eaf": number("effect_allele_frequency"),
                              "p_reported": number("p_value"), "chr": retained.chromosome,
                              "pos": number("base_pair_location"), "build": "GRCh38",
                              "se_source": np.where(from_ci, "reconstructed_from_rounded_OR_95CI", "reported_standard_error")})
        elif name == "finngen_r12":
            number = lambda key: pd.to_numeric(retained[key], errors="coerce")
            n = pd.DataFrame({"rsid": retained.rsids, "EA": retained.alt, "OA": retained.ref,
                              "beta": number("beta"), "se": number("sebeta"), "eaf": number("af_alt"),
                              "p_reported": number("pval"), "chr": retained["#chrom"], "pos": number("pos"),
                              "build": "GRCh38", "se_source": "reported_standard_error"})
        else:
            number = lambda key: pd.to_numeric(retained[key], errors="coerce")
            n = pd.DataFrame({"rsid": retained.hm_rsid, "EA": retained.hm_effect_allele, "OA": retained.hm_other_allele,
                              "beta": number("hm_beta"), "se": number("standard_error"), "eaf": np.nan,
                              "p_reported": number("p_value"), "chr": retained.hm_chrom,
                              "pos": number("hm_pos"), "build": "GRCh38",
                              "se_source": "reported_standard_error_invariant_under_allele_flip"})
        n["rsid"] = n["rsid"].str.replace(";", ",", regex=False).str.split(",")
        n = n.explode("rsid")
        n = n[n.rsid.isin(rsids)]
        n["usable_effect"] = (np.isfinite(n.beta) & np.isfinite(n.se) & (n.se > 0))
        target = out / f"{name}_normalized.tsv"
        n.to_csv(target, sep="\t", index=False)
        inventories.append({"outcome": name, "source_rows": count, "retained_source_rows": len(retained),
                            "normalized_rows": len(n), "usable_effect_rows": int(n.usable_effect.sum()),
                            "standard_error_sources": n.se_source.value_counts().to_dict(),
                            "source": str(source.relative_to(ROOT)), "source_sha256": digest(source),
                            "script_sha256": digest(Path(__file__)),
                            "normalized_sha256": digest(target)})
        print(name, count, len(n), flush=True)
    record = {"input_exposure_sha256": digest(exposure), "requested_rsids": len(rsids),
              "script_sha256": digest(Path(__file__)),
              "join": "explicit rsID and allele matching; not cross-build coordinate joining",
              "outcomes": inventories}
    record_path.write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outcomes", nargs="+", choices=["mvp", "finngen_r12", "stuart"], default=["mvp", "finngen_r12"])
    main(parser.parse_args().outcomes)
