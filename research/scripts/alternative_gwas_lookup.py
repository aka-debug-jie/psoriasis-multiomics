"""Public summary-GWAS checks in place of the omitted UKB individual-data branch.

These are per-variant associations, not PGS, survival estimates or a trans-ancestry
protein MR. MVP GRCh38 coordinates are explicitly converted with the UCSC
hg38ToHg19 chain before matching the GRCh37 alternative datasets.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import zipfile
import numpy as np
import pandas as pd
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def orient(ea, oa, other_ea, other_oa, beta, allow_complement=True):
    a = [str(x).upper() for x in [ea, oa, other_ea, other_oa]]
    if any(x not in "ACGT" or len(x) != 1 for x in a):
        return np.nan, "non_SNV"
    if a[0] == a[1] or a[2] == a[3]:
        return np.nan, "invalid_alleles"
    if set(a[:2]) in [set("AT"), set("CG")]:
        return np.nan, "palindromic_unresolved"
    comp = str.maketrans("ACGT", "TGCA")
    pair = a[2:]
    if a[:2] == pair:
        return beta, "aligned"
    if a[:2] == pair[::-1]:
        return -beta, "swapped"
    if not allow_complement:
        return np.nan, "allele_mismatch"
    pair = [x.translate(comp) for x in pair]
    if a[:2] == pair:
        return beta, "strand_aligned"
    if a[:2] == pair[::-1]:
        return -beta, "strand_swapped"
    return np.nan, "allele_mismatch"


def compare(rows, baseline, name, out):
    paired = rows.merge(baseline, on="rsid", suffixes=("_alt", "_mvp"))
    aligned = [orient(r.EA_mvp, r.OA_mvp, r.EA_alt, r.OA_alt, r.beta_alt, allow_complement=False)
               for r in paired.itertuples()]
    paired["beta_alt_aligned_to_mvp"] = [a[0] for a in aligned]
    paired["harmonization"] = [a[1] for a in aligned]
    usable = paired[np.isfinite(paired.beta_alt_aligned_to_mvp) & np.isfinite(paired.se_alt) & (paired.se_alt > 0)].copy()
    # Multiple aliases/alleles must not create duplicate independent observations.
    unique_effects = usable.drop_duplicates(["rsid", "beta_alt_aligned_to_mvp", "se_alt", "chr", "pos"])
    conflicts = unique_effects.rsid[unique_effects.rsid.duplicated(keep=False)].unique()
    paired["conflicting_outcome_match"] = paired.rsid.isin(conflicts)
    paired.to_csv(out / f"{name}_all_matches.tsv", sep="\t", index=False)
    usable = unique_effects[~unique_effects.rsid.isin(conflicts)].copy()
    usable["direction_agrees_with_mvp"] = np.sign(usable.beta_alt_aligned_to_mvp) == np.sign(usable.beta_mvp)
    usable["nonzero_direction"] = (usable.beta_alt_aligned_to_mvp != 0) & (usable.beta_mvp != 0)
    usable["evidence_type"] = "public_summary_variant_association_not_PGS_Cox"
    usable.to_csv(out / f"{name}_variant_comparison.tsv", sep="\t", index=False)
    directional = usable[usable.nonzero_direction]
    return {"matched_variants": len(usable), "ambiguous_variants_excluded": len(conflicts),
            "nonzero_direction_comparisons": len(directional),
            "same_direction": int(directional.direction_agrees_with_mvp.sum()),
            "alternative_nominal_p_below_0_05": int((usable.p_alt < .05).sum())}


def main():
    out = ROOT / "research/results/public_alternatives"
    out.mkdir(parents=True, exist_ok=True)
    status_file = out / "run.json"
    started = datetime.now(timezone.utc).isoformat()
    status_file.write_text(json.dumps({"status": "running", "started_utc": started}) + "\n")
    instruments_file = ROOT / "research/results/mr_available/cohort_wald.tsv"
    mvp_file = ROOT / "research/data/processed/outcomes/mvp_normalized.tsv"
    chain_file = ROOT / "research/data/metadata/hg38ToHg19.over.chain.gz"
    japan_file = ROOT / "research/data/raw/alternatives/hum0197.v25.PsV-gwas.v1.zip"
    psa_file = ROOT / "research/data/raw/alternatives/pmid35507331_psa_gwas_hg19.txt.gz"
    instruments = pd.read_csv(instruments_file, sep="\t")
    selected = set(instruments.rsid)
    original_mvp = pd.read_csv(mvp_file, sep="\t", dtype={"chr": str})
    original_mvp = original_mvp[original_mvp.rsid.isin(selected)].copy()
    baseline = original_mvp[["rsid", "EA", "OA", "beta", "se", "p_reported"]].drop_duplicates()
    baseline = baseline.rename(columns={"p_reported": "p"})
    if baseline.rsid.duplicated().any():
        raise ValueError("Nonunique baseline effects require resolution")
    lifter = LiftOver(str(chain_file))
    bridge_rows, excluded = [], []
    for r in original_mvp[["rsid", "chr", "pos"]].drop_duplicates().itertuples(index=False):
        chrom = {"23": "X", "24": "Y", "MT": "M"}.get(r.chr, r.chr)
        mapped = lifter.convert_coordinate("chr" + chrom.removeprefix("chr"), int(r.pos) - 1) or []
        if len(mapped) != 1:
            excluded.append({"rsid": r.rsid, "source_chr": r.chr, "source_pos": r.pos, "mapping_count": len(mapped)})
            continue
        new_chr, new_pos, strand, score = mapped[0]
        bridge_rows.append({"rsid": r.rsid, "chr": new_chr.removeprefix("chr"), "pos": new_pos + 1,
                            "liftover_strand": strand, "source_chr": r.chr, "source_pos": int(r.pos)})
    bridge = pd.DataFrame(bridge_rows)
    bridge["position_key"] = bridge["chr"] + ":" + bridge.pos.astype(str)
    pd.DataFrame(excluded, columns=["rsid", "source_chr", "source_pos", "mapping_count"]).to_csv(out / "liftover_excluded.tsv", sep="\t", index=False)
    bridge.to_csv(out / "GRCh37_variant_bridge.tsv", sep="\t", index=False)
    keys = set(bridge.position_key)
    japan_parts, japan_rows = [], 0
    with zipfile.ZipFile(japan_file) as archive:
        name = next(n for n in archive.namelist() if n.endswith(".txt"))
        with archive.open(name) as handle:
            for block in pd.read_csv(handle, sep="\t", dtype=str, chunksize=100000):
                japan_rows += len(block)
                key = block.MarkerName.str.split("_", n=1).str[0].str.removeprefix("chr")
                keep = key.isin(keys)
                if keep.any():
                    piece = block.loc[keep].copy(); piece["position_key"] = key[keep]
                    japan_parts.append(piece)
    japan = pd.concat(japan_parts, ignore_index=True).merge(bridge, on="position_key")
    japan.to_csv(out / "japan_original_rows.tsv", sep="\t", index=False)
    j = pd.DataFrame({"rsid": japan.rsid, "chr": japan["chr"], "pos": japan.pos,
                      "EA": japan.Allele2.str.upper(), "OA": japan.Allele1.str.upper(),
                      "beta": pd.to_numeric(japan.Effect), "se": pd.to_numeric(japan.StdErr),
                      "p": pd.to_numeric(japan["P-value"])})
    complement = str.maketrans("ACGT", "TGCA")
    negative = japan.liftover_strand == "-"
    for column in ["EA", "OA"]:
        j.loc[negative, column] = j.loc[negative, column].str.translate(complement)
    japan_summary = compare(j, baseline, "japan_psoriasis", out)
    print("Japan", japan_rows, japan_summary, flush=True)
    psa_parts, psa_rows = [], 0
    for block in pd.read_csv(psa_file, sep="\t", dtype=str, chunksize=100000):
        psa_rows += len(block)
        keep = block.rsid.isin(selected)
        if keep.any(): psa_parts.append(block.loc[keep].copy())
    psa = pd.concat(psa_parts, ignore_index=True)
    psa.to_csv(out / "psa_original_rows.tsv", sep="\t", index=False)
    psa = psa.merge(bridge, on="rsid")
    coordinate_ok = (psa.chromosome.str.removeprefix("chr") == psa["chr"]) & (pd.to_numeric(psa.position) == psa.pos)
    psa.loc[~coordinate_ok].to_csv(out / "psa_coordinate_mismatches.tsv", sep="\t", index=False)
    psa = psa.loc[coordinate_ok].reset_index(drop=True)
    p = pd.DataFrame({"rsid": psa.rsid, "chr": psa.chromosome, "pos": psa.position,
                      "EA": psa.effect_allele.str.upper(), "OA": psa.non_effect_allele.str.upper(),
                      "beta": pd.to_numeric(psa.beta), "se": pd.to_numeric(psa.se), "p": pd.to_numeric(psa.p_value)})
    negative = psa.liftover_strand == "-"
    for column in ["EA", "OA"]:
        p.loc[negative, column] = p.loc[negative, column].str.translate(complement)
    psa_summary = compare(p, baseline, "european_psa", out)
    print("PsA", psa_rows, psa_summary, flush=True)
    produced = list(out.glob("*.tsv"))
    input_files = [instruments_file, mvp_file, chain_file, japan_file, psa_file,
                   ROOT / "research/references/Dictionary_hum0197v25.docx", ROOT / "research/references/psa_figshare_metadata.json"]
    summary = {"status": "public_summary_lookup_completed", "started_utc": started,
               "finished_utc": datetime.now(timezone.utc).isoformat(), "selected_unique_variants": len(selected),
               "GRCh37_bridged_variants": bridge.rsid.nunique(),
               "coordinate_mapping": "UCSC hg38ToHg19; GWAS 1-based -> chain 0-based -> GWAS 1-based; unique mappings only; negative-strand alleles complemented",
               "japan_psoriasis": dict(japan_summary, source_rows=japan_rows, cases=1415, controls=3968,
                                        interpretation="cross_ancestry_variant_lookup; not European protein MR"),
               "european_psa": dict(psa_summary, source_rows=psa_rows, cases=3609, controls=9192,
                                    interpretation="PsA case-control lookup; overlap with European discovery/replication cohorts not ruled out"),
               "individual_PGS_or_Cox_performed": False,
               "replacement_for_original_PGS_score_slots": False,
               "input_hashes": {str(p.relative_to(ROOT)): digest(p) for p in input_files},
               "output_hashes": {str(p.relative_to(ROOT)): digest(p) for p in produced},
               "script_sha256": digest(Path(__file__))}
    status_file.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        path = ROOT / "research/results/public_alternatives/run.json"
        if path.exists():
            state = json.loads(path.read_text()); state.update(status="failed", error=str(exc))
            path.write_text(json.dumps(state, indent=2) + "\n")
        raise
