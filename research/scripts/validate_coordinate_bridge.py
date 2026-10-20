"""Check three coordinate conversions against independent Ensembl GRCh37 records."""
from pathlib import Path
import hashlib
import json
import pandas as pd
import requests
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]


def main():
    chain = ROOT / "research/data/metadata/hg38ToHg19.over.chain.gz"
    source = ROOT / "research/data/processed/outcomes/mvp_normalized.tsv"
    lo = LiftOver(str(chain))
    mvp = pd.read_csv(source, sep="\t", dtype={"chr": str})
    checks, inputs = [], [chain, source]
    for rsid in ["rs1133071", "rs17783344", "rs5498"]:
        rows = mvp[mvp.rsid == rsid]
        assert len(rows) == 1
        row = rows.iloc[0]
        mapping = lo.convert_coordinate("chr" + row["chr"], int(row.pos) - 1)
        assert len(mapping) == 1
        path = ROOT / "research/data/metadata" / f"{rsid}_GRCh37.json"
        url = f"https://grch37.rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
        if not path.exists():
            response = requests.get(url, timeout=30); response.raise_for_status()
            path.write_text(json.dumps(response.json(), indent=2) + "\n")
        inputs.append(path)
        reference = json.loads(path.read_text())
        coordinates = {(x["seq_region_name"], x["start"]) for x in reference["mappings"] if x["assembly_name"] == "GRCh37"}
        target = (mapping[0][0].removeprefix("chr"), mapping[0][1] + 1)
        assert target in coordinates
        checks.append({"rsid": rsid, "hg38_chr": row["chr"], "hg38_position": int(row.pos),
                       "hg19_chr": target[0], "hg19_position": target[1], "matches_Ensembl_GRCh37": True, "reference_url": url})
    result = {"checked": len(checks), "checks": checks, "scope": "three named variants; not full-genome independent validation",
              "input_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT / "research/records/coordinate_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print("PASS: three GRCh38-to-GRCh37 conversions match independent Ensembl records")


if __name__ == "__main__":
    main()
