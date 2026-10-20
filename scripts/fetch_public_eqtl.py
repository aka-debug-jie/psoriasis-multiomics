"""Fetch two fixed eQTL Catalogue regions for the public teaching package."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
import zlib


SOURCE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
    "QTS000015/QTD000356/QTD000356.all.tsv.gz"
)
GENES = {
    "ENSG00000164308": "ERAP2",
    "ENSG00000107798": "LIPA",
}


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def get(url: str, **kwargs):
    import requests

    os.environ.setdefault("CURL_CA_BUNDLE", requests.certs.where())
    for attempt in range(2):
        try:
            response = requests.get(url, timeout=60, **kwargs)
            response.raise_for_status()
            return response
        except requests.RequestException:
            if attempt:
                raise
            time.sleep(2)
    raise AssertionError("unreachable")


def metadata_for(gene: str, metadata_dir: Path | None) -> dict:
    if metadata_dir is not None:
        return json.loads((metadata_dir / f"{gene}.json").read_text())
    response = get(
        f"https://rest.ensembl.org/lookup/id/{gene}",
        params={"content-type": "application/json"},
    )
    return response.json()


def source_header() -> list[str]:
    response = get(SOURCE_URL, headers={"Range": "bytes=0-65535"})
    if response.status_code != 206:
        raise ValueError("eQTL source did not honor the header range request")
    text = zlib.decompressobj(31).decompress(response.content).decode("utf-8")
    first = text.splitlines()[0].lstrip("#")
    return first.split("\t")


def fetch_gene(gene: str, meta: dict, header: list[str], index: Path, output: Path) -> dict:
    import pysam

    if meta.get("id") != gene:
        raise ValueError(f"metadata ID does not match {gene}")
    chrom = str(meta["seq_region_name"])
    tss = int(meta["start"] if int(meta["strand"]) == 1 else meta["end"])
    gene_column = header.index("gene_id")
    output.parent.mkdir(parents=True, exist_ok=True)

    fd, temporary = tempfile.mkstemp(prefix=f".{gene}.", suffix=".tsv.tmp", dir=output.parent)
    count = 0
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
            stream.write("\t".join(header) + "\n")
            with pysam.TabixFile(SOURCE_URL, index=str(index)) as tabix:
                for line in tabix.fetch(chrom, max(0, tss - 1_000_001), tss + 1_000_000):
                    fields = line.split("\t")
                    if fields[gene_column] == gene:
                        stream.write(line + "\n")
                        count += 1
            stream.flush()
            os.fsync(stream.fileno())
        if count == 0:
            raise ValueError(f"no regional records found for {gene}")
        os.replace(temporary, output)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
    return {"gene_id": gene, "symbol": GENES[gene], "tss": tss, "chromosome": chrom, "records": count}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Directory for gene TSVs, metadata, and receipt")
    parser.add_argument("--metadata-dir", type=Path, help="Optional directory containing pinned Ensembl <gene_id>.json files")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="public-eqtl-") as temporary_dir:
        index = Path(temporary_dir) / "QTD000356.all.tsv.gz.tbi"
        index.write_bytes(get(SOURCE_URL + ".tbi").content)
        header = source_header()
        results = []
        for gene in GENES:
            meta = metadata_for(gene, args.metadata_dir)
            result = fetch_gene(gene, meta, header, index, args.output / f"{gene}.tsv")
            metadata_path = args.output / f"{gene}.json"
            metadata_path.write_text(json.dumps(meta, indent=2) + "\n")
            result["tsv_sha256"] = digest(args.output / f"{gene}.tsv")
            result["metadata_sha256"] = digest(metadata_path)
            results.append(result)
        receipt = {
            "study_id": "QTS000015",
            "dataset_id": "QTD000356",
            "source_url": SOURCE_URL,
            "index_url": SOURCE_URL + ".tbi",
            "index_sha256": hashlib.sha256(index.read_bytes()).hexdigest(),
            "genes": results,
        }
        (args.output / "fetch_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
