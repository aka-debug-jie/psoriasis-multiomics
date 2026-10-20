"""Acquire original reference documents and extract author-reported tables.

These tables are comparison inputs, not a proteome-wide discovery dataset.
Run from any directory with the system Python (requests and beautifulsoup4).
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import io
import json
import subprocess
import zipfile
import xml.etree.ElementTree as ET

import requests
from bs4 import BeautifulSoup

if __package__:
    from .acquisition_state import file_lock, read_records, update_record, atomic_write, url_identity
else:
    from acquisition_state import file_lock, read_records, update_record, atomic_write, url_identity

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "research/references"
META = ROOT / "research/data/metadata"
RECORDS = ROOT / "research/records"
MANIFEST = RECORDS / "reference_manifest.json"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url, path, **kwargs):
    with file_lock(path):
        return fetch_locked(url, path, **kwargs)


def fetch_locked(url, path, **kwargs):
    records = read_records(MANIFEST)
    request_url = requests.Request("GET", url, params=kwargs.get("params")).prepare().url
    previous = next((x for x in records if x["path"] == str(path.relative_to(ROOT))), None)
    if previous and url_identity(request_url) != url_identity(previous.get("request_url", previous["url"])):
        raise ValueError(f"Reference source changed: {path}; preserve it under a separately versioned path")
    if path.exists() and previous:
        data = path.read_bytes()
        if sha256(data) != previous["sha256"]:
            raise ValueError(f"Previously acquired file changed: {path}")
        return data
    if path.exists():
        raise ValueError(f"Existing reference has no acquisition record: {path}")
    response = requests.get(url, timeout=60, **kwargs)
    response.raise_for_status()
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write(path, response.content)
    update_record(MANIFEST, {"path": str(path.relative_to(ROOT)), "url": response.url,
                    "request_url": request_url,
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                    "bytes": len(response.content), "sha256": sha256(response.content)}, key="path")
    return response.content


def write_tsv(path, rows):
    with path.open("w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(rows)


def main():
    for p in (REF, META, RECORDS):
        p.mkdir(parents=True, exist_ok=True)
    paper = ROOT / "btag347.pdf"
    (RECORDS / "paper_identity.json").write_text(json.dumps({
        "path": "btag347.pdf", "doi": "10.1093/bioinformatics/btag347",
        "sha256": sha256(paper.read_bytes()), "bytes": paper.stat().st_size
    }, indent=2) + "\n")
    subprocess.run(["pdftotext", "-layout", str(paper), str(REF / "paper.txt")], check=True)
    article = fetch("https://oup.silverchair-cdn.com/article-minimal/8697682", REF / "article.html")
    soup = BeautifulSoup(article, "html.parser")
    main_table = soup.find("table")
    reported_rows = [[c.get_text(" ", strip=True) for c in row.find_all(["td", "th"], recursive=False)]
                     for row in main_table.find_all("tr")[2:]]
    if len(reported_rows) != 29 or any(len(row) != 16 for row in reported_rows):
        raise ValueError("Unexpected Table 1 structure; inspect before extracting")
    header = ["target", "pwas", "sensitivity", "smr_blood", "smr_skin", "pgs_psoriasis", "pgs_psa",
              "bulk", "single_cell", "coloc", "ppi", "author_tier_score", "prior_gwas", "prior_observational",
              "prior_mr", "drug_phase", "evidence_origin"]
    write_tsv(REF / "author_main_table1.tsv", [header] + [row + ["author_reported"] for row in reported_rows])
    url = next(a["href"] for a in soup.find_all("a", href=True)
               if "supplementary_data.zip" in a["href"])
    archive = fetch(url, REF / "btag347_supplementary_data.zip")
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    table_index = []
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        for filename in ("Supplement Figures.docx", "Supplement Tables.docx"):
            data = z.read(filename)
            (REF / filename).write_bytes(data)
            with zipfile.ZipFile(io.BytesIO(data)) as doc:
                body = ET.fromstring(doc.read("word/document.xml")).find("w:body", ns)
            paragraphs, count, caption = [], 0, ""
            for node in body:
                if node.tag.endswith("}p"):
                    text = "".join(t.text or "" for t in node.findall(".//w:t", ns))
                    if text.strip():
                        paragraphs.append(text)
                        if text.startswith("Table S"):
                            caption = text
                elif node.tag.endswith("}tbl"):
                    count += 1
                    rows = [["".join(t.text or "" for t in cell.findall(".//w:t", ns))
                             for cell in row.findall("w:tc", ns)]
                            for row in node.findall("w:tr", ns)]
                    target = REF / f"author_table_S{count}.tsv"
                    write_tsv(target, rows)
                    table_index.append({"table": f"S{count}", "caption": caption,
                                        "rows_including_header": len(rows),
                                        "path": str(target.relative_to(ROOT)),
                                        "sha256": sha256(target.read_bytes()),
                                        "evidence_origin": "author_reported"})
            (REF / (Path(filename).stem + ".txt")).write_text("\n\n".join(paragraphs) + "\n")
    (REF / "table_index.json").write_text(json.dumps(table_index, indent=2) + "\n")
    fetch("https://zenodo.org/api/records/19692128", REF / "zenodo_19692128.json")
    for accession in ("GSE54456", "GSE173706"):
        fetch("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi", META / f"{accession}_family.soft",
              params={"acc": accession, "targ": "all", "form": "text", "view": "full"})
    print(json.dumps({"reference_files": len(read_records(MANIFEST)), "tables": table_index,
                      "zenodo_access": json.loads((REF / "zenodo_19692128.json").read_text())["metadata"]["access_right"]}, indent=2))


if __name__ == "__main__":
    main()
