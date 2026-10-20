"""Losslessly extract the two original Eldjarn pQTL worksheets; no IV selection."""
from pathlib import Path
import csv
import hashlib
import json
import openpyxl

ROOT = Path(__file__).resolve().parents[2]


def main():
    source = ROOT / "research/references/eldjarn_2023_supplement.xlsx"
    out = ROOT / "research/data/processed/exposures"
    out.mkdir(parents=True, exist_ok=True)
    book = openpyxl.load_workbook(source, read_only=True, data_only=True)
    inventory = []
    for sheet in ["ST19_somascan_pQTLs_smpnorm", "ST21_olink_ukb_bi_pQTLs"]:
        target = out / f"{sheet}.tsv"
        iterator = book[sheet].iter_rows(values_only=True)
        notes = [next(iterator) for _ in range(3)]
        headers = [str(x).replace("\n", " ") if x is not None else "" for x in next(iterator)]
        count = 0
        with target.open("w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(headers)
            for row in iterator:
                if all(x is None for x in row):
                    continue
                writer.writerow(row)
                count += 1
        inventory.append({"sheet": sheet, "rows": count, "columns": headers,
                          "notes": [list(r) for r in notes], "output": str(target.relative_to(ROOT)),
                          "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                          "scope": "source table extraction; no variant selection or SE reconstruction"})
        print(sheet, count, flush=True)
    record = {"source": str(source.relative_to(ROOT)),
              "sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "worksheets": inventory,
              "openpyxl_version": openpyxl.__version__}
    (ROOT / "research/records/exposure_table_inventory.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
