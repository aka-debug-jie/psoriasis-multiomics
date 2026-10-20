#!/usr/bin/env python3
"""Create unique-header tables from the byte-preserved release data package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from itertools import zip_longest
from pathlib import Path


RENAMES = {
    "mr_mvp_cohorts.tsv": {
        "neg_log10_p": ("exposure_neg_log10_p", "mr_neg_log10_p"),
        "F_statistic": ("F_exposure_source", "F_exposure_recomputed"),
    },
    "mr_finngen_cohorts.tsv": {
        "neg_log10_p": ("exposure_neg_log10_p", "mr_neg_log10_p"),
        "F_statistic": ("F_exposure_source", "F_exposure_recomputed"),
    },
    "mr_stuart_cohorts.tsv": {
        "neg_log10_p": ("exposure_neg_log10_p", "mr_neg_log10_p"),
        "F_statistic": ("F_exposure_source", "F_exposure_recomputed"),
    },
}
SENTINEL = object()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        try:
            return next(csv.reader(stream, delimiter="\t"))
        except StopIteration as exc:
            raise ValueError(f"TSV is empty and has no header: {path}") from exc


def transformed_header(source_name: str, header: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    duplicates = {name: count for name, count in Counter(header).items() if count > 1}
    expected = RENAMES.get(source_name, {})
    if duplicates != {name: 2 for name in expected}:
        raise ValueError(f"Unexpected duplicate header(s) in {source_name}: {duplicates}; expected {expected}")
    names_by_field = {name: iter(new_names) for name, new_names in expected.items()}
    new_header = []
    mapping = []
    for index, old_name in enumerate(header, start=1):
        new_name = next(names_by_field[old_name]) if old_name in names_by_field else old_name
        new_header.append(new_name)
        mapping.append({"original_name": old_name, "index": str(index), "new_name": new_name})
    if len(set(new_header)) != len(new_header):
        raise ValueError(f"Headers are still duplicated after approved renames in {source_name}: {new_header}")
    return new_header, mapping


def rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        next(reader)
        yield from reader


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2], help="project root (defaults to this checkout)")
    args = parser.parse_args()
    research = args.root.resolve() / "research"
    data_root = research / "release" / "data"
    data_manifest_path = data_root / "manifest.json"
    if not data_manifest_path.is_file():
        raise FileNotFoundError(f"Release data manifest is missing: {data_manifest_path}")
    data_manifest = json.loads(data_manifest_path.read_text(encoding="utf-8"))
    files = data_manifest.get("files", [])
    if len(files) != 30:
        raise ValueError(f"Expected 30 curated data tables in data/manifest.json, found {len(files)}")
    listed = {entry["output"] for entry in files}
    if set(RENAMES) - listed:
        raise ValueError(f"Required duplicate-header tables missing from data manifest: {sorted(set(RENAMES) - listed)}")

    output_root = research / "release" / "tables"
    own_manifest_path = output_root / "manifest.json"
    owned: set[str] = set()
    if own_manifest_path.exists():
        owned = set(json.loads(own_manifest_path.read_text(encoding="utf-8")).get("outputs", []))
    expected_outputs = listed | {"field_mapping.tsv", "manifest.json"}
    collisions = sorted(name for name in expected_outputs if (output_root / name).exists() and name not in owned)
    if collisions:
        raise FileExistsError("Refusing to overwrite outputs not owned by this manifest: " + ", ".join(collisions))

    output_root.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".release-tables-", dir=output_root.parent))
    table_records = []
    mapping_rows = []
    try:
        for entry in files:
            relative = entry["output"]
            source = data_root / relative
            if not source.is_file():
                raise FileNotFoundError(f"Manifest-listed release table is missing: {source}")
            header = read_header(source)
            output_header, field_maps = transformed_header(relative, header)
            destination = stage / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            if relative in RENAMES:
                with source.open("r", encoding="utf-8-sig", newline="") as src, destination.open("w", encoding="utf-8", newline="") as dst:
                    reader = csv.reader(src, delimiter="\t")
                    next(reader)
                    writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
                    writer.writerow(output_header)
                    for record_number, row in enumerate(reader, start=1):
                        if len(row) != len(header):
                            raise ValueError(f"Malformed row {record_number} in {relative}: expected {len(header)} cells, got {len(row)}")
                        writer.writerow(row)
            else:
                shutil.copyfile(source, destination)

            source_hash = sha256(source)
            output_hash = sha256(destination)
            if relative not in RENAMES and source_hash != output_hash:
                raise ValueError(f"Byte copy changed unexpectedly: {relative}")
            original_rows = rows(source)
            output_rows = rows(destination)
            row_count = 0
            for row_number, (old, new) in enumerate(zip_longest(original_rows, output_rows, fillvalue=SENTINEL), start=1):
                if old is SENTINEL or new is SENTINEL:
                    raise ValueError(f"Record count changed in {relative} at data row {row_number}")
                if old != new:
                    raise ValueError(f"Cell values changed in {relative} at data row {row_number}")
                if len(old) != len(header) or len(new) != len(output_header):
                    raise ValueError(f"Column count changed in {relative} at data row {row_number}")
                row_count += 1
            table_records.append({
                "output": relative,
                "source": f"release/data/{relative}",
                "source_sha256": source_hash,
                "output_sha256": output_hash,
                "rows": row_count,
                "columns": len(output_header),
            })
            for original_field, output_field in zip(field_maps, output_header):
                mapping_rows.append({
                    "table": relative,
                    **original_field,
                    "source_sha256": source_hash,
                    "output_sha256": output_hash,
                })

        with (stage / "field_mapping.tsv").open("w", encoding="utf-8", newline="") as stream:
            fields = ["table", "original_name", "index", "new_name", "source_sha256", "output_sha256"]
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(mapping_rows)
        output_manifest = {
            "format": "curated-release-tables-v1",
            "source_manifest": "release/data/manifest.json",
            "outputs": sorted(expected_outputs),
            "files": table_records,
            "mapping": {
                "file": "field_mapping.tsv",
                "sha256": sha256(stage / "field_mapping.tsv"),
                "field_records": len(mapping_rows),
            },
        }
        (stage / "manifest.json").write_text(json.dumps(output_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        output_root.mkdir(parents=True, exist_ok=True)
        for name in sorted(expected_outputs):
            destination = output_root / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(stage / name, destination)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
    print(f"Curated {len(table_records)} tables into {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
