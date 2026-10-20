#!/usr/bin/env python3
"""Create a metadata-only audit of the two local GEO SOFT files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


FIELDS = (
    "gsm",
    "title",
    "donor_id",
    "subject_status",
    "tissue_type",
    "supplementary_url",
    "selected",
    "selection_reason",
)


def sample_blocks(text: str) -> list[dict[str, list[str]]]:
    """Parse repeated !Sample_* fields from SOFT sample blocks."""
    blocks: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            if current is not None:
                blocks.append(current)
            current = {"gsm": [line.split("=", 1)[1].strip()]}
        elif current is not None and line.startswith("!Sample_") and " = " in line:
            key, value = line[1:].split(" = ", 1)
            current.setdefault(key, []).append(value.strip())
    if current is not None:
        blocks.append(current)
    return blocks


def _first(block: dict[str, list[str]], key: str) -> str:
    values = block.get(key, [])
    return values[0] if values else ""


def _characteristics(block: dict[str, list[str]]) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in block.get("Sample_characteristics_ch1", []):
        if ":" in item:
            key, value = item.split(":", 1)
            values[key.strip().casefold()] = value.strip()
    return values


def _record(block: dict[str, list[str]], *, sc: bool) -> dict[str, str | bool]:
    chars = _characteristics(block)
    donor = chars.get("subject id", "")
    donor = re.sub(r"^donor\s+", "", donor, flags=re.IGNORECASE).strip()
    tissue = chars.get("tissue type", "")
    status = chars.get("subject status", "")

    if sc:
        if tissue.casefold() == "normal skin" and status.casefold() == "healthy donor":
            selected, reason = True, "healthy donor with Normal skin tissue"
        elif tissue.casefold() == "psoriatic skin" and status.casefold() == "psoriasis patient":
            selected, reason = True, "psoriasis patient with Psoriatic skin tissue"
        elif tissue.casefold() == "peripheral normal skin" and status.casefold() == "psoriasis patient":
            selected, reason = False, "peripheral normal skin excluded by selection rule"
        else:
            selected, reason = False, "not selected: characteristics do not match an inclusion rule"
    else:
        selected, reason = True, "included in bulk metadata audit"

    urls = block.get("Sample_supplementary_file_1", [])
    url = next((value for value in urls if value and value.upper() != "NONE"), "")
    return {
        "gsm": _first(block, "gsm"),
        "title": _first(block, "Sample_title"),
        "donor_id": donor,
        "subject_status": status,
        "tissue_type": tissue,
        "supplementary_url": url,
        "selected": selected,
        "selection_reason": reason,
    }


def _write_tsv(path: Path, rows: list[dict[str, str | bool]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = dict(row)
            output["selected"] = "true" if row["selected"] else "false"
            writer.writerow(output)


def _processing_text(text: str) -> str:
    return "\n".join(
        line.split(" = ", 1)[1]
        for line in text.splitlines()
        if line.startswith("!Sample_data_processing = ")
    )


def _validate_processing(blocks: list[dict[str, list[str]]], *, unit: str) -> None:
    for block in blocks:
        evidence = _processing_text("\n".join(
            f"!{key} = {value}"
            for key, values in block.items()
            if key == "Sample_data_processing"
            for value in values
        ))
        normalized = evidence.casefold()
        if unit == "RPKM":
            valid = "rpkm" in normalized
        elif unit == "raw gene counts":
            valid = "raw gene counts" in normalized
        else:
            raise ValueError(f"Unsupported expected unit: {unit}")
        if not valid:
            gsm = _first(block, "gsm") or "<unknown GSM>"
            raise ValueError(f"{gsm}: expected {unit!r} in Sample_data_processing; found {evidence!r}")


def audit(metadata_dir: Path, output_dir: Path) -> dict[str, object]:
    bulk_path = metadata_dir / "GSE54456_family.soft"
    sc_path = metadata_dir / "GSE173706_family.soft"
    bulk_text = bulk_path.read_text(encoding="utf-8")
    sc_text = sc_path.read_text(encoding="utf-8")
    bulk_blocks = sample_blocks(bulk_text)
    sc_blocks = sample_blocks(sc_text)
    _validate_processing(bulk_blocks, unit="RPKM")
    _validate_processing(sc_blocks, unit="raw gene counts")
    bulk = [_record(block, sc=False) for block in bulk_blocks]
    sc = [_record(block, sc=True) for block in sc_blocks]
    selected = [row for row in sc if row["selected"]]

    bulk_tissues = Counter(str(row["tissue_type"]).casefold() for row in bulk)
    sc_tissues = Counter(str(row["tissue_type"]).casefold() for row in sc)
    bulk_expected = {"normal skin": 82, "lesional psoriatic skin": 92}
    if len(bulk) != 174 or dict(bulk_tissues) != bulk_expected:
        raise ValueError(f"GSE54456 sample labels/counts do not match expected audit: {len(bulk)} {dict(bulk_tissues)}")
    expected_sc = {"normal skin": 8, "peripheral normal skin": 11, "psoriatic skin": 14}
    if len(sc) != 33 or dict(sc_tissues) != expected_sc or len(selected) != 22:
        raise ValueError(f"GSE173706 sample labels/counts do not match expected audit: {len(sc)} {dict(sc_tissues)} selected={len(selected)}")
    if any(not row["donor_id"] for row in sc):
        raise ValueError("GSE173706 contains a sample without a characteristics-derived subject id")

    warnings = []
    mismatch = next((row for row in sc if row["gsm"] == "GSM5277178"), None)
    if mismatch and "3696" in str(mismatch["title"]) and "PN-30696" in str(mismatch["title"]):
        warnings.append(
            "GSM5277178 retains SOFT title donor 3696 and filename token PN-30696 as recorded; values conflict and were not corrected."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_dir / "samples_bulk.tsv", bulk)
    _write_tsv(output_dir / "samples_sc.tsv", sc)
    _write_tsv(output_dir / "sc_selected_samples.tsv", selected)

    summary: dict[str, object] = {
        "datasets": {
            "GSE54456": {
                "samples": len(bulk),
                "tissue_type_counts": {"normal skin": 82, "lesional psoriatic skin": 92},
                "supplementary_matrix_unit": "RPKM",
                "unit_evidence": "SOFT Sample_data_processing states normalized to reads per kilobase per million mapped reads (RPKM) and describes the supplementary matrix as RPKM data matrix.",
            },
            "GSE173706": {
                "samples": len(sc),
                "selected": len(selected),
                "selected_tissue_type_counts": {"Normal skin": 8, "Psoriatic skin": 14},
                "excluded_tissue_type_counts": {"Peripheral normal skin": 11},
                "supplementary_matrix_unit": "raw gene counts",
                "unit_evidence": "SOFT Sample_data_processing describes each supplementary CSV as raw gene counts for every gene and every cell.",
            },
        },
        "warnings": warnings,
        "input_sha256": {
            "GSE54456_family.soft": hashlib.sha256((metadata_dir / "GSE54456_family.soft").read_bytes()).hexdigest(),
            "GSE173706_family.soft": hashlib.sha256((metadata_dir / "GSE173706_family.soft").read_bytes()).hexdigest(),
        },
        "scope": "SOFT metadata parsing and supplementary matrix unit audit only; no statistical analysis performed.",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    project_root = Path(__file__).resolve().parents[2]
    parser.add_argument("--metadata-dir", type=Path, default=project_root / "research/data/metadata")
    parser.add_argument("--output-dir", type=Path, default=project_root / "research/results/geo_audit")
    args = parser.parse_args()
    summary = audit(args.metadata_dir, args.output_dir)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
