#!/usr/bin/env python3
"""Download official ChEMBL mechanisms for phase-4 psoriasis indications."""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import requests

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "research/data/raw/functional/chembl_psoriasis_indications.json"
MECHANISMS = ROOT / "research/data/raw/functional/chembl_mechanisms"
TARGETS = ROOT / "research/data/raw/functional/chembl_targets"
MECHANISM_URL = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
TARGET_URL = "https://www.ebi.ac.uk/chembl/api/data/target/{}.json"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def atomic_bytes(path: Path, content: bytes):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_bytes(content)
    tmp.replace(path)


def phase4_parents(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("drug_indications"), list):
        raise ValueError("input lacks drug_indications list")
    parents = set()
    for row in payload["drug_indications"]:
        try:
            if float(row["max_phase_for_ind"]) == 4.0:
                parent = row["parent_molecule_chembl_id"]
                if not parent:
                    raise ValueError("empty parent_molecule_chembl_id on phase-4 row")
                parents.add(parent)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid phase field/parent in indication row: {row}") from exc
    return sorted(parents)


def parse_json_response(content: bytes, expected_key: str):
    try:
        payload = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"response is not valid JSON: {exc}") from exc
    if not isinstance(payload, dict) or expected_key not in payload:
        raise ValueError(f"response JSON lacks expected key {expected_key!r}")
    return payload


def request_bytes(url: str, params=None):
    last_error = None
    for attempt in range(2):
        try:
            response = requests.get(url, params=params, timeout=60)
            response.raise_for_status()
            return response.url, response.headers.get("Date"), response.content
        except requests.RequestException as exc:
            last_error = exc
            if attempt == 0:
                time.sleep(1)
    raise RuntimeError(f"request failed after one retry: {last_error}")


def valid_mechanism_cache(folder: Path, molecule: str):
    manifest_path = folder / "manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("status") != "complete" or manifest.get("molecule_chembl_id") != molecule:
            return None
        all_rows = []
        for i, item in enumerate(manifest["responses"]):
            path = folder / item["file"]
            content = path.read_bytes()
            if digest(content) != item["sha256"]:
                return None
            payload = parse_json_response(content, "mechanisms")
            if not isinstance(payload["mechanisms"], list) or not isinstance(payload.get("page_meta"), dict):
                return None
            if i < len(manifest["responses"]) - 1:
                if payload["page_meta"].get("next") != manifest["responses"][i + 1]["requested_url"]:
                    return None
            elif payload["page_meta"].get("next"):
                return None
            all_rows.extend(payload["mechanisms"])
        if len(all_rows) != manifest.get("mechanism_count"):
            return None
        return manifest, all_rows
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def fetch_mechanisms(molecule: str, limit: int):
    folder = MECHANISMS / molecule
    cached = valid_mechanism_cache(folder, molecule)
    if cached:
        manifest, mechanisms = cached
        return {"molecule": molecule, "manifest": manifest, "mechanisms": mechanisms, "cache": True}
    folder.mkdir(parents=True, exist_ok=True)
    params = {"molecule_chembl_id": molecule, "limit": limit}
    request_url = requests.Request("GET", MECHANISM_URL, params=params).prepare().url
    visited = set()
    responses = []
    mechanisms = []
    while request_url:
        if request_url in visited:
            raise RuntimeError(f"pagination loop for {molecule}: {request_url}")
        visited.add(request_url)
        actual_url, date_header, content = request_bytes(request_url)
        payload = parse_json_response(content, "mechanisms")
        if not isinstance(payload["mechanisms"], list) or not isinstance(payload.get("page_meta"), dict):
            raise ValueError(f"invalid mechanisms/page_meta structure for {molecule}: {actual_url}")
        page = len(responses)
        file_name = f"page_{page:04d}.json"
        atomic_bytes(folder / file_name, content)
        responses.append({"file": file_name, "requested_url": request_url, "response_url": actual_url,
                          "date": date_header, "sha256": digest(content),
                          "row_count": len(payload["mechanisms"])})
        mechanisms.extend(payload["mechanisms"])
        next_url = payload["page_meta"].get("next")
        request_url = urljoin(actual_url, next_url) if next_url else None
    if any(not isinstance(row, dict) or not row.get("molecule_chembl_id") for row in mechanisms):
        raise ValueError(f"mechanism rows missing molecule_chembl_id for {molecule}")
    manifest = {"status": "complete", "molecule_chembl_id": molecule,
                "source_url": MECHANISM_URL, "limit": limit, "responses": responses,
                "mechanism_count": len(mechanisms), "finished_utc": datetime.now(timezone.utc).isoformat()}
    atomic_json(folder / "manifest.json", manifest)
    return {"molecule": molecule, "manifest": manifest, "mechanisms": mechanisms, "cache": False}


def valid_target_cache(folder: Path, target_id: str):
    manifest_path, body_path = folder / "manifest.json", folder / "target.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        content = body_path.read_bytes()
        if manifest.get("status") != "complete" or manifest.get("target_chembl_id") != target_id:
            return None
        if digest(content) != manifest.get("sha256"):
            return None
        payload = json.loads(content)
        if not isinstance(payload, dict) or payload.get("target_chembl_id") != target_id:
            return None
        return manifest, payload
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def fetch_target(target_id: str):
    folder = TARGETS / target_id
    cached = valid_target_cache(folder, target_id)
    if cached:
        return {"target": target_id, "manifest": cached[0], "cache": True}
    folder.mkdir(parents=True, exist_ok=True)
    url = TARGET_URL.format(target_id)
    actual_url, date_header, content = request_bytes(url)
    try:
        payload = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"target response is not valid JSON for {target_id}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("target_chembl_id") != target_id:
        raise ValueError(f"target response id mismatch for {target_id}: {actual_url}")
    body = folder / "target.json"
    atomic_bytes(body, content)
    manifest = {"status": "complete", "target_chembl_id": target_id, "requested_url": url,
                "response_url": actual_url, "date": date_header, "sha256": digest(content),
                "finished_utc": datetime.now(timezone.utc).isoformat()}
    atomic_json(folder / "manifest.json", manifest)
    return {"target": target_id, "manifest": manifest, "cache": False}


def run(input_path=INPUT, mechanism_dir=MECHANISMS, target_dir=TARGETS, workers=3, limit=1000):
    global MECHANISMS, TARGETS
    MECHANISMS, TARGETS = mechanism_dir, target_dir
    parents = phase4_parents(input_path)
    if len(parents) != 45:
        raise ValueError(f"expected 45 unique phase-4 parent molecules, found {len(parents)}")
    mechanism_results, errors = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(3, workers)) as pool:
        futures = {pool.submit(fetch_mechanisms, molecule, limit): molecule for molecule in parents}
        for future in concurrent.futures.as_completed(futures):
            molecule = futures[future]
            try:
                mechanism_results.append(future.result())
            except Exception as exc:
                errors.append({"stage": "mechanism", "id": molecule, "error": repr(exc)})
    target_ids = sorted({row["target_chembl_id"] for result in mechanism_results
                         for row in result["mechanisms"] if row.get("target_chembl_id")})
    unresolved = [{"molecule_chembl_id": result["molecule"], "mechanism": row}
                  for result in mechanism_results for row in result["mechanisms"]
                  if not row.get("target_chembl_id")]
    unresolved_path = mechanism_dir / "unresolved_target_mechanisms.json"
    mechanism_dir.mkdir(parents=True, exist_ok=True)
    atomic_json(unresolved_path, unresolved)
    target_results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(3, workers)) as pool:
        futures = {pool.submit(fetch_target, target_id): target_id for target_id in target_ids}
        for future in concurrent.futures.as_completed(futures):
            target_id = futures[future]
            try:
                target_results.append(future.result())
            except Exception as exc:
                errors.append({"stage": "target", "id": target_id, "error": repr(exc)})
    record = {
        "status": "complete" if not errors and len(mechanism_results) == len(parents) and len(target_results) == len(target_ids) else "incomplete",
        "finished_utc": datetime.now(timezone.utc).isoformat(), "input": str(input_path),
        "input_sha256": file_digest(input_path), "selection": "float(max_phase_for_ind)==4.0; unique parent_molecule_chembl_id",
        "phase4_parent_count": len(parents), "mechanism_parent_success_count": len(mechanism_results),
        "mechanism_rows": sum(len(r["mechanisms"]) for r in mechanism_results),
        "unresolved_target_mechanism_count": len(unresolved),
        "unresolved_target_mechanisms_file": unresolved_path.name,
        "mechanism_pages": sum(len(r["manifest"]["responses"]) for r in mechanism_results),
        "mechanism_cache_hits": sum(bool(r["cache"]) for r in mechanism_results),
        "unique_target_ids": len(target_ids), "target_success_count": len(target_results),
        "target_cache_hits": sum(bool(r["cache"]) for r in target_results),
        "errors": errors,
        "mechanism_manifests": {r["molecule"]: r["manifest"] for r in mechanism_results},
        "target_manifests": {r["target"]: r["manifest"] for r in target_results},
        "script_sha256": file_digest(Path(__file__)),
        "unresolved_target_mechanisms_sha256": file_digest(unresolved_path),
    }
    mechanism_dir.mkdir(parents=True, exist_ok=True); target_dir.mkdir(parents=True, exist_ok=True)
    atomic_json(mechanism_dir / "run.json", record)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--mechanism-dir", type=Path, default=MECHANISMS)
    parser.add_argument("--target-dir", type=Path, default=TARGETS)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--limit", type=int, default=1000)
    args = parser.parse_args()
    if args.workers < 1 or args.limit < 1:
        parser.error("--workers and --limit must be positive")
    result = run(args.input, args.mechanism_dir, args.target_dir, args.workers, args.limit)
    print(json.dumps({k: v for k, v in result.items() if k not in ("mechanism_manifests", "target_manifests")}, indent=2))
    if result["status"] != "complete":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
