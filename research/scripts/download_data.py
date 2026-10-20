"""Download only explicitly listed public inputs; retain hashes and original URLs."""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import argparse
import csv
import hashlib
import subprocess
import requests

if __package__:
    from .acquisition_state import file_lock, read_records, update_record, check_identity, failure_record
else:
    from acquisition_state import file_lock, read_records, update_record, check_identity, failure_record

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "research/records/download_manifest.json"


def file_hash(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_content(path, filename):
    with path.open("rb") as f:
        prefix = f.read(32768)
    if b"checking your browser" in prefix.lower() or b"preparing to download" in prefix.lower():
        raise ValueError("Server returned an interstitial page, not the requested resource")
    if filename.endswith(".gz") and prefix[:2] != b"\x1f\x8b":
        raise ValueError("Expected gzip input, received another format")
    if filename.endswith((".zip", ".xlsx", ".docx")) and prefix[:2] != b"PK":
        raise ValueError("Expected ZIP-based input, received another format")


def acquire(row, records=None, transport="requests"):
    with file_lock(ROOT / row["path"]):
        return acquire_locked(row, transport)


def acquire_locked(row, transport):
    target = ROOT / row["path"]
    old = next((x for x in read_records(MANIFEST) if x["id"] == row["id"]), None)
    if old:
        check_identity(row, old)
    if target.exists():
        if not old or old["status"] != "downloaded":
            raise ValueError(f"Existing input needs provenance/checksum verification before reuse: {target}")
        validate_content(target, target.name)
        if file_hash(target) != old["sha256"]:
            raise ValueError(f"Downloaded input changed: {target}")
        expected = row.get("expected_bytes", "")
        if expected and target.stat().st_size != int(expected):
            raise ValueError(f"Cached size mismatch: {target}")
        print(f"CACHED {row['id']}", flush=True)
        return old
    target.parent.mkdir(parents=True, exist_ok=True)
    temp = target.with_suffix(target.suffix + ".part")
    if not old and (temp.exists() or Path(str(temp) + ".aria2").exists()):
        raise ValueError(f"Untracked partial input cannot be resumed: {temp}")
    reservation = dict(row, status="downloading", started_utc=datetime.now(timezone.utc).isoformat())
    if old and old.get("retrieved_utc"):
        reservation["retrieved_utc"] = old["retrieved_utc"]
    update_record(MANIFEST, reservation, source=row)
    last_error = ""
    result = None
    for attempt in range(1 if transport == "aria2" else 2):
        try:
            if transport == "aria2":
                subprocess.run(["aria2c", "--file-allocation=none", "--continue=true", "--max-tries=2",
                                "--retry-wait=3", "--connect-timeout=30", "--timeout=90", "--max-connection-per-server=8",
                                "--split=8", "--min-split-size=8M", "--auto-file-renaming=false", "--summary-interval=60",
                                "--console-log-level=warn", "--dir", str(target.parent), "--out", temp.name, row["url"]], check=True)
            else:
                with requests.get(row["url"], stream=True, timeout=(20, 90)) as r:
                    r.raise_for_status()
                    with temp.open("wb") as f:
                        for chunk in r.iter_content(1024 * 1024):
                            f.write(chunk)
            expected = row.get("expected_bytes", "")
            validate_content(temp, target.name)
            if expected and temp.stat().st_size != int(expected):
                raise ValueError(f"Size mismatch: {temp.stat().st_size} != {expected}")
            temp.replace(target)
            result = dict(row, status="downloaded", bytes=target.stat().st_size,
                          sha256=file_hash(target), retrieved_utc=datetime.now(timezone.utc).isoformat())
            break
        except (requests.RequestException, OSError, ValueError, subprocess.CalledProcessError) as exc:
            last_error = str(exc)
            if isinstance(exc, requests.HTTPError) and exc.response.status_code in (401, 403):
                break
    if result is None:
        result = failure_record(row, old, last_error, datetime.now(timezone.utc).isoformat())
    update_record(MANIFEST, result, source=row)
    print(f"{result['status'].upper()} {row['id']} {result.get('bytes', result.get('error'))}", flush=True)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--group")
    p.add_argument("--ids", nargs="+")
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--transport", choices=["requests", "aria2"], default="requests")
    args = p.parse_args()
    rows = list(csv.DictReader((ROOT / "research/config/public_sources.tsv").open(), delimiter="\t"))
    if not args.group and not args.ids:
        p.error("Specify --group or --ids")
    if args.group:
        rows = [r for r in rows if r["group"] in args.group.split(",")]
    if args.ids:
        rows = [r for r in rows if r["id"] in args.ids]
    if not rows:
        raise ValueError("No resources for requested group")
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        results = list(executor.map(lambda row: acquire(row, transport=args.transport), rows))
    if any(r["status"] != "downloaded" for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
