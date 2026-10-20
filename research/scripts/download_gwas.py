"""Resume EBI GWAS downloads, verifying every requested final file and its source."""
from contextlib import ExitStack
from pathlib import Path
import csv
import hashlib
import re
import signal
import subprocess
from datetime import datetime, timezone

if __package__:
    from .acquisition_state import file_lock, read_records, update_record, check_identity, failure_record
else:
    from acquisition_state import file_lock, read_records, update_record, check_identity, failure_record

ROOT = Path(__file__).resolve().parents[2]


def hashes(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
            md5.update(block)
    return sha.hexdigest(), md5.hexdigest()


def verify(path, row, expected_md5):
    if not path.is_file():
        raise ValueError(f"Requested file is missing: {path}")
    size = path.stat().st_size
    if size != int(row["expected_bytes"]):
        raise ValueError(f"Size mismatch: {path}: {size} != {row['expected_bytes']}")
    sha, md5 = hashes(path)
    if md5 != expected_md5:
        raise ValueError(f"Official MD5 mismatch: {path}")
    return sha, md5


def acquire_all(rows):
    manifest = ROOT / "research/records/download_manifest.json"
    checksums, previous, verified, errors, lines = {}, {}, {}, {}, []
    for row in rows:
        identifier = row["id"]
        if identifier == "stuart_harmonised":
            metadata = (ROOT / "research/data/metadata/stuart_harmonised.yaml").read_text()
            expected = re.search(r"data_file_md5sum:\s*([0-9a-f]{32})", metadata).group(1)
        else:
            name = "mvp_md5.txt" if identifier == "mvp_gwas" else "stuart_md5.txt"
            pairs = [line.split() for line in (ROOT / "research/data/metadata" / name).read_text().splitlines()]
            expected = next(x[0] for x in pairs if x[-1].lstrip("*") == Path(row["path"]).name)
        checksums[identifier] = expected
        old = next((r for r in read_records(manifest) if r["id"] == identifier), None)
        previous[identifier] = old
        # An identity conflict is not a failed retrieval of the old resource.
        if old:
            check_identity(row, old)
        target = ROOT / row["path"]
        if target.exists():
            try:
                verified[identifier] = verify(target, row, expected)
                if old and old.get("sha256") and old["sha256"] != verified[identifier][0]:
                    raise ValueError(f"Previously recorded SHA256 mismatch: {target}")
            except ValueError as exc:
                verified.pop(identifier, None)
                errors[identifier] = str(exc)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        reservation = dict(row, status="downloading", started_utc=datetime.now(timezone.utc).isoformat())
        if old and old.get("retrieved_utc"):
            reservation["retrieved_utc"] = old["retrieved_utc"]
        update_record(manifest, reservation, source=row)
        lines.extend([row["url"], f" dir={target.parent}", f" out={target.name}.part", f" checksum=md5={expected}"])
    code = 0
    if lines:
        spec = ROOT / "research/records/aria2_gwas_inputs.txt"
        spec.write_text("\n".join(lines) + "\n")
        command = ["aria2c", "--input-file", str(spec), "--max-concurrent-downloads=2", "--max-connection-per-server=8",
                   "--split=8", "--min-split-size=1M", "--continue=true", "--file-allocation=none",
                   "--auto-file-renaming=false", "--max-tries=2", "--retry-wait=3", "--connect-timeout=30",
                   "--timeout=30", "--lowest-speed-limit=10K", "--console-log-level=warn", "--summary-interval=60", "--enable-color=false"]
        child = subprocess.Popen(command)
        def stop(signum, frame):
            child.send_signal(signum)
            child.wait()
            raise SystemExit(128 + signum)
        saved = {sig: signal.signal(sig, stop) for sig in (signal.SIGINT, signal.SIGTERM)}
        try:
            code = child.wait()
        finally:
            for sig, handler in saved.items():
                signal.signal(sig, handler)
    failures = 0
    for row in rows:
        identifier = row["id"]
        target = ROOT / row["path"]
        partial = Path(str(target) + ".part")
        old = previous[identifier]
        try:
            if identifier in errors:
                raise ValueError(errors[identifier])
            if identifier not in verified:
                if Path(str(partial) + ".aria2").exists():
                    raise ValueError(f"Download is incomplete: {partial}")
                verified[identifier] = verify(partial, row, checksums[identifier])
                partial.replace(target)
            sha, md5 = verified[identifier]
            now = datetime.now(timezone.utc).isoformat()
            result = dict(row, status="downloaded", bytes=target.stat().st_size, sha256=sha,
                          md5=md5, official_md5_verified=True, verified_utc=now)
            if old and old.get("retrieved_utc"):
                result["retrieved_utc"] = old["retrieved_utc"]
            # Existing files can have an unknown original acquisition date.
            elif any(line == f" out={target.name}.part" for line in lines):
                result["retrieved_utc"] = now
            update_record(manifest, result, source=row)
            print("VERIFIED", identifier, md5, flush=True)
        except (OSError, ValueError) as exc:
            failures += 1
            result = failure_record(row, old, str(exc), datetime.now(timezone.utc).isoformat())
            update_record(manifest, result, source=row)
            print("FAILED", identifier, str(exc), flush=True)
    return code or (1 if failures else 0)


def main():
    rows = list(csv.DictReader((ROOT / "research/config/public_sources.tsv").open(), delimiter="\t"))
    ids = {"mvp_gwas", "stuart_gwas", "stuart_harmonised"}
    rows = sorted([row for row in rows if row["id"] in ids], key=lambda row: row["path"])
    if {row["id"] for row in rows} != ids:
        raise ValueError("GWAS catalog must include all three requested source records")
    with ExitStack() as stack:
        for row in rows:
            stack.enter_context(file_lock(ROOT / row["path"]))
        code = acquire_all(rows)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
