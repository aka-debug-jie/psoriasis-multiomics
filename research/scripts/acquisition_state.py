"""Project-local file locks and atomic updates for acquisition records."""
from contextlib import contextmanager
from pathlib import Path
from urllib.parse import urlsplit, parse_qsl
import fcntl
import json
import os
import tempfile


@contextmanager
def file_lock(path):
    """All download entry points lock the same destination before touching it."""
    lock_path = Path(str(path) + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def read_records(path):
    return json.loads(path.read_text()) if path.exists() else []


def update_record(path, record, key="id", source=None):
    # Reload under an interprocess lock: a stale worker cannot erase other IDs.
    with file_lock(path):
        records = read_records(path)
        previous = next((r for r in records if r[key] == record[key]), None)
        if previous and source:
            check_identity(source, previous)
        records = [r for r in records if r[key] != record[key]] + [record]
        atomic_write(path, (json.dumps(records, ensure_ascii=False, indent=2) + "\n").encode())
    return records


def atomic_write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def url_identity(url):
    parts = urlsplit(url)
    return (parts.scheme.lower(), parts.netloc.lower(), parts.path,
            tuple(sorted(parse_qsl(parts.query, keep_blank_values=True), key=lambda pair: pair[0])))


def failure_record(row, old, error, checked_utc):
    result = dict(row, status="failed", error=error, checked_utc=checked_utc)
    if old:
        if old.get("retrieved_utc"):
            result["retrieved_utc"] = old["retrieved_utc"]
        evidence = old.get("previous_verification")
        if old.get("status") == "downloaded":
            evidence = {k: old[k] for k in ("sha256", "md5", "bytes", "verified_utc", "official_md5_verified") if k in old}
        if evidence:
            result["previous_verification"] = evidence
    return result


def check_identity(row, old):
    """Reject changed catalog metadata before reusing or replacing a resource."""
    changes = []
    for key, value in row.items():
        previous = old.get(key, "")
        same = url_identity(value) == url_identity(previous) if key == "url" else str(value) == str(previous)
        if not same:
            changes.append(key)
    if changes:
        raise ValueError(f"Source identity changed ({', '.join(changes)}) for {row.get('id', row.get('path'))}; use a separately versioned input path and record")
