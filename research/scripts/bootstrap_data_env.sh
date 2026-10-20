#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."

ENV_PREFIX="$PWD/research/envs/data"
REQUIREMENTS="$PWD/research/requirements-data.txt"
PYTHON_VERSION_RECORD="$PWD/research/records/python_version.txt"

validate_python_env() {
  "$1" - "$REQUIREMENTS" "$PYTHON_VERSION_RECORD" <<'PY'
import importlib.metadata
import re
import sys
from pathlib import Path

requirements_path, version_record_path = map(Path, sys.argv[1:])
recorded = version_record_path.read_text(encoding="utf-8").strip()
match = re.fullmatch(r"Python (\d+\.\d+\.\d+)", recorded)
if not match or match.group(1) != "3.12.7":
    raise SystemExit(f"Expected recorded Python version 3.12.7, found {recorded!r}")
expected_python = tuple(map(int, match.group(1).split(".")))
if sys.version_info[:3] != expected_python:
    raise SystemExit(f"Python version mismatch: expected {recorded}, found {sys.version.split()[0]}")

errors = []
for raw_line in requirements_path.read_text(encoding="utf-8").splitlines():
    line = raw_line.strip()
    if not line or line.startswith("#"):
        continue
    pinned = re.fullmatch(r"([A-Za-z0-9_.-]+)==([^\s;]+)", line)
    if not pinned:
        errors.append(f"requirement is not an exact version pin: {line}")
        continue
    name, expected = pinned.groups()
    try:
        actual = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        errors.append(f"{name} missing (expected {expected})")
    else:
        if actual != expected:
            errors.append(f"{name} expected {expected}, found {actual}")
if errors:
    raise SystemExit("Python environment does not match lock:\n" + "\n".join(errors))
PY
}

if [[ -x "$ENV_PREFIX/bin/python" ]]; then
  validate_python_env "$ENV_PREFIX/bin/python"
else
  [[ -f "$REQUIREMENTS" && -f "$PYTHON_VERSION_RECORD" ]] || {
    echo "Required Python lock files missing: $REQUIREMENTS or $PYTHON_VERSION_RECORD" >&2
    exit 1
  }
  expected_version="$(sed -n 's/^Python //p' "$PYTHON_VERSION_RECORD")"
  [[ "$expected_version" == "3.12.7" ]] || {
    echo "Expected recorded Python version 3.12.7, found: $expected_version" >&2
    exit 1
  }
  PYTHON_BOOTSTRAP="${PYTHON_BOOTSTRAP:-python3.12}"
  command -v "$PYTHON_BOOTSTRAP" >/dev/null || {
    echo "Python $expected_version interpreter not found: $PYTHON_BOOTSTRAP" >&2
    exit 1
  }
  bootstrap_version="$("$PYTHON_BOOTSTRAP" --version 2>&1)"
  [[ "$bootstrap_version" == "Python $expected_version" ]] || {
    echo "Bootstrap interpreter mismatch: expected Python $expected_version, found $bootstrap_version" >&2
    exit 1
  }
  "$PYTHON_BOOTSTRAP" -m venv "$ENV_PREFIX"
  export PIP_CACHE_DIR="$PWD/research/cache/pip"
  "$ENV_PREFIX/bin/python" -m pip install -r "$REQUIREMENTS"
  validate_python_env "$ENV_PREFIX/bin/python"
fi
