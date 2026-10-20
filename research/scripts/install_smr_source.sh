#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SOURCE="$ROOT/research/tools/SMR-source"
BUILD="$ROOT/research/tools/smr-local"
RECORDS="$ROOT/research/records/advance_20261017"
REPOSITORY="https://github.com/JianYang-Lab/SMR.git"
COMMIT="eedb1a62ba4c66ab0b09acd0205876955a803f18"
CXX="${CXX:-g++}"
BUILD_JOBS="${BUILD_JOBS:-2}"
SOURCE_SNAPSHOT="$BUILD/source-tree"
BUILD_DIR="$BUILD/build"
INSTALL_DIR="$BUILD/install"

command -v git >/dev/null
command -v cmake >/dev/null
command -v "$CXX" >/dev/null
command -v tar >/dev/null
command -v python3 >/dev/null
mkdir -p "$RECORDS"

if [[ ! -e "$SOURCE" ]]; then
  mkdir -p "$SOURCE"
  git -C "$SOURCE" init
  git -C "$SOURCE" remote add origin "$REPOSITORY"
  git -C "$SOURCE" fetch --depth 1 origin "$COMMIT"
  git -C "$SOURCE" checkout --detach FETCH_HEAD
elif [[ ! -d "$SOURCE/.git" ]]; then
  echo "SMR source path exists without a git checkout; refusing to overwrite it" >&2
  exit 1
fi
if [[ "$(git -C "$SOURCE" rev-parse HEAD)" != "$COMMIT" ]]; then
  echo "SMR source checkout is not pinned commit $COMMIT" >&2
  exit 1
fi
if [[ -n "$(git -C "$SOURCE" status --porcelain)" ]]; then
  echo "SMR source checkout is not clean; refusing to build modified source" >&2
  exit 1
fi

mkdir -p "$SOURCE_SNAPSHOT"
git -C "$SOURCE" archive "$COMMIT" | tar -x -C "$SOURCE_SNAPSHOT"
export CPM_SOURCE_CACHE="$BUILD/cpm-source-cache"
export CXX
cmake -S "$SOURCE_SNAPSHOT" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
  -DSMR_BLAS_BACKEND=none \
  -DBUILD_WITH_MIMALLOC=OFF \
  -DCPM_SOURCE_CACHE="$CPM_SOURCE_CACHE" \
  > "$RECORDS/smr_build_configure.log" 2>&1
cmake --build "$BUILD_DIR" --parallel "$BUILD_JOBS" \
  > "$RECORDS/smr_build.log" 2>&1
cmake --install "$BUILD_DIR" \
  > "$RECORDS/smr_build_install.log" 2>&1

BINARY="$INSTALL_DIR/bin/smr"
[[ -x "$BINARY" ]] || { echo "SMR binary was not installed at $BINARY" >&2; exit 1; }
set +e
"$BINARY" > "$RECORDS/smr_build_help.txt" 2>&1
HELP_STATUS=$?
set -e
grep -q 'Version 1.4.3' "$RECORDS/smr_build_help.txt" || {
  echo "Could not read SMR version from no-argument help output" >&2
  exit 1
}
grep -q 'flags include:' "$RECORDS/smr_build_help.txt" || {
  echo "Could not read SMR option list from no-argument help output" >&2
  exit 1
}

SMR_ROOT="$ROOT" SMR_SOURCE="$SOURCE" SMR_BINARY="$BINARY" SMR_COMMIT="$COMMIT" \
SMR_CXX="$CXX" SMR_HELP_STATUS="$HELP_STATUS" python3 - <<'PY'
import hashlib
import json
import os
import subprocess
from pathlib import Path

root = Path(os.environ["SMR_ROOT"])
source = Path(os.environ["SMR_SOURCE"])
binary = Path(os.environ["SMR_BINARY"])
records = root / "research/records/advance_20261017"
cxx = os.environ["SMR_CXX"]
compiler_output = subprocess.run([cxx, "--version"], check=True, text=True, capture_output=True).stdout
cmake_output = subprocess.run(["cmake", "--version"], check=True, text=True, capture_output=True).stdout
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
help_path = records / "smr_build_help.txt"
text = help_path.read_text(encoding="utf-8", errors="replace")
version_line = next(line.strip("* ") for line in text.splitlines() if "Version " in line)
record = {
    "repository": "https://github.com/JianYang-Lab/SMR.git",
    "commit": os.environ["SMR_COMMIT"],
    "source_worktree_commit": subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"], check=True, text=True, capture_output=True
    ).stdout.strip(),
    "version_output": version_line,
    "binary_path": str(binary.relative_to(root)),
    "binary_sha256": sha(binary),
    "compiler_path": subprocess.run(["which", cxx], check=True, text=True, capture_output=True).stdout.strip(),
    "compiler_version": compiler_output.splitlines()[0],
    "cmake_version": cmake_output.splitlines()[0],
    "build_type": "Release",
    "blas_backend": "none",
    "mimalloc": False,
    "help_command": "smr (no arguments; upstream prints the option list and exits nonzero)",
    "help_exit_status": int(os.environ["SMR_HELP_STATUS"]),
    "help_output": str(help_path.relative_to(root)),
    "help_sha256": sha(help_path),
    "scope": "Local tool build and version/help inspection only; no SMR/HEIDI analysis was run.",
}
(records / "smr_build.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
PY

printf 'Built %s\n' "$BINARY"
printf 'Source commit %s\n' "$COMMIT"
printf 'Binary SHA256 '
sha256sum "$BINARY" | cut -d ' ' -f 1
