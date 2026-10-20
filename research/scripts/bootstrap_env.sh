#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."

export CONDA_PKGS_DIRS="$PWD/research/cache/conda"
export XDG_CACHE_HOME="$PWD/research/cache"
export TMPDIR="$PWD/research/cache/tmp"
mkdir -p "$CONDA_PKGS_DIRS" "$TMPDIR"

R_PREFIX="$PWD/research/envs/r442"
CONDA_LOCK="$PWD/research/records/conda_explicit.txt"
RENV_LOCK="$PWD/research/renv.lock"

validate_r_environment() {
  "$1" --vanilla - <<'RS'
if (!requireNamespace("jsonlite", quietly=TRUE)) stop("jsonlite is required to validate renv.lock")
lock <- jsonlite::fromJSON("research/renv.lock", simplifyVector=FALSE)
expected <- vapply(lock$Packages, function(x) x$Version, character(1))
actual <- installed.packages()
missing <- setdiff(names(expected), rownames(actual))
present <- intersect(names(expected), rownames(actual))
drift <- present[expected[present] != actual[present, "Version"]]
if (as.character(getRversion()) != lock$R$Version) {
  stop(sprintf("R version mismatch: expected %s, found %s", lock$R$Version, getRversion()))
}
if (length(missing)) stop(paste("Locked R packages missing:", paste(missing, collapse=", ")))
if (length(drift)) {
  details <- sprintf("%s expected %s found %s", drift, expected[drift], actual[drift, "Version"])
  stop(paste("Locked R package version mismatch:", paste(details, collapse="; ")))
}
RS
}

if [[ -x "$R_PREFIX/bin/Rscript" ]]; then
  validate_r_environment "$R_PREFIX/bin/Rscript"
else
  if [[ ! -f "$CONDA_LOCK" || ! -f "$RENV_LOCK" ]]; then
    echo "Required R lock file missing: $CONDA_LOCK or $RENV_LOCK" >&2
    exit 1
  fi
  command -v conda >/dev/null || { echo "conda is required to restore the locked R environment" >&2; exit 1; }
  export RENV_PATHS_ROOT="$PWD/research/cache/renv-root"
  export RENV_PATHS_CACHE="$PWD/research/cache/renv"
  conda create --yes --prefix "$R_PREFIX" --file "$CONDA_LOCK"
  "$R_PREFIX/bin/Rscript" --vanilla -e \
    'renv::restore(project="research", library=.Library, lockfile="research/renv.lock", prompt=FALSE)'
  validate_r_environment "$R_PREFIX/bin/Rscript"
fi
