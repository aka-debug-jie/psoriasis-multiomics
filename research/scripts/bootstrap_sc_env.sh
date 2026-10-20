#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
export CONDA_PKGS_DIRS="$PWD/research/cache/conda"
export XDG_CACHE_HOME="$PWD/research/cache"
SC_PREFIX="$PWD/research/envs/r442_sc"
if [[ ! -x "$SC_PREFIX/bin/Rscript" ]]; then
  conda create --yes --prefix "$SC_PREFIX" --file research/records/conda_sc_explicit.txt
fi
"$SC_PREFIX/bin/Rscript" --vanilla - <<'RS'
stopifnot(as.character(getRversion()) == "4.4.2")
expected <- read.csv("research/records/sc_R_packages.csv", stringsAsFactors=FALSE)
actual <- installed.packages()
stopifnot(all(expected$Package %in% rownames(actual)),
          all(expected$Version == actual[expected$Package, "Version"]))
stopifnot(as.character(packageVersion("Seurat")) == "5.3.0", as.character(packageVersion("Matrix")) == "1.7.5")
RS
"$SC_PREFIX/bin/Rscript" --vanilla research/scripts/install_presto.R
"$SC_PREFIX/bin/Rscript" --vanilla research/scripts/install_coloc.R
