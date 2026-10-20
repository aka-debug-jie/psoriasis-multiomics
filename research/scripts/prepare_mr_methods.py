#!/usr/bin/env python3
"""Fetch pinned official TwoSampleMR method sources and run a small source smoke check."""

import hashlib
import json
import subprocess
import tempfile
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPOSITORY = "https://github.com/MRCIEU/TwoSampleMR.git"
TAG = "v0.6.22"
COMMIT = "d8b07182da23fb984157d9335fefcb83e894a944"
METHODS_DIR = ROOT / "research/tools/TwoSampleMR-methods"
RECORDS_DIR = ROOT / "research/records/advance_20261017"
FILES = ("R/mr.R", "R/mr_mode.R", "LICENSE")
REQUIRED_FUNCTIONS = ("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def verify_tag():
    result = subprocess.run(
        ["git", "ls-remote", REPOSITORY, f"refs/tags/{TAG}"],
        check=True, capture_output=True, text=True,
    )
    entries = [line.split() for line in result.stdout.splitlines() if line.strip()]
    tag_hashes = [entry[0] for entry in entries if len(entry) >= 2 and entry[1] == f"refs/tags/{TAG}"]
    if tag_hashes != [COMMIT]:
        raise RuntimeError(f"Official tag {TAG} does not resolve to pinned commit {COMMIT}: {entries}")
    return tag_hashes[0]


def fetch_sources():
    METHODS_DIR.mkdir(parents=True, exist_ok=True)
    source_records = []
    for relative_path in FILES:
        url = f"https://raw.githubusercontent.com/MRCIEU/TwoSampleMR/{COMMIT}/{relative_path}"
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read()
        target = METHODS_DIR / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if target.read_bytes() != data:
                raise RuntimeError(f"Refusing to overwrite modified or mismatched source file: {target}")
        else:
            target.write_bytes(data)
        source_records.append({
            "path": str(target.relative_to(ROOT)),
            "source_url": url,
            "commit": COMMIT,
            "sha256": sha256(data),
            "bytes": len(data),
        })
    return source_records


def run_r_smoke(records_dir):
    rscript = ROOT / "research/envs/r442_sc/bin/Rscript"
    if not rscript.is_file():
        raise FileNotFoundError(f"Project Rscript not found: {rscript}")
    output_path = records_dir / "twosamplemr_validation.json"
    r_script = r'''args <- commandArgs(trailingOnly = TRUE)
source(args[[1]], local = .GlobalEnv)
source(args[[2]], local = .GlobalEnv)
required <- c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
missing <- required[!vapply(required, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)]
if (length(missing)) stop("Missing functions: ", paste(missing, collapse = ", "))
parameters <- default_parameters()
signature <- function(name) {
  f <- get(name, envir = .GlobalEnv)
  fml <- formals(f)
  arg_names <- names(fml)
  args <- vapply(seq_along(arg_names), function(i) {
    value <- paste(deparse(fml[[i]]), collapse = "")
    if (!nzchar(trimws(value))) arg_names[[i]] else paste0(arg_names[[i]], " = ", value)
  }, character(1))
  paste0(name, "(", paste(args, collapse = ", "), ")")
}
set.seed(20261017)
b_exp <- c(0.08, 0.10, 0.12, 0.09, 0.11)
b_out <- c(0.040, 0.055, 0.050, 0.045, 0.060)
se_exp <- c(0.015, 0.018, 0.020, 0.016, 0.017)
se_out <- c(0.020, 0.021, 0.019, 0.018, 0.020)
results <- list(
  mr_ivw = mr_ivw(b_exp, b_out, se_exp, se_out, parameters),
  mr_weighted_median = mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters),
  mr_weighted_mode = mr_weighted_mode(b_exp, b_out, se_exp, se_out, parameters),
  mr_egger_regression = mr_egger_regression(b_exp, b_out, se_exp, se_out, parameters)
)
scalar_result <- lapply(results, function(x) {
  fields <- intersect(c("b", "se", "pval", "nsnp", "b_i", "se_i", "pval_i"), names(x))
  lapply(x[fields], function(value) as.numeric(value[[1]]))
})
for (name in names(scalar_result)) {
  vals <- unlist(scalar_result[[name]][intersect(c("b", "se", "pval"), names(scalar_result[[name]]))])
  if (length(vals) < 3 || any(!is.finite(vals))) stop("Non-finite smoke result for ", name)
}
record <- list(
  source_validation = "source_success",
  required_functions = as.list(required),
  signatures = as.list(setNames(vapply(required, signature, character(1)), required)),
  default_parameters = parameters,
  smoke = list(seed = 20261017L, instruments = length(b_exp), input = list(b_exp=b_exp, b_out=b_out, se_exp=se_exp, se_out=se_out), results=scalar_result),
  external_packages_loaded = character(),
  analysis_performed = FALSE
)
jsonlite::write_json(record, args[[3]], pretty=TRUE, auto_unbox=TRUE, digits=NA)
cat(jsonlite::toJSON(record, pretty=TRUE, auto_unbox=TRUE, digits=NA), "\n")
'''
    with tempfile.TemporaryDirectory(prefix="twosamplemr-smoke-") as temporary_dir:
        helper = Path(temporary_dir) / "smoke.R"
        helper.write_text(r_script, encoding="utf-8")
        result = subprocess.run(
            [str(rscript), "--vanilla", str(helper),
             str(METHODS_DIR / "R/mr.R"), str(METHODS_DIR / "R/mr_mode.R"), str(output_path)],
            cwd=ROOT, capture_output=True, text=True,
        )
        if result.returncode:
            raise RuntimeError(f"R source/smoke validation failed:\n{result.stdout}\n{result.stderr}")
        (records_dir / "twosamplemr_smoke.log").write_text(result.stdout, encoding="utf-8")
    return json.loads(output_path.read_text(encoding="utf-8"))


def main():
    commit = verify_tag()
    source_records = fetch_sources()
    RECORDS_DIR.mkdir(parents=True, exist_ok=True)
    validation = run_r_smoke(RECORDS_DIR)
    manifest = {
        "repository": REPOSITORY,
        "tag": TAG,
        "commit": commit,
        "source_files": source_records,
        "r_validation": str((RECORDS_DIR / "twosamplemr_validation.json").relative_to(ROOT)),
        "source_scope": "Only upstream R/mr.R, R/mr_mode.R, and LICENSE were fetched; source was not edited.",
        "analysis_performed": False,
    }
    (RECORDS_DIR / "twosamplemr_sources.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "tag": TAG,
        "commit": commit,
        "source_sha256": {Path(item["path"]).name: item["sha256"] for item in source_records},
        "functions": validation["signatures"],
        "default_parameters": validation["default_parameters"],
        "smoke_methods": list(validation["smoke"]["results"]),
        "analysis_performed": False,
    }, indent=2))


if __name__ == "__main__":
    main()
