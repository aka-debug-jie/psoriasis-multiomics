#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
library_path <- normalizePath(file.path(root, "research/envs/r442_sc/lib/R/library"), mustWork = TRUE)
runtime_bin <- file.path(root, "research/envs/r442_sc/bin")
manifest_path <- file.path(root, "research/records/download_manifest.json")
archive_expected <- file.path(root, "research/tools/sources/presto-b5df6ee.tar.gz")
expected_package_version <- "1.1.0"
expected_commit <- "b5df6ee6097eb62522f2e557aa93ac21bab2d05f"

if (!file.exists(file.path(runtime_bin, "x86_64-conda-linux-gnu-c++"))) {
  stop("r442_sc C++ compiler wrapper is missing; no package changes made")
}
Sys.setenv(PATH = paste(runtime_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
Sys.setenv(TMPDIR = file.path(root, "research/cache/tmp"))
.libPaths(c(library_path, .Library))

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
entry <- Filter(function(x) identical(x$id, "presto_1_1_0_source"), manifest)
if (length(entry) != 1L) stop("Expected exactly one presto_1_1_0_source download manifest entry")
entry <- entry[[1]]
if (!identical(entry$status, "downloaded") || !identical(entry$path, "research/tools/sources/presto-b5df6ee.tar.gz")) {
  stop("Presto archive path/status differs from the existing download manifest")
}
if (!grepl(expected_commit, entry$url, fixed = TRUE)) stop("Download manifest URL does not pin the expected Presto commit")
if (!file.exists(archive_expected)) stop("Pinned Presto source archive is missing")
archive_sha256 <- digest::digest(file = archive_expected, algo = "sha256", serialize = FALSE)
if (!identical(archive_sha256, entry$sha256)) {
  stop(sprintf("Pinned source SHA256 mismatch: manifest=%s file=%s", entry$sha256, archive_sha256))
}

description_member <- grep("(^|/)DESCRIPTION$", utils::untar(archive_expected, list = TRUE), value = TRUE)
if (length(description_member) != 1L) stop("Pinned source archive must contain exactly one DESCRIPTION")
description_dir <- tempfile("presto-description-")
dir.create(description_dir)
on.exit(unlink(description_dir, recursive = TRUE), add = TRUE)
utils::untar(archive_expected, files = description_member, exdir = description_dir)
description <- read.dcf(file.path(description_dir, description_member))
unlink(description_dir, recursive = TRUE)
if (description[1, "Package"] != "presto" || description[1, "Version"] != expected_package_version) {
  stop("Pinned source DESCRIPTION is not presto 1.1.0")
}

installed <- installed.packages(lib.loc = c(library_path, .Library))
installed_keep <- !duplicated(rownames(installed))
installed_versions <- setNames(installed[installed_keep, "Version"], rownames(installed)[installed_keep])

parse_dependencies <- function(desc) {
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(desc))
  text <- paste(desc[1, fields], collapse = ",")
  parts <- trimws(unlist(strsplit(gsub("\n", " ", text), ",", fixed = TRUE)))
  dependencies <- list()
  for (part in parts) {
    hit <- regmatches(part, regexec("^([A-Za-z][A-Za-z0-9.]*)\\s*(?:\\((>=|<=|==|>|<)\\s*([^)]*)\\))?$", part))[[1]]
    if (length(hit)) dependencies[[length(dependencies) + 1L]] <- list(
      name = hit[[2]], operator = if (length(hit) >= 3L) hit[[3]] else "",
      version = if (length(hit) >= 4L) hit[[4]] else ""
    )
  }
  dependencies
}

satisfies <- function(actual, operator, expected) {
  if (!nzchar(operator)) return(TRUE)
  comparison <- utils::compareVersion(actual, expected)
  switch(operator, ">=" = comparison >= 0, "<=" = comparison <= 0,
         ">" = comparison > 0, "<" = comparison < 0,
         "==" = comparison == 0, FALSE)
}

dependency_versions <- list()
for (dependency in parse_dependencies(description)) {
  name <- dependency$name
  actual <- if (name == "R") {
    as.character(getRversion())
  } else if (name %in% names(installed_versions)) {
    installed_versions[[name]]
  } else {
    NA_character_
  }
  if (is.na(actual)) stop("Missing Presto dependency; refusing dependency installation: ", name)
  if (!satisfies(actual, dependency$operator, dependency$version)) {
    stop(sprintf("Installed Presto dependency conflict: %s %s; requires %s %s", name, actual,
                 dependency$operator, dependency$version))
  }
  dependency_versions[[name]] <- actual
}

presto_before <- installed.packages(lib.loc = library_path)
presto_rows <- which(presto_before[, "Package"] == "presto")
status <- "installed"
if (length(presto_rows)) {
  current_version <- presto_before[presto_rows[[1]], "Version"]
  if (current_version != expected_package_version) {
    stop(sprintf("Target library already has Presto %s; refusing to replace it with %s", current_version, expected_package_version))
  }
  status <- "already_installed_same_version_skipped"
} else {
  utils::install.packages(archive_expected, lib = library_path, repos = NULL,
                          type = "source", dependencies = FALSE, quiet = FALSE)
}

presto_after <- installed.packages(lib.loc = library_path)
after_versions <- setNames(presto_after[, "Version"], presto_after[, "Package"])
changed <- intersect(names(presto_before[, "Version"]), names(after_versions))
changed <- changed[presto_before[match(changed, presto_before[, "Package"]), "Version"] != after_versions[changed]]
missing_before <- setdiff(presto_before[, "Package"], names(after_versions))
if (length(changed) || length(missing_before)) {
  stop("Existing r442_sc packages changed during Presto setup: ", paste(c(changed, missing_before), collapse = ", "))
}
if (!requireNamespace("presto", quietly = TRUE, lib.loc = library_path) ||
    as.character(utils::packageVersion("presto", lib.loc = library_path)) != expected_package_version) {
  stop("Post-setup verification failed for Presto 1.1.0")
}

result <- list(
  status = status,
  package = "presto",
  version = as.character(utils::packageVersion("presto", lib.loc = library_path)),
  library = library_path,
  source_archive = "research/tools/sources/presto-b5df6ee.tar.gz",
  source_commit = expected_commit,
  manifest_sha256 = entry$sha256,
  actual_source_sha256 = archive_sha256,
  compiler = Sys.which("x86_64-conda-linux-gnu-c++"),
  compiler_path_added_for_project_setup = runtime_bin,
  dependencies = dependency_versions,
  existing_package_versions_changed = FALSE,
  analysis_performed = FALSE
)
cat(jsonlite::toJSON(result, pretty = TRUE, auto_unbox = TRUE), "\n")
