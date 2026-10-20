#!/usr/bin/env Rscript

target_versions <- c(
  coloc = "5.2.3",
  viridis = "0.6.5",
  susieR = "0.14.2",
  mixsqp = "0.3-54",
  reshape = "0.8.10"
)
install_order <- c("mixsqp", "reshape", "viridis", "susieR", "coloc")
cran <- "https://cran.r-project.org"

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
library_path <- normalizePath(file.path(root, "research/envs/r442_sc/lib/R/library"), mustWork = TRUE)
records <- file.path(root, "research/records/advance_20261017")
cache <- file.path(root, "research/tools/coloc-cache")
dir.create(records, recursive = TRUE, showWarnings = FALSE)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
runtime_bin <- file.path(root, "research/envs/r442_sc/bin")
if (!file.exists(file.path(runtime_bin, "x86_64-conda-linux-gnu-c++"))) {
  stop("The r442_sc C++ compiler wrapper is missing; refusing to install source packages")
}
Sys.setenv(PATH = paste(runtime_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
.libPaths(c(library_path, .Library))
options(repos = c(CRAN = cran))

baseline_path <- file.path(records, "coloc_install_before.tsv")
current_target <- installed.packages(lib.loc = library_path)
current_target_table <- data.frame(Package = current_target[, "Package"], Version = current_target[, "Version"], stringsAsFactors = FALSE)
if (file.exists(baseline_path)) {
  baseline_table <- read.delim(baseline_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
} else {
  baseline_table <- current_target_table
  write.table(baseline_table, baseline_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
baseline_versions <- setNames(baseline_table$Version, baseline_table$Package)
current_versions <- setNames(current_target_table$Version, current_target_table$Package)
missing_baseline <- setdiff(names(baseline_versions), names(current_versions))
changed_baseline <- intersect(names(baseline_versions), names(current_versions))
changed_baseline <- changed_baseline[baseline_versions[changed_baseline] != current_versions[changed_baseline]]
if (length(missing_baseline) || length(changed_baseline)) {
  stop("An original target-library package is missing or changed: ", paste(c(missing_baseline, changed_baseline), collapse = ", "))
}
before_all <- installed.packages()
before_versions <- setNames(before_all[!duplicated(rownames(before_all)), "Version"], rownames(before_all))

if ("coloc" %in% current_target_table$Package && current_target_table$Version[match("coloc", current_target_table$Package)] != target_versions[["coloc"]]) {
  stop("The target library already contains coloc at a different version; refusing to replace it")
}

parse_dependencies <- function(description) {
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(description))
  text <- paste(description[1, fields], collapse = ",")
  if (!nzchar(text)) return(list())
  parts <- trimws(unlist(strsplit(gsub("\n", " ", text), ",", fixed = TRUE)))
  dependencies <- list()
  for (part in parts) {
    hit <- regmatches(part, regexec("^([A-Za-z][A-Za-z0-9.]*)\\s*(?:\\((>=|<=|==|>|<)\\s*([^)]*)\\))?$", part))[[1]]
    if (length(hit)) {
      dependencies[[length(dependencies) + 1L]] <- list(
        name = hit[[2]], operator = if (length(hit) >= 3L) hit[[3]] else "",
        version = if (length(hit) >= 4L) hit[[4]] else ""
      )
    }
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

download_source <- function(package, version) {
  archive <- file.path(cache, sprintf("%s_%s.tar.gz", package, version))
  if (file.exists(archive)) return(archive)
  urls <- c(
    sprintf("%s/src/contrib/%s_%s.tar.gz", cran, package, version),
    sprintf("%s/src/contrib/Archive/%s/%s_%s.tar.gz", cran, package, package, version)
  )
  for (url in urls) {
    ok <- tryCatch({
      utils::download.file(url, archive, mode = "wb", quiet = TRUE)
      TRUE
    }, warning = function(w) FALSE, error = function(e) FALSE)
    if (ok && file.exists(archive) && file.info(archive)$size > 0) return(archive)
    if (file.exists(archive)) unlink(archive)
  }
  stop("Could not download official CRAN source for ", package, " ", version)
}

read_description <- function(archive, package) {
  member <- paste0(package, "/DESCRIPTION")
  extracted <- tempfile("coloc-desc-")
  dir.create(extracted)
  on.exit(unlink(extracted, recursive = TRUE), add = TRUE)
  utils::untar(archive, files = member, exdir = extracted)
  description_path <- file.path(extracted, member)
  if (!file.exists(description_path)) stop("Source DESCRIPTION not found for ", package)
  read.dcf(description_path)
}

planned_packages <- setdiff(install_order, baseline_table$Package)
install_packages <- setdiff(install_order, current_target_table$Package)
archives <- setNames(character(length(planned_packages)), planned_packages)
descriptions <- list()
for (package in planned_packages) {
  version <- target_versions[[package]]
  archives[[package]] <- download_source(package, version)
  descriptions[[package]] <- read_description(archives[[package]], package)
  if (descriptions[[package]][1, "Version"] != version) {
    stop(package, " source archive version does not match pinned version ", version)
  }
}

planned_versions <- c(before_versions, target_versions)
plan_rows <- list()
for (package in planned_packages) {
  description <- descriptions[[package]]
  for (dependency in parse_dependencies(description)) {
    name <- dependency$name
    if (name == "R") {
      actual <- as.character(getRversion())
    } else if (name %in% names(before_versions)) {
      actual <- before_versions[[name]]
    } else if (name %in% names(target_versions)) {
      actual <- target_versions[[name]]
    } else {
      stop("Missing dependency not covered by the safe pinned plan: ", name, " (required by ", package, ")")
    }
    if (!satisfies(actual, dependency$operator, dependency$version)) {
      stop(sprintf("Dependency conflict: %s %s is installed/planned, but %s requires %s %s",
                   name, actual, package, dependency$operator, dependency$version))
    }
    plan_rows[[length(plan_rows) + 1L]] <- data.frame(
      package = package, package_version = target_versions[[package]],
      dependency = name, dependency_version = actual,
      constraint = paste(dependency$operator, dependency$version),
      dependency_source = if (name == "R") "R_runtime" else if (name %in% names(before_versions)) "already_installed" else "pinned_in_this_plan",
      stringsAsFactors = FALSE
    )
  }
}

plan_path <- file.path(records, "coloc_install_plan.tsv")
if (length(plan_rows)) {
  write.table(do.call(rbind, plan_rows), plan_path, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  writeLines("package\tpackage_version\tdependency\tdependency_version\tconstraint\tdependency_source", plan_path)
}

source_rows <- list()
for (package in planned_packages) {
  archive <- archives[[package]]
  source_rows[[length(source_rows) + 1L]] <- data.frame(
    package = package,
    version = target_versions[[package]],
    source_url = if (package == "coloc") paste0(cran, "/src/contrib/coloc_5.2.3.tar.gz") else paste0(cran, "/src/contrib/", basename(archive)),
    source_file = normalizePath(archive),
    sha256 = digest::digest(file = archive, algo = "sha256", serialize = FALSE),
    stringsAsFactors = FALSE
  )
}
source_table <- if (length(source_rows)) do.call(rbind, source_rows) else data.frame()
write.table(source_table, file.path(records, "coloc_install_sources.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (package in install_order) {
  if (package %in% install_packages) {
    utils::install.packages(archives[[package]], lib = library_path, repos = NULL,
                            type = "source", dependencies = FALSE, quiet = FALSE)
  }
}

after_target <- installed.packages(lib.loc = library_path)
after_versions <- setNames(after_target[, "Version"], after_target[, "Package"])
missing_old <- setdiff(names(baseline_versions), names(after_versions))
changed_old <- intersect(names(baseline_versions), names(after_versions))
changed_old <- changed_old[baseline_versions[changed_old] != after_versions[changed_old]]
if (length(missing_old) || length(changed_old)) {
  stop("An existing target-library package changed during install: ", paste(c(missing_old, changed_old), collapse = ", "))
}
if (!requireNamespace("coloc", quietly = TRUE, lib.loc = library_path) ||
    as.character(utils::packageVersion("coloc", lib.loc = library_path)) != target_versions[["coloc"]]) {
  stop("coloc 5.2.3 failed post-install namespace/version verification")
}

new_packages <- setdiff(after_target[, "Package"], baseline_table$Package)
new_table <- data.frame(Package = new_packages, Version = after_versions[new_packages], stringsAsFactors = FALSE)
write.table(new_table, file.path(records, "coloc_install_new.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Package = after_target[, "Package"], Version = after_target[, "Version"], stringsAsFactors = FALSE),
            file.path(records, "coloc_install_after.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary <- list(
  status = "installed_and_verified",
  r_version = as.character(getRversion()),
  library = library_path,
  coloc_version = as.character(utils::packageVersion("coloc", lib.loc = library_path)),
  newly_installed = as.list(setNames(new_table$Version, new_table$Package)),
  source_records = "research/records/advance_20261017/coloc_install_sources.tsv",
  dependency_policy = "Only pinned packages missing from the target library were installed; all existing target-library versions were checked unchanged.",
  analysis_performed = FALSE
)
jsonlite::write_json(summary, file.path(records, "coloc_install.json"), pretty = TRUE, auto_unbox = TRUE)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
