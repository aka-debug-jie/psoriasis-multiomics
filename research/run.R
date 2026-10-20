#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."))
value <- function(key, default = NULL) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", key)
  args[i + 1]
}
module <- value("--module", "status")
config_path <- value("--config", file.path(root, "research/config/analysis.json"))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Use research/envs/r442/bin/Rscript after bootstrap_env.sh")
config <- jsonlite::fromJSON(config_path)
if (module == "status") {
  status_path <- file.path(root, "research/records/module_status.tsv")
  if (!file.exists(status_path)) stop("Module audit has not been recorded yet")
  print(read.delim(status_path), row.names = FALSE)
  quit(status = 0)
}
if (!module %in% c("supplement_check", "geo_audit", "sc_qc", "mr_available", "bulk_ncbi", "sc_embedding", "sc_de", "sc_pseudobulk")) stop("No runnable implementation for module: ", module, "; see --module status")
if (module %in% c("geo_audit", "sc_qc")) {
  python_script <- file.path(root, "research/scripts", if (module == "geo_audit") "audit_geo.py" else "sc_qc.py")
  python_args <- shQuote(python_script)
  if (module == "sc_qc") python_args <- c(python_args, "--config", shQuote(config_path))
  python <- file.path(root, "research/envs/data/bin/python")
  if (!file.exists(python)) stop("Run research/scripts/bootstrap_data_env.sh first")
  exit_status <- system2(python, python_args)
  if (exit_status != 0L) stop("Python preprocessing failed for ", module, "; exit ", exit_status)
  quit(status = 0)
}
source(file.path(root, "research/scripts/mr_core.R"))
source(file.path(root, "research/scripts", paste0(module, ".R")))
output <- value("--output", file.path(root, "research/results", module))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
jsonlite::write_json(list(status = "running", module = module, started_utc = started),
                     file.path(output, "run.json"), pretty = TRUE, auto_unbox = TRUE)
set.seed(config$random_seed)
sha256 <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " +")[[1]][1]
sources <- c(config_path, script, file.path(root, "research/scripts/mr_core.R"),
             file.path(root, "research/scripts", paste0(module, ".R")))
if (module == "supplement_check") sources <- c(sources, file.path(root, "research/references", c(paste0("author_table_S", c(2, 5, 6), ".tsv"), "author_main_table1.tsv")))
if (module == "mr_available") sources <- c(sources, file.path(root, "research/data/processed", c("exposures/eldjarn_all_cis_unselected.tsv", "outcomes/mvp_normalized.tsv", "outcomes/finngen_r12_normalized.tsv")))
stuart_input <- file.path(root, "research/data/processed/outcomes/stuart_normalized.tsv")
if (module == "mr_available" && file.exists(stuart_input)) sources <- c(sources, stuart_input)
if (module == "bulk_ncbi") sources <- c(sources, file.path(root, c(config$bulk$counts, config$bulk$annotation,
                                                                 config$bulk$sample_metadata, "research/references/author_main_table1.tsv")))
if (module == "sc_embedding") sources <- c(sources, file.path(root, c(config$approved_parameters,
  config$sc$counts, config$sc$features, config$sc$cells, config$sc$input_manifest)))
if (module == "sc_de") sources <- c(sources, file.path(root, c("research/scripts/sc_de_helpers.R",
  config$sc_de$checkpoint, config$sc_de$embedding_manifest, config$sc_de$features,
  config$sc_de$annotation, config$sc_de$annotation_evidence, config$sc_de$author_targets)))
if (module == "sc_pseudobulk") sources <- c(sources, file.path(root, c(config$pseudobulk$counts,
  config$pseudobulk$samples, config$pseudobulk$features, config$pseudobulk$cell_DE)))
result <- tryCatch(get(module)(root, config, output), error = function(e) list(status = "failed", error = conditionMessage(e)))
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
outputs <- if (identical(result$status, "failed")) character() else file.path(output, c(result$generated_files, "sessionInfo.txt"))
manifest <- list(module = module, started_utc = started,
                 finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                 result = result, config = config,
                 inputs = setNames(as.list(vapply(sources, sha256, character(1))), sources),
                 outputs = setNames(as.list(vapply(outputs, sha256, character(1))), basename(outputs)))
jsonlite::write_json(manifest, file.path(output, "run.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
print(result)
if (identical(result$status, "failed")) quit(status = 1)
