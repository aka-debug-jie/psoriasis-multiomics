#!/usr/bin/env Rscript

# Small, explicit dispatcher for the existing research modules.
args <- commandArgs(trailingOnly = TRUE)
script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_arg), ".."), mustWork = TRUE)
rscript <- file.path(root, "research/envs/r442_sc/bin/Rscript")
python <- file.path(root, "research/envs/data/bin/python")
smr_tissues <- c("eQTLGen_Blood", "GTEx_Skin_Not_Sun_Exposed_Suprapubic",
                 "GTEx_Skin_Sun_Exposed_Lower_leg", "GTEx_Whole_Blood")
old_modules <- c("status", "supplement_check", "geo_audit", "sc_qc", "mr_available",
                 "bulk_ncbi", "sc_embedding", "sc_de", "sc_pseudobulk")
modules <- c(old_modules, "mr_substitute", "steiger", "smr", "coloc_mrna", "coloc_protein",
             "multi_iv", "functional", "drugs", "clinical", "covariance_resolution", "final_research")

parse_args <- function(x) {
  result <- list(module = NULL, resume = FALSE, list = FALSE)
  i <- 1L
  while (i <= length(x)) {
    if (x[i] == "--module") {
      if (i == length(x) || startsWith(x[i + 1L], "--")) stop("--module requires a name")
      result$module <- x[i + 1L]; i <- i + 2L
    } else if (x[i] == "--resume") {
      result$resume <- TRUE; i <- i + 1L
    } else if (x[i] == "--list") {
      result$list <- TRUE; i <- i + 1L
    } else stop("Unknown argument: ", x[i])
  }
  result
}

command_for <- function(module) {
  if (module %in% old_modules) {
    old_args <- c(file.path(root, "research/run.R"), "--module", module)
    configs<-c(bulk_ncbi='completion.json',sc_embedding='sc_analysis.json',sc_de='sc_de.json',sc_pseudobulk='sc_pseudobulk.json')
    if(module %in% names(configs)) old_args <- c(old_args,'--config',file.path(root,'research/config',configs[[module]]))
    program <- if(module %in% c('sc_embedding','sc_de','sc_pseudobulk')) rscript else file.path(root,'research/envs/r442/bin/Rscript')
    return(list(program=program,args=old_args))
  }
  if (module == "mr_substitute") return(list(program = rscript, args = file.path(root, "research/scripts/run_substitute_mr.R")))
  if (module == "steiger") return(list(program = rscript, args = file.path(root, "research/scripts/run_steiger_sensitivity.R")))
  if (module == "smr") return(lapply(smr_tissues, function(tissue)
    list(program = python, args = c(file.path(root, "research/scripts/run_smr.py"), "--tissue", tissue))))
  if (module == "coloc_mrna") return(list(program = rscript, args = c(file.path(root, "research/scripts/run_coloc.R"),
    "--config", file.path(root, "research/config/coloc_mrna_substitute.json"))))
  if (module == "coloc_protein") return(list(program = rscript, args = c(file.path(root, "research/scripts/run_protein_coloc.R"), "--substitute")))
  if (module == "multi_iv") return(list(program = rscript, args = c(file.path(root, "research/scripts/run_multi_instrument.R"), "--substitute")))
  if (module == "functional") return(list(program = python, args = c(file.path(root, "research/scripts/run_functional_enrichment.py"), "--scope", "substitute")))
  if (module == "drugs") return(list(program = python, args = file.path(root, "research/scripts/run_drug_network.py")))
  if (module == "clinical") return(list(program = rscript, args = file.path(root, "research/scripts/run_clinical_mr.R")))
  if (module == "covariance_resolution") return(lapply(c("resolve_shared_outcome.py", "diagnose_mvp_conventions.py", "assess_givw_heterogeneity.py"),
    function(script) list(program = python, args = file.path(root, "research/scripts", script))))
  if (module == "final_research") return(list(list(program=python,args=file.path(root,"research/scripts/finalize_research.py")),
    list(program=rscript,args=file.path(root,"research/scripts/plot_final_research.R"))))
  stop("Unknown module: ", module)
}

result_dirs <- function(module) {
  if (module == "mr_substitute") return(file.path(root, "research/results/mr_fenland_substitute"))
  if (module == "steiger") return(file.path(root, "research/results/steiger_substitute"))
  if (module == "smr") return(file.path(root, "research/results/smr", smr_tissues))
  if (module == "coloc_mrna") return(file.path(root, "research/results/coloc_mrna_substitute"))
  if (module == "coloc_protein") return(file.path(root, "research/results/coloc_protein_substitute"))
  if (module == "multi_iv") return(file.path(root, "research/results/multi_instrument_substitute"))
  if (module == "functional") return(file.path(root, "research/results/functional_substitute"))
  if (module == "drugs") return(file.path(root, "research/results/drug_network_substitute"))
  if (module == "clinical") return(file.path(root, "research/results/clinical_mr_substitute"))
  if (module == "covariance_resolution") return(file.path(root, "research/results", c("covariance_resolution", "mvp_convention_resolution", "covariance_resolution/heterogeneity")))
  if (module == "final_research") return(file.path(root,c("research/final","research/final/figures")))
  if (module %in% old_modules) {
    if (module == "status") return(character())
    if (module == "supplement_check") return(file.path(root, "research/results/supplement_check"))
    return(file.path(root, "research/results", module))
  }
  character()
}

hash_file <- function(path) {
  output <- tryCatch(system2("sha256sum", c("--", shQuote(path)), stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (length(output) != 1L || !grepl("^[0-9a-fA-F]{64}\\s", output[1])) return(NA_character_)
  sub("\\s.*$", "", output[1])
}

hash_entries <- function(x) {
  if (!is.list(x) || length(x) == 0L || is.null(names(x)) || any(!nzchar(names(x)))) return(NULL)
  values <- vapply(x, function(value) {
    if (length(value) != 1L || is.null(value)) return(NA_character_)
    as.character(value)
  }, character(1))
  if (anyNA(values) || any(!grepl("^[0-9a-fA-F]{64}$", values))) return(NULL)
  values
}

resolve_entry <- function(name, kind, result_dir) {
  if (grepl("^/", name)) return(name)
  candidates <- if (kind == "output") c(file.path(result_dir, name), file.path(root, name)) else {
    if (startsWith(name, "research/")) c(file.path(root, name), file.path(root, sub("^research/", "", name)))
    else c(file.path(root, "research", name), file.path(root, name))
  }
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[1] else candidates[1]
}

check_manifest <- function(module, manifest_path) {
  if (!file.exists(manifest_path)) return(list(ok = FALSE, reason = paste("missing manifest", manifest_path)))
  manifest <- tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = FALSE), error = function(e) NULL)
  status <- if(!is.null(manifest$status))manifest$status else manifest$result$status
  if (is.null(manifest) || is.null(status) || status %in% c("running", "failed", "incomplete"))
    return(list(ok = FALSE, reason = paste("invalid or unfinished manifest", manifest_path)))
  if(!is.null(manifest$script_sha256)) {
    command <- command_for(module)
    if(is.null(command$program)) command <- command[[1]]
    if(!identical(hash_file(command$args[[1]]),manifest$script_sha256)) return(list(ok=FALSE,reason='analysis script SHA256 mismatch'))
  }
  input_hashes <- hash_entries(manifest$inputs)
  output_hashes <- hash_entries(manifest$outputs)
  if (is.null(input_hashes) || is.null(output_hashes))
    return(list(ok = FALSE, reason = paste("manifest lacks named SHA256 inputs/outputs", manifest_path)))
  dir <- dirname(manifest_path)
  for (name in names(input_hashes)) {
    path <- resolve_entry(name, "input", dir)
    observed <- if (file.exists(path)) hash_file(path) else NA_character_
    if (is.na(observed) || !identical(tolower(observed), tolower(input_hashes[[name]])))
      return(list(ok = FALSE, reason = paste("input hash missing/mismatch", name)))
  }
  for (name in names(output_hashes)) {
    path <- resolve_entry(name, "output", dir)
    observed <- if (file.exists(path)) hash_file(path) else NA_character_
    if (is.na(observed) || !identical(tolower(observed), tolower(output_hashes[[name]])))
      return(list(ok = FALSE, reason = paste("output hash missing/mismatch", name)))
  }
  list(ok = TRUE, reason = "all named input/output SHA256 values match")
}

resume_status <- function(module) {
  dirs <- result_dirs(module)
  if (length(dirs) == 0L) return(list(ok = FALSE, reason = "module has no configured result manifest"))
  checks <- lapply(dirs, function(dir) check_manifest(module, file.path(dir, "run.json")))
  bad <- which(!vapply(checks, `[[`, logical(1), "ok"))
  if (length(bad)) return(list(ok = FALSE, reason = paste(basename(dirs[bad[1]]), checks[[bad[1]]]$reason)))
  list(ok = TRUE, reason = paste(length(checks), "manifest(s): all named input/output SHA256 values match"))
}

append_log <- function(path, lines) {
  cat(paste0(lines, collapse = "\n"), "\n", file = path, append = TRUE, sep = "")
}

run_command <- function(module, command, log_path) {
  args <- vapply(command$args, shQuote, character(1))
  line <- paste(shQuote(command$program), paste(args, collapse = " "))
  append_log(log_path, c(paste0("[", format(Sys.time(), tz = "UTC", usetz = TRUE), "] module=", module), line))
  message("Running module ", module)
  captured <- tryCatch(system2(command$program, args, stdout = TRUE, stderr = TRUE), error = function(e) {
    structure(conditionMessage(e), status = 127L)
  })
  append_log(log_path, as.character(captured))
  status <- attr(captured, "status")
  if (is.null(status)) status <- 0L
  if (length(status) != 1L || is.na(status)) status <- 1L
  as.integer(status)
}

parsed <- parse_args(args)
if (parsed$list) {
  if (!is.null(parsed$module) || parsed$resume) stop("--list cannot be combined with --module or --resume")
  cat("module\tcommand\n")
  for (module in modules) {
    cmd <- command_for(module)
    if (is.list(cmd) && !is.null(cmd$program)) {
      lines <- paste(basename(cmd$program), paste(basename(cmd$args), collapse = " "))
      description <- lines
    } else {
      lines <- vapply(cmd, function(x) paste(basename(x$program), paste(basename(x$args), collapse = " ")), character(1))
      description <- paste(unique(lines), collapse = " ; ")
    }
    cat(module, "\t", description, "\n", sep = "")
  }
  quit(status = 0L)
}
if (is.null(parsed$module)) stop("Usage: Rscript research/execute.R --module NAME [--resume] or --list")
module <- parsed$module
if (!module %in% modules) stop("Unknown module: ", module)

log_dir <- file.path(root, "research/records/dispatch")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(log_dir, paste0(module, ".log"))
if (parsed$resume) {
  check <- resume_status(module)
  append_log(log_path, c(paste0("[", format(Sys.time(), tz = "UTC", usetz = TRUE), "] module=", module, " resume=true"), check$reason))
  if (check$ok) {
    message("Skipping module ", module, " on --resume: ", check$reason)
    quit(status = 0L)
  }
  message("Resume verification failed for ", module, "; running: ", check$reason)
}

commands <- command_for(module)
if (!is.null(commands$program)) commands <- list(commands)
status <- 0L
for (command in commands) {
  status <- run_command(module, command, log_path)
  if (status != 0L) break
}
quit(status = status)
