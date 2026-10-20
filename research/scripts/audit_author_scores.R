#!/usr/bin/env Rscript

score_slots <- c("sensitivity", "smr_blood", "smr_skin", "pgs_psoriasis",
                 "pgs_psa", "bulk", "single_cell", "coloc", "ppi")

score_slot <- function(value, reference, slot) {
  value <- trimws(as.character(value))
  reference <- trimws(as.character(reference))
  if (is.na(value) || !nzchar(value)) return(0L)
  if (slot %in% c("coloc", "ppi")) return(as.integer(grepl("√", value, fixed = TRUE)))
  if (reference == "+") {
    if (slot %in% c("bulk", "single_cell")) return(as.integer(grepl("↑", value, fixed = TRUE)))
    return(as.integer(grepl("+", value, fixed = TRUE)))
  }
  if (reference == "−") {
    if (slot %in% c("bulk", "single_cell")) return(as.integer(grepl("↓", value, fixed = TRUE)))
    return(as.integer(grepl("−", value, fixed = TRUE)))
  }
  stop("Invalid PWAS reference direction: ", reference)
}

parse_author_score <- function(value) {
  match <- regexec("\\(([0-9]+)\\)", as.character(value))
  parts <- regmatches(as.character(value), match)[[1]]
  if (length(parts) != 2L) stop("Could not parse author_tier_score: ", value)
  as.integer(parts[[2]])
}

recalculate_scores <- function(data) {
  required <- c("target", "pwas", score_slots, "author_tier_score")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing required author table columns: ", paste(missing, collapse = ", "))
  if (any(!data$pwas %in% c("+", "−"))) stop("Every row must have pwas direction '+' or '−'")
  slot_values <- vapply(seq_len(nrow(data)), function(i) {
    vapply(score_slots, function(slot) score_slot(data[[slot]][i], data$pwas[i], slot), integer(1))
  }, integer(length(score_slots)))
  slot_values <- t(slot_values)
  colnames(slot_values) <- score_slots
  author_score <- vapply(data$author_tier_score, parse_author_score, integer(1))
  recalculated <- rowSums(slot_values)
  per_slot <- data.frame(
    target = rep(data$target, each = length(score_slots)),
    slot = rep(score_slots, times = nrow(data)),
    contribution = as.vector(t(slot_values)),
    author_score = rep(author_score, each = length(score_slots)),
    recalculated_score = rep(recalculated, each = length(score_slots)),
    difference = rep(recalculated - author_score, each = length(score_slots)),
    evidence_origin = "author_reported_arithmetic_only",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(per_slot = per_slot, author_score = author_score,
       recalculated_score = recalculated, difference = recalculated - author_score)
}

sha256_file <- function(path) {
  result <- system2("sha256sum", args = shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!length(result) || (!is.null(status) && status != 0L)) {
    stop("sha256sum failed for ", path, ": ", paste(result, collapse = "\n"))
  }
  sub("[[:space:]].*$", "", result[[1]])
}

run_audit <- function(input_path, output_dir, script_path) {
  raw <- readBin(input_path, what = "raw", n = file.info(input_path)$size)
  input_sha256 <- sha256_file(input_path)
  data <- read.delim(text = rawToChar(raw), sep = "\t", check.names = FALSE,
                     stringsAsFactors = FALSE, quote = "", na.strings = character())
  result <- recalculate_scores(data)
  if (nrow(data) != 29L) stop("Expected 29 author rows; found ", nrow(data))
  if (any(result$difference != 0L)) {
    failed <- which(result$difference != 0L)
    details <- paste(sprintf("%s author=%d recalculated=%d difference=%d",
                             data$target[failed], result$author_score[failed],
                             result$recalculated_score[failed], result$difference[failed]),
                     collapse = "; ")
    stop("Author score arithmetic mismatch: ", details)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(result$per_slot, file.path(output_dir, "author_score_slots.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  summary <- list(
    evidence_origin = "author_reported_arithmetic_only",
    rows = nrow(data),
    slots = as.list(score_slots),
    matched_rows = sum(result$difference == 0L),
    mismatched_rows = sum(result$difference != 0L),
    author_score_range = as.list(range(result$author_score)),
    recalculated_score_range = as.list(range(result$recalculated_score)),
    input_sha256 = input_sha256,
    script_sha256 = sha256_file(script_path),
    output_sha256 = sha256_file(file.path(output_dir, "author_score_slots.tsv")),
    scope = paste("Arithmetic recheck of author-reported scores only; not evidence reconstruction",
                  "and not a source for adding UKB scores.")
  )
  jsonlite::write_json(summary, file.path(output_dir, "summary.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
  project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
  input_path <- if (length(args) >= 1L) args[[1]] else file.path(project_root, "research/references/author_main_table1.tsv")
  output_dir <- if (length(args) >= 2L) args[[2]] else file.path(project_root, "research/results/author_score_audit")
  run_audit(input_path, output_dir, script_path)
}
