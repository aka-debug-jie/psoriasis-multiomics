script_arg <- grep("^--file=", commandArgs(), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(test_path), "..", "scripts", "audit_author_scores.R"), local = TRUE)

fixture <- data.frame(
  target = "fixture",
  pwas = "+",
  sensitivity = "+/ −",
  smr_blood = "−",
  smr_skin = "+",
  pgs_psoriasis = "−/+",
  pgs_psa = "",
  bulk = "↑",
  single_cell = "↓",
  coloc = "√",
  ppi = "",
  author_tier_score = "Tier 1(5)",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
scored <- recalculate_scores(fixture)
stopifnot(identical(as.integer(scored$per_slot$contribution), c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 0L)))
stopifnot(identical(unname(scored$author_score), 5L))
stopifnot(unname(scored$recalculated_score) == 5)
stopifnot(unname(scored$difference) == 0)
cat("author score rule fixture: OK\n")
