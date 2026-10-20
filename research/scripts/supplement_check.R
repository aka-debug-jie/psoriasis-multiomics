supplement_check <- function(root, config, output) {
  input <- file.path(root, "research/references/author_table_S2.tsv")
  x <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(x[c("Protein", "SNP", "dat")]),
            !anyDuplicated(x[c("Protein", "dat")]))
  estimates <- wald_ratio(x$beta.exposure, x$se.exposure, x$beta.outcome, x$se.outcome)
  per_cohort <- cbind(x, estimates)
  per_cohort$evidence_origin <- "recalculated_from_author_selected_summary_rows"
  per_cohort$allele_provenance <- "EA/NEA supplied already harmonized in Table S2; not independently verified"
  write.table(per_cohort, file.path(output, "wald_by_cohort.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  combined <- do.call(rbind, lapply(split(per_cohort, per_cohort$Protein), function(g) {
    cbind(protein = g$Protein[1], cohorts = paste(g$dat, collapse = ";"), fixed_effect(g$beta, g$se))
  }))
  rownames(combined) <- NULL
  combined$evidence_origin <- "arithmetic_check_only_not_discovery"
  main <- read.delim(file.path(root, "research/references/author_main_table1.tsv"), check.names = FALSE)
  main_key <- ifelse(combined$protein == "IL12A", "IL12", combined$protein)
  reported <- main$pwas[match(main_key, main$target)]
  expected_sign <- ifelse(reported == "+", 1, ifelse(grepl("−", reported, fixed = TRUE), -1, NA_real_))
  combined$author_main_direction <- reported
  combined$direction_agrees_with_main_table <- sign(combined$beta) == expected_sign
  combined$passes_reported_discovery_threshold <- combined$p < config$discovery$reported_p_threshold
  write.table(combined, file.path(output, "fixed_effect_by_protein.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  s5 <- read.delim(file.path(root, "research/references/author_table_S5.tsv"), check.names = FALSE)
  eligible <- abs(s5$logFC) >= .25 & s5$adjusted_pvalue < .05
  s5$passes_paper_threshold_using_adjusted_pvalue <- eligible
  s5$BH_recalculated_over_26_reported_rows <- p.adjust(s5$p.value, method = "BH")
  s5$evidence_origin <- "author_reported_reclassified_not_recomputed_expression"
  write.table(s5, file.path(output, "author_sc_threshold_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  s6 <- read.delim(file.path(root, "research/references/author_table_S6.tsv"), check.names = FALSE)
  pp <- as.matrix(s6[paste0("PP.H", 0:4, ".abf")])
  stopifnot(all(is.finite(pp)), all(pp >= 0 & pp <= 1))
  s6$posterior_sum <- rowSums(pp)
  s6$evidence_origin <- "author_reported_probability_audit_not_coloc_rerun"
  write.table(s6, file.path(output, "author_coloc_probability_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  pdf(file.path(output, "arithmetic_check_forest.pdf"), width = 8, height = 9)
  par(mar = c(6.5, 7, 4, 2))
  n <- nrow(combined)
  plot(combined$beta, seq_len(n), xlim = range(log(c(combined$CI_low, combined$CI_high))),
       yaxt = "n", ylab = "", xlab = "Log odds ratio per SD (author-selected inputs)",
       main = "Table S2: Wald ratios and fixed-effect arithmetic check", pch = 19)
  segments(log(combined$CI_low), seq_len(n), log(combined$CI_high), seq_len(n))
  axis(2, at = seq_len(n), labels = combined$protein, las = 2, cex.axis = .8)
  abline(v = 0, lty = 2, col = "grey50")
  mtext("Not an independent discovery or replication analysis", side = 1, line = 5, cex = .8)
  dev.off()
  list(status = "calculation_checked", rows = nrow(x), proteins = nrow(combined),
       cohorts = as.list(setNames(as.integer(table(x$dat)), names(table(x$dat)))), min_F = min(estimates$F_statistic),
       sc_table_rows = nrow(s5), sc_table_genes = length(unique(s5$Gene)),
       sc_threshold_rows = sum(eligible), sc_threshold_genes = length(unique(s5$Gene[eligible])),
       s2_main_direction_conflicts = combined$protein[which(!combined$direction_agrees_with_main_table)],
       s2_meta_not_passing_main_threshold = combined$protein[which(!combined$passes_reported_discovery_threshold)],
       max_coloc_posterior_rounding_error = max(abs(rowSums(pp) - 1)),
       independent_reproduction = FALSE,
       generated_files = c("wald_by_cohort.tsv", "fixed_effect_by_protein.tsv", "author_sc_threshold_audit.tsv", "author_coloc_probability_audit.tsv", "arithmetic_check_forest.pdf"),
       missing = c("full protein universe", "independent exposure/outcome harmonization", "numeric source for published forest comparison"))
}
