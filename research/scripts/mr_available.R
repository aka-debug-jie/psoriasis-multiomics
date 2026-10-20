# Independent calculations on available sources, kept separate from the formal
# three-cohort discovery until its data coverage and selection rules are resolved.
mr_available <- function(root, config, output) {
  unique_effect <- function(d, key, columns, name) {
    groups <- split(seq_len(nrow(d)), d[[key]])
    conflicting <- unlist(lapply(groups, function(i) if (nrow(unique(d[i, columns, drop = FALSE])) > 1) i else integer()))
    write.table(d[conflicting, , drop = FALSE], file.path(output, paste0(name, "_ambiguous_matches.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    if (length(conflicting)) d <- d[-conflicting, , drop = FALSE]
    d[!duplicated(d[[key]]), , drop = FALSE]
  }
  x <- read.delim(file.path(root, "research/data/processed/exposures/eldjarn_all_cis_unselected.tsv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
  x$exposure_row_id <- seq_len(nrow(x))
  # Expand only explicit aliases. Rows lacking rsIDs are not assigned proxies.
  x <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    ids <- strsplit(x$variant[i], ",", fixed = TRUE)[[1]]
    ids <- ids[grepl("^rs[0-9]+$", ids)]
    if (!length(ids)) return(NULL)
    row <- x[rep(i, length(ids)), ]; row$rsid <- ids; row
  }))
  y <- read.delim(file.path(root, "research/data/processed/outcomes/mvp_normalized.tsv"), stringsAsFactors = FALSE)
  y <- y[as.logical(y$usable_effect), ]
  paired <- merge(x, y, by = "rsid", suffixes = c(".exposure", ".outcome"))
  orientation <- lapply(seq_len(nrow(paired)), function(i) with(paired[i, ],
    orient_alleles(EA.exposure, OA.exposure, EA.outcome, OA.outcome, beta.outcome, eaf.outcome)))
  paired$harmonization <- vapply(orientation, `[[`, character(1), "status")
  paired$beta_outcome_aligned <- vapply(orientation, `[[`, numeric(1), "beta")
  valid <- paired[is.finite(paired$beta_outcome_aligned) & is.finite(paired$se_reconstructed) &
                    paired$se_reconstructed > 0 & paired$beta.exposure != 0 & !is.na(paired$uniprot), ]
  valid <- unique_effect(valid, "exposure_row_id", c("beta_outcome_aligned", "se", "chr.outcome", "pos.outcome"), "mvp")
  # Freeze selection before examining outcome significance: strongest published
  # marginal cis association that has an unambiguous usable outcome match.
  valid <- valid[order(-valid$neg_log10_p, valid$assay, valid$rsid), ]
  valid <- valid[!duplicated(valid[c("cohort", "uniprot")]), ]
  effects <- wald_ratio(valid$beta.exposure, valid$se_reconstructed, valid$beta_outcome_aligned, valid$se)
  names(effects)[1:2] <- c("beta_mr", "se_mr")
  per_cohort <- cbind(valid, effects)
  per_cohort$evidence_origin <- "independent_available_source_reconstruction"
  per_cohort$analysis_status <- "partial_coverage_not_formal_discovery"
  write.table(per_cohort, file.path(output, "cohort_wald.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  combined <- do.call(rbind, lapply(split(per_cohort, per_cohort$uniprot), function(g) {
    cbind(uniprot = g$uniprot[1], gene = paste(unique(g$gene), collapse = ";"),
          cohorts = paste(g$cohort, collapse = ";"), fixed_effect(g$beta_mr, g$se_mr))
  }))
  combined$passes_paper_reported_threshold <- combined$p < config$discovery$reported_p_threshold
  combined$observed_coverage_bonferroni_p <- .05 / nrow(combined)
  combined$passes_observed_coverage_threshold <- combined$p < .05 / nrow(combined)
  combined$analysis_status <- "partial_coverage_not_formal_discovery"
  write.table(combined, file.path(output, "available_meta.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(as.data.frame(table(paired$harmonization)), file.path(output, "harmonization_counts.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  instruments <- valid[c("cohort", "assay", "gene", "uniprot", "rsid", "EA.exposure", "OA.exposure", "beta.exposure", "se_reconstructed", "exposure_row_id")]
  replication_counts <- list()
  generated <- c("cohort_wald.tsv", "available_meta.tsv", "harmonization_counts.tsv", "mvp_ambiguous_matches.tsv", "available_source_volcano.pdf")
  for (outcome in c("finngen_r12", "stuart")) {
  input <- file.path(root, "research/data/processed/outcomes", paste0(outcome, "_normalized.tsv"))
  if (!file.exists(input)) next
  prefix <- if (outcome == "finngen_r12") "finngen" else "stuart"
  fy <- read.delim(input, stringsAsFactors = FALSE)
  fy <- fy[as.logical(fy$usable_effect), ]
  rep <- merge(instruments, fy, by = "rsid")
  rep_orientation <- lapply(seq_len(nrow(rep)), function(i) with(rep[i, ],
    orient_alleles(EA.exposure, OA.exposure, EA, OA, beta, eaf)))
  rep$harmonization <- vapply(rep_orientation, `[[`, character(1), "status")
  rep$beta_outcome_aligned <- vapply(rep_orientation, `[[`, numeric(1), "beta")
  rep <- rep[is.finite(rep$beta_outcome_aligned), ]
  rep <- unique_effect(rep, "exposure_row_id", c("beta_outcome_aligned", "se", "chr", "pos"), prefix)
  rep_effects <- wald_ratio(rep$beta.exposure, rep$se_reconstructed, rep$beta_outcome_aligned, rep$se)
  names(rep_effects)[1:2] <- c("beta_mr", "se_mr")
  rep <- cbind(rep, rep_effects)
  write.table(rep, file.path(output, paste0(prefix, "_same_instrument_wald.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  rep_meta <- do.call(rbind, lapply(split(rep, rep$uniprot), function(g) {
    cbind(uniprot = g$uniprot[1], gene = paste(unique(g$gene), collapse = ";"), fixed_effect(g$beta_mr, g$se_mr))
  }))
  write.table(rep_meta, file.path(output, paste0(prefix, "_available_meta.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  replication_counts[[outcome]] <- nrow(rep_meta)
  generated <- c(generated, paste0(prefix, c("_ambiguous_matches.tsv", "_same_instrument_wald.tsv", "_available_meta.tsv")))
  }
  pdf(file.path(output, "available_source_volcano.pdf"), width = 7, height = 6)
  plot(combined$beta, combined$neg_log10_p, pch = 16, cex = .5,
       xlab = "Combined log odds ratio", ylab = "-log10(P)",
       main = "Available deCODE and UKB-PPP inputs")
  abline(h = -log10(config$discovery$reported_p_threshold), lty = 2, col = "red")
  mtext("Partial coverage; formal discovery filters not completed", side = 3, cex = .8)
  dev.off()
  list(status = "partial_available_data_calculated", proteins = nrow(combined), cohort_estimates = nrow(per_cohort),
       missing = c("Fenland full exposure inputs", "author's full IV selection manifest", "MHC/Steiger formal inclusion decisions", if (length(replication_counts) < 2) "complete replication"),
       formal_candidate_list_produced = FALSE,
       finngen_proteins = replication_counts$finngen_r12,
       replication_proteins = replication_counts, generated_files = generated,
       caveats = c("Palindromic and non-SNV alleles unresolved and excluded", "MVP SE may be inferred from rounded confidence intervals"))
}
