# Same-series NCBI reprocessed counts; not the submitter's historical count matrix.
bulk_ncbi <- function(root, config, output) {
  settings <- config$bulk
  stopifnot(settings$gene_filter == "remove_all_zero_only", settings$design == "~ group",
            settings$contrast == "psoriasis - normal", settings$model == "limma_voom_lmFit_eBayes",
            settings$p_adjust == "BH_across_all_tested_genes")
  count_file <- file.path(root, settings$counts)
  raw <- data.table::fread(cmd = paste("gzip -dc", shQuote(count_file)), data.table = FALSE, check.names = FALSE)
  stopifnot(!anyDuplicated(raw$GeneID))
  counts <- as.matrix(raw[, -1]); rownames(counts) <- as.character(raw$GeneID)
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)))
  metadata <- read.delim(file.path(root, settings$sample_metadata), stringsAsFactors = FALSE)
  metadata$included_in_counts <- metadata$gsm %in% colnames(counts)
  write.table(metadata, file.path(output, "sample_inclusion.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  samples <- metadata[match(colnames(counts), metadata$gsm), ]
  stopifnot(!anyNA(samples$gsm), !anyDuplicated(samples$gsm),
            all(samples$tissue_type %in% c("normal skin", "lesional psoriatic skin")))
  group <- factor(ifelse(samples$tissue_type == "normal skin", "normal", "psoriasis"), levels = c("normal", "psoriasis"))
  # No undocumented low-expression cutoff: only all-zero rows are removed.
  keep <- rowSums(counts) > 0
  y <- edgeR::DGEList(counts[keep, , drop = FALSE], group = group)
  y <- edgeR::calcNormFactors(y, method = settings$normalization)
  design <- model.matrix(~ group)
  pdf(file.path(output, "voom_mean_variance.pdf"), width = 7, height = 6)
  v <- limma::voom(y, design, plot = TRUE)
  dev.off()
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  result <- limma::topTable(fit, coef = "grouppsoriasis", number = Inf, adjust.method = "BH", sort.by = "none")
  result$GeneID <- rownames(result)
  annotation <- read.delim(gzfile(file.path(root, settings$annotation)), check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  match_id <- match(result$GeneID, as.character(annotation$GeneID))
  result$Symbol <- annotation$Symbol[match_id]
  result$annotation_status <- ifelse(is.na(match_id), "unmapped", "current_NCBI_symbol")
  result$significant <- result$adj.P.Val < settings$fdr & abs(result$logFC) >= settings$abs_log2fc
  result$contrast <- "psoriasis - normal"
  result$source_mode <- settings$source_mode
  write.table(result, file.path(output, "all_gene_differential_expression.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  sample_qc <- data.frame(gsm = samples$gsm, group = group, library_size = y$samples$lib.size,
                          TMM_factor = y$samples$norm.factors)
  write.table(sample_qc, file.path(output, "library_normalization.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  author <- read.delim(file.path(root, "research/references/author_main_table1.tsv"), check.names = FALSE)
  target_comparison <- do.call(rbind, lapply(seq_len(nrow(author)), function(i) {
    target <- author$target[i]
    symbol <- if (target == "IL12") "IL12A" else if (target == "IL23") "IL23A" else target
    ids <- annotation$GeneID[annotation$Symbol == symbol]
    if (!length(ids)) ids <- annotation$GeneID[grepl(paste0("(^|\\|)", symbol, "(\\||$)"), annotation$Synonyms)]
    rows <- result[result$GeneID %in% as.character(ids), ]
    if (nrow(rows) != 1) return(data.frame(paper_target = target, requested_gene = symbol, GeneID = NA,
                                         Symbol = NA, logFC = NA, adj.P.Val = NA, significant = NA,
                                         paper_bulk = author$bulk[i], mapping = "missing_or_ambiguous"))
    data.frame(paper_target = target, requested_gene = symbol, rows[, c("GeneID", "Symbol", "logFC", "adj.P.Val", "significant")],
               paper_bulk = author$bulk[i], mapping = "unique_GeneID")
  }))
  write.table(target_comparison, file.path(output, "author_target_comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  pdf(file.path(output, "bulk_volcano.pdf"), width = 7, height = 6)
  plot(result$logFC, -log10(pmax(result$adj.P.Val, .Machine$double.xmin)), pch = 16, cex = .35,
       col = ifelse(result$significant, "#B23A48", "#777777"), xlab = "log2 fold change: psoriasis - normal",
       ylab = "-log10(BH-adjusted P)", main = sprintf("GSE54456: NCBI reprocessed counts (%d samples)", ncol(counts)))
  abline(v = c(-1, 1) * settings$abs_log2fc, h = -log10(settings$fdr), lty = 2, col = "grey40")
  dev.off()
  list(status = "same_study_reprocessed_bulk_completed", raw_genes = nrow(counts), tested_genes = sum(keep),
       samples = ncol(counts), psoriasis = sum(group == "psoriasis"), normal = sum(group == "normal"),
       missing_from_GEO = metadata$gsm[!metadata$included_in_counts],
       significant_genes = sum(result$significant), author_targets_significant = sum(target_comparison$significant, na.rm = TRUE),
       original_author_counts_recovered = FALSE, annotation = "current_NCBI_gene_symbols; count alignment remains NCBI Release109.20190905",
       settings = settings, versions = list(edgeR = as.character(packageVersion("edgeR")), limma = as.character(packageVersion("limma"))),
       generated_files = c("sample_inclusion.tsv", "all_gene_differential_expression.tsv", "library_normalization.tsv",
                           "author_target_comparison.tsv", "voom_mean_variance.pdf", "bulk_volcano.pdf"))
}
