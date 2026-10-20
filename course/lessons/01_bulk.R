source("R/common.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
})

output <- lesson_output("01_bulk")
counts_file <- "data/bulk/counts.tsv.gz"
samples_file <- "data/bulk/samples.tsv"
annotation_file <- "data/bulk/gene_symbols.tsv"

fdr_cutoff <- 0.05
logfc_cutoff <- 1
candidate_symbols <- c("IL12B", "SGSH", "TNFAIP3", "DDX58", "IL7R", "PRSS53", "UBLCP1")

# 1. 读取原始整数计数并按样本顺序连接临床分组。
raw <- read_tsv(counts_file)
stopifnot(!anyDuplicated(raw$GeneID))
counts <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- as.character(raw$GeneID)
stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)))

metadata <- read_tsv(samples_file)
metadata$included_in_counts <- metadata$gsm %in% colnames(counts)
write_tsv(metadata, file.path(output, "sample_inclusion.tsv"))
samples <- metadata[match(colnames(counts), metadata$gsm), , drop = FALSE]
stopifnot(!anyNA(samples$gsm), !anyDuplicated(samples$gsm))
stopifnot(all(samples$tissue_type %in% c("normal skin", "lesional psoriatic skin")))
group <- factor(ifelse(samples$tissue_type == "normal skin", "normal", "psoriasis"),
                levels = c("normal", "psoriasis"))

# 2. 只移除全零基因；不增加低表达阈值。
keep <- rowSums(counts) > 0
y <- DGEList(counts = counts[keep, , drop = FALSE], group = group)

# 3. TMM 标准化、voom 均值-方差建模和 limma 经验贝叶斯检验。
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~ group)
pdf(file.path(output, "voom_mean_variance.pdf"), width = 7, height = 6)
v <- voom(y, design, plot = TRUE)
dev.off()
fit <- eBayes(lmFit(v, design))
result <- topTable(fit, coef = "grouppsoriasis", number = Inf,
                   adjust.method = "BH", sort.by = "none")
result$GeneID <- rownames(result)

# 4. 用当前 NCBI 符号注释结果；候选列表同时支持符号别名。
annotation <- read_tsv(annotation_file)
annotation$GeneID <- as.character(annotation$GeneID)
match_id <- match(result$GeneID, annotation$GeneID)
result$Symbol <- annotation$Symbol[match_id]
result$annotation_status <- ifelse(is.na(match_id), "unmapped", "current_NCBI_symbol")
result$significant <- result$adj.P.Val < fdr_cutoff & abs(result$logFC) >= logfc_cutoff
result$contrast <- "psoriasis - normal"
result$source_mode <- "NCBI_reprocessed_171_of_174"
write_tsv(result, file.path(output, "all_gene_differential_expression.tsv"))

sample_qc <- data.frame(gsm = samples$gsm, group = group,
                        library_size = y$samples$lib.size,
                        TMM_factor = y$samples$norm.factors)
write_tsv(sample_qc, file.path(output, "library_normalization.tsv"))
qc_plot <- ggplot(sample_qc, aes(x = group, y = library_size / 1e6, colour = group)) +
  geom_boxplot(width = .48, outlier.shape = NA, alpha = .12, show.legend = FALSE) +
  geom_jitter(width = .12, height = 0, alpha = .68, size = 1.5) +
  scale_colour_manual(values = c(normal = "#2C7FB8", psoriasis = "#C44536")) +
  labs(title = "Bulk sample library sizes", subtitle = "171 NCBI-reprocessed samples",
       x = NULL, y = "Library size (millions of reads)") + theme_course()
save_plot(qc_plot, output, "sample_qc", width = 6.5, height = 5)

volcano <- ggplot(result, aes(x = logFC, y = -log10(pmax(adj.P.Val, .Machine$double.xmin)))) +
  geom_point(aes(colour = significant), alpha = .62, size = 1.15) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = 2, colour = "#68737D") +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, colour = "#68737D") +
  scale_colour_manual(values = c(`FALSE` = "#AAB2B8", `TRUE` = "#B23A48"),
                      breaks = c(FALSE, TRUE), labels = c("Other genes", "Meets both cutoffs"),
                      name = "Differential expression") +
  labs(title = "Psoriatic lesional skin vs normal skin",
       subtitle = "All 37,065 tested genes; Benjamini-Hochberg adjustment",
       x = "log2 fold change (psoriasis - normal)", y = "-log10 adjusted P") + theme_course()
save_plot(volcano, output, "bulk_volcano", width = 7.4, height = 6)

candidate_rows <- lapply(candidate_symbols, function(target) {
  alias_hit <- !is.na(annotation$Synonyms) & vapply(annotation$Synonyms, function(x) {
    target %in% strsplit(x, "|", fixed = TRUE)[[1]]
  }, logical(1))
  ids <- unique(annotation$GeneID[annotation$Symbol == target | alias_hit])
  if (length(ids) != 1L) stop("候选基因无法唯一映射：", target)
  row <- result[result$GeneID == ids, , drop = FALSE]
  i <- match(ids, rownames(fit$coefficients))
  se <- fit$stdev.unscaled[i, "grouppsoriasis"] * sqrt(fit$s2.post[i])
  critical <- qt(.975, df = fit$df.total[i])
  data.frame(candidate = target, GeneID = ids, Symbol = row$Symbol,
             logFC = row$logFC, lower95 = row$logFC - critical * se,
             upper95 = row$logFC + critical * se, adj.P.Val = row$adj.P.Val,
             significant = row$significant)
})
candidates <- do.call(rbind, candidate_rows)
write_tsv(candidates, file.path(output, "candidate_effects.tsv"))
candidate_plot <- ggplot(candidates, aes(x = logFC, y = reorder(candidate, logFC), colour = significant)) +
  geom_vline(xintercept = 0, colour = "#68737D", linewidth = .4) +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y",
                width = .18, linewidth = .7) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(`FALSE` = "#2C7FB8", `TRUE` = "#B23A48"),
                      breaks = c(FALSE, TRUE), labels = c("Below both cutoffs", "Meets both cutoffs"),
                      name = "Differential expression") +
  labs(title = "Effect estimates for seven candidate genes",
       subtitle = "Moderated 95% confidence intervals; psoriasis - normal",
       x = "log2 fold change", y = NULL) + theme_course()
save_plot(candidate_plot, output, "candidate_effects", width = 7.2, height = 4.8)

cat(sprintf("samples=%d; normal=%d; psoriasis=%d; raw_genes=%d; tested_genes=%d; significant=%d\n",
            ncol(counts), sum(group == "normal"), sum(group == "psoriasis"), nrow(counts),
            sum(keep), sum(result$significant)))
