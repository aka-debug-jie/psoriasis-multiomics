source("R/common.R")
source("R/genetics.R")
suppressPackageStartupMessages(library(ggplot2))

output <- lesson_output("11_clinical_mr")
family <- read_tsv("data/clinical/family_17.tsv")
endpoints <- read_tsv("data/clinical/endpoints.tsv")
exposures <- read_tsv("data/genetics/exposures.tsv")
formal_candidates <- read_tsv("data/genetics/primary_candidates.tsv")
ld <- load_ld()

stopifnot(length(unique(family$uniprot)) == 17L,
          nrow(endpoints) == 22L,
          nrow(formal_candidates) == 7L,
          all(as.character(formal_candidates$primary_selected) == "True"))
exposures <- exposures[exposures$uniprot %in% family$uniprot, , drop = FALSE]
stopifnot(setequal(unique(exposures$uniprot), family$uniprot))
plan <- freeze_instruments(exposures, ld, max_abs_r = .95)
write_tsv(plan, file.path(output, "frozen_instruments.tsv"))

fits <- vector("list", nrow(endpoints))
harmonized <- vector("list", nrow(endpoints))
for (i in seq_len(nrow(endpoints))) {
  code <- endpoints$phenocode[i]
  y <- read_tsv(file.path("data/clinical", paste0(code, ".tsv")))
  fit <- analyse_outcome(plan, y, ld)
  fit$results$phenocode <- code
  fit$harmonized$phenocode <- code
  fits[[i]] <- fit$results
  harmonized[[i]] <- fit$harmonized
}

# The original multiplicity family is fixed before selecting seven formal candidates.
all_results <- do.call(rbind, fits)
stopifnot(nrow(all_results) == 374L,
          setequal(unique(all_results$phenocode), endpoints$phenocode),
          setequal(unique(all_results$uniprot), family$uniprot))
all_results$BH_original_family_374 <- p.adjust(all_results$p, method = "BH", n = 374L)
all_results <- merge(all_results, endpoints, by = "phenocode", all.x = TRUE,
                     sort = FALSE, suffixes = c("", "_endpoint"))
all_results <- all_results[match(paste(rep(family$uniprot, each = nrow(endpoints)),
                                        rep(endpoints$phenocode, times = nrow(family))),
                                 paste(all_results$uniprot, all_results$phenocode)), ]
write_tsv(all_results, file.path(output, "clinical_all_374.tsv"))

formal_ids <- unique(formal_candidates$uniprot)
candidate_results <- all_results[all_results$uniprot %in% formal_ids, , drop = FALSE]
stopifnot(nrow(candidate_results) == 154L,
          setequal(unique(candidate_results$uniprot), formal_ids))
candidate_results$BH_candidate_family_exploratory <- p.adjust(
  candidate_results$p, method = "BH", n = 154L)
candidate_results$candidate_family_size <- 154L
write_tsv(candidate_results, file.path(output, "clinical_candidates_154.tsv"))
write_tsv(do.call(rbind, harmonized), file.path(output, "harmonized_374.tsv"))

plot_data <- candidate_results
plot_data$gene <- factor(plot_data$gene,
                         levels = rev(formal_candidates$gene[order(formal_candidates$gene)]))
plot_data$phenocode <- factor(plot_data$phenocode, levels = endpoints$phenocode)
endpoint_labels <- vapply(sub(" [(]PheCode.*", "", endpoints$phenotype), function(x) paste(strwrap(x, width=38), collapse="\n"), "")
names(endpoint_labels) <- endpoints$phenocode
plot_data$marker <- ifelse(!is.na(plot_data$BH_original_family_374) &
                             plot_data$BH_original_family_374 < .05, "BH 374 < 0.05",
                           ifelse(!is.na(plot_data$BH_candidate_family_exploratory) &
                                    plot_data$BH_candidate_family_exploratory < .05,
                                  "Exploratory only", "Not below either threshold"))
plot <- ggplot(plot_data, aes(x = phenocode, y = gene, fill = beta)) +
  geom_tile(colour = "white", linewidth = .45) +
  geom_point(data = plot_data[plot_data$marker != "Not below either threshold", ], aes(shape=marker), size=1.8,colour="#263746") + scale_shape_manual(values=c("BH 374 < 0.05"=16,"Exploratory only"=1),name=NULL) +
  scale_fill_gradient2(low = "#2878A5", mid = "#F7F8F8", high = "#B64D4D",
                       midpoint = 0, na.value = "#D9E2E7", name = "Log odds ratio") +
  scale_x_discrete(labels = endpoint_labels) + coord_flip() +
  labs(title = "Clinical associations across 22 outcomes",
       subtitle = "17-protein BH is primary; open circles mark seven-candidate exploratory BH only",
       x = NULL, y = NULL) + theme_course() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, size = 9),
        panel.grid = element_blank(), legend.position = "bottom")
save_plot(plot, output, "clinical_effect_heatmap", width = 10, height = 9.5)

cat(sprintf("Full family: %d; seven-candidate comparisons: %d; BH374 significant: %d\n",
 nrow(all_results),nrow(candidate_results),sum(candidate_results$BH_original_family_374<.05,na.rm=TRUE)))
