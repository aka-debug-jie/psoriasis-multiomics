source("R/common.R")
suppressPackageStartupMessages(library(ggplot2))

output <- lesson_output("12_evidence")
score_slots <- c("sensitivity", "smr_blood", "smr_skin", "pgs_psoriasis",
                 "pgs_psa", "bulk", "single_cell", "coloc", "ppi")

# This is the recorded author-table rule. The author total is parsed for comparison only.
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

author <- read_tsv("data/evidence/author_score_input.tsv")
stopifnot(nrow(author) == 29L, all(author$pwas %in% c("+", "−")))
slot_matrix <- t(vapply(seq_len(nrow(author)), function(i) {
  vapply(score_slots, function(slot) score_slot(author[[slot]][i], author$pwas[i], slot), integer(1))
}, integer(length(score_slots))))
colnames(slot_matrix) <- score_slots
author_reported <- vapply(author$author_tier_score, parse_author_score, integer(1))
recalculated <- rowSums(slot_matrix)
difference <- recalculated - author_reported
if (any(difference != 0L)) {
  failed <- which(difference != 0L)
  stop("Author score mismatch: ", paste(author$target[failed], collapse = ", "))
}

# Keep the reported number and the slot arithmetic in separate output columns.
comparison <- data.frame(target = author$target, pwas = author$pwas,
                         author_reported_tier_score = author$author_tier_score,
                         author_reported_total = author_reported,
                         recalculated_slot_total = recalculated, difference = difference)
write_tsv(comparison, file.path(output, "author_score_comparison.tsv"))

direction_class <- function(value, reference, slot) {
  value <- trimws(as.character(value))
  if (is.na(value) || !nzchar(value)) return("Not reported (NA)")
  if (slot %in% c("coloc", "ppi")) {
    return(if (grepl("√", value, fixed = TRUE)) "Scored support" else "Not scored")
  }
  if (slot %in% c("bulk", "single_cell")) {
    expected <- if (reference == "+") "↑" else "↓"
    opposite <- if (reference == "+") "↓" else "↑"
  } else {
    expected <- if (reference == "+") "+" else "−"
    opposite <- if (reference == "+") "−" else "+"
  }
  has_expected <- grepl(expected, value, fixed = TRUE)
  has_opposite <- grepl(opposite, value, fixed = TRUE)
  if (has_expected && has_opposite) return("Mixed directions (scored by rule)")
  if (has_expected) return("Scored support")
  if (has_opposite) return("Direction conflict")
  "Not scored"
}

slot_rows <- do.call(rbind, lapply(seq_len(nrow(author)), function(i) {
  data.frame(target = author$target[i], pwas = author$pwas[i], slot = score_slots,
             raw_flag = vapply(score_slots, function(slot) as.character(author[[slot]][i]), character(1)),
             contribution = as.integer(slot_matrix[i, ]),
             display_class = vapply(score_slots, function(slot) {
               direction_class(author[[slot]][i], author$pwas[i], slot)
             }, character(1)), stringsAsFactors = FALSE)
}))
write_tsv(slot_rows, file.path(output, "author_score_slots.tsv"))
slot_rows$target <- factor(slot_rows$target, levels = rev(author$target))
slot_rows$slot <- factor(slot_rows$slot, levels = score_slots)
slot_rows$display_class <- factor(slot_rows$display_class, levels = c(
  "Scored support", "Direction conflict", "Mixed directions (scored by rule)",
  "Not scored", "Not reported (NA)"))
author_plot <- ggplot(slot_rows, aes(x = slot, y = target, fill = display_class)) +
  geom_tile(colour = "white", linewidth = .35) +
  geom_text(aes(label = ifelse(display_class == "Not reported (NA)", "NA",
                               ifelse(display_class == "Direction conflict", "!",
                                      ifelse(display_class == "Mixed directions (scored by rule)", "±",
                                             ifelse(contribution == 1L, "1", "0"))))),
            size = 2.6, colour = "#243746") +
  scale_fill_manual(values = c("Scored support" = "#4C9A83",
    "Direction conflict" = "#D87955", "Mixed directions (scored by rule)" = "#8873A6",
    "Not scored" = "#D9E2E7", "Not reported (NA)" = "#F7F8F8"),
    drop = FALSE, name = NULL) +
  labs(title = "Author-reported evidence flags and slot contributions",
       subtitle = "Colors classify the flags; the table retains raw values. Mixed signs can score under the source rule",
       x = NULL, y = NULL) + theme_course() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(),
        legend.position = "bottom")
save_plot(author_plot, output, "author_score_heatmap", width = 9, height = 8)

# The evidence matrix uses the project's recorded availability/support categories.
candidate_order <- c("IL12B", "SGSH", "TNFAIP3", "DDX58", "IL7R", "PRSS53", "UBLCP1")
primary <- read_tsv("data/evidence/primary_candidates_input.tsv")
evidence <- read_tsv("data/evidence/candidate_evidence_input.tsv")
context <- read_tsv("data/evidence/candidate_context_input.tsv")
stopifnot(nrow(primary) == 7L, setequal(primary$gene, candidate_order), nrow(evidence) == 63L)
stopifnot(all(context$UKB_status == "omitted_by_user_not_replaced_by_summary_GWAS"))
stopifnot(all(is.na(context$UKB_PGS_psoriasis)), all(is.na(context$UKB_PGS_PsA)),
          all(is.na(context$reconstructed_total_score)))
primary$gene <- factor(primary$gene, levels = candidate_order)
evidence$gene <- factor(evidence$gene, levels = rev(candidate_order))
context$gene <- factor(context$gene, levels = candidate_order)
write_tsv(evidence, file.path(output, "candidate_evidence.tsv"))
candidate_context <- merge(primary, context, by = "gene", all.x = TRUE, sort = FALSE)
candidate_context <- candidate_context[match(candidate_order, candidate_context$gene), ]
write_tsv(candidate_context, file.path(output, "candidate_context.tsv"))

evidence$status <- factor(evidence$status, levels = c("supported", "not_supported_in_tested_data",
  "partial_no_observed_support", "not_evaluable"))
component_labels <- c("Bulk" = "Bulk", "Cell Wilcoxon" = "Cells",
  "Pseudobulk" = "Pseudobulk",
  "GTEx_Skin_Not_Sun_Exposed_Suprapubic" = "Skin: unexposed",
  "GTEx_Skin_Sun_Exposed_Lower_leg" = "Skin: exposed",
  "GTEx_Whole_Blood" = "Blood: GTEx", "eQTLGen_Blood" = "Blood: eQTLGen",
  "coloc_mRNA" = "mRNA", "coloc_protein" = "Protein")
component_order <- names(component_labels)[names(component_labels) %in% unique(as.character(evidence$component))]
evidence$component <- factor(evidence$component, levels = component_order)
candidate_plot <- ggplot(evidence, aes(x = component, y = gene, fill = status)) +
  geom_tile(colour = "white", linewidth = .45) +
  geom_text(data = evidence[evidence$status == "not_evaluable", ], label = "×",
            colour = "#77848D", size = 3) +
  facet_grid(. ~ domain, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("supported" = "#4C9A83",
    "not_supported_in_tested_data" = "#D9E2E7",
    "partial_no_observed_support" = "#B8A9CC", "not_evaluable" = "#F7F8F8"),
    labels = c("Observed support", "No support in tested data",
               "Partial coverage; no observed support", "Not evaluable (NA)"),
    drop = FALSE, name = NULL) +
  scale_x_discrete(labels = component_labels) +
  labs(title = "Evidence coverage for the seven reselected candidates",
       subtitle = "Recorded evidence statuses; NA remains distinct from tested absence",
       x = NULL, y = NULL) + theme_course() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        panel.grid = element_blank(), legend.position = "bottom")
save_plot(candidate_plot, output, "candidate_evidence_heatmap", width = 12, height = 5.2)

cat(sprintf("author_rows=%d; slots=%d; arithmetic_matches=%d; candidate_genes=%d; evidence_rows=%d; UKB_scores_missing=%d\n",
  nrow(author), length(score_slots), sum(difference == 0L), length(candidate_order), nrow(evidence),
  sum(is.na(context$UKB_PGS_psoriasis) & is.na(context$UKB_PGS_PsA))))
