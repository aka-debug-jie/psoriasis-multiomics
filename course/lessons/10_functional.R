# 第 10 课：过度代表分析、蛋白互作网络与药物关联
# 工作目录设为 teaching/ 后运行：Rscript lessons/10_functional.R

source("R/common.R")
required <- c("ggplot2", "igraph")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("请先安装课程依赖：", paste(missing, collapse = ", "))

out <- lesson_output("10_functional")
data_dir <- "data/functional"
input <- function(name) file.path(data_dir, name)
candidates <- read_tsv(input("candidates.tsv"))
background <- read_tsv(input("opportunity_background.tsv"))
term_gene <- read_tsv(input("term_gene.tsv"))
annotated_input <- read_tsv(input("annotated_background.tsv"))

# 1. 先定义研究问题中的候选集与机会背景，再看哪些基因有本体注释。
query <- unique(as.character(candidates$GeneID))
background_genes <- unique(as.character(background$GeneID))
stopifnot(length(query) == 7L, length(background_genes) == 2043L,
          all(query %in% background_genes),
          all(term_gene$family %in% c("BP", "MF", "CC", "Reactome")))

# 2. 对每个家族分别用超几何分布计算右尾概率，并在家族内部校正 BH。
# N=该家族有注释的背景基因数；K=term 在注释背景中的基因数；
# n=候选集中有该家族任一注释的基因数；k=候选与 term 的交集。
families <- c("BP", "MF", "CC", "Reactome")
results <- vector("list", length(families))
names(results) <- families
for (family in families) {
  rows <- term_gene[term_gene$family == family, , drop = FALSE]
  rows$GeneID <- as.character(rows$GeneID)
  term_ids <- unique(as.character(rows$term))
  annotated <- as.character(annotated_input$GeneID[annotated_input$family == family])
  query_annotated <- intersect(query, annotated)
  n <- length(query_annotated)
  N <- length(annotated)
  one <- lapply(term_ids, function(term_id) {
    term_rows <- rows[rows$term == term_id, , drop = FALSE]
    members <- unique(term_rows$GeneID)
    hits <- sort(intersect(members, query_annotated))
    K <- length(members)
    k <- length(hits)
    data.frame(
      family = family,
      term = term_id,
      description = as.character(term_rows$description[1]),
      overlap = k,
      query_annotated = n,
      term_background = K,
      background_annotated = N,
      p_value = stats::phyper(k - 1L, K, N - K, n, lower.tail = FALSE),
      gene_ids = paste(hits, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, one)
  table$fdr_bh <- stats::p.adjust(table$p_value, method = "BH")
  table <- table[order(table$fdr_bh, table$p_value, table$term), , drop = FALSE]
  rownames(table) <- NULL
  results[[family]] <- table
  write_tsv(table, file.path(out, paste0("enrichment_", family, ".tsv")))
}
enrichment <- do.call(rbind, results)
write_tsv(enrichment, file.path(out, "enrichment_all_families.tsv"))

# 3. 每个家族单独统计检验项与显著项；保留零交集 term。
comparison <- do.call(rbind,lapply(families,function(f) data.frame(family=f,
  tested_terms=nrow(results[[f]]),significant=sum(results[[f]]$fdr_bh<.05))))
write_tsv(comparison,file.path(out,'family_summary.tsv'))
# 用 q = k - 1 计算右尾 P(X >= k)，零交集的 P 为 1。
show <- do.call(rbind,lapply(results,function(z) head(z[z$overlap>0,],6)))
show$label <- factor(paste(show$family,show$description,sep=': '),levels=rev(paste(show$family,show$description,sep=': ')))
p <- ggplot2::ggplot(show,ggplot2::aes(-log10(fdr_bh),label,size=overlap,colour=family))+
  ggplot2::geom_vline(xintercept=-log10(.05),linetype=2,colour='grey65')+
  ggplot2::geom_point()+ggplot2::labs(x='-log10(BH FDR)',y=NULL,title='Functional enrichment',size='Hits',colour='Family')+theme_course()
save_plot(p,out,'functional_enrichment',11,8)

# 4. Draw the finalized STRING induced network. Edges are functional associations;
# they do not establish direct binding or a drug effect.
candidate_names <- as.character(candidates$gene_symbol)
network <- read_tsv(input("string_induced_edges.tsv"))
stopifnot(nrow(network) == 77L, all(network$ncbiTaxonId == 9606))
graph <- igraph::graph_from_data_frame(
  data.frame(from = network$preferredName_A, to = network$preferredName_B,
             score = as.numeric(network$score)), directed = FALSE
)
set.seed(20261017)
layout <- igraph::layout_with_fr(graph, weights = igraph::E(graph)$score)
layout <- igraph::norm_coords(layout,xmin=-1,xmax=1,ymin=-1,ymax=1)
vertices <- data.frame(name = igraph::V(graph)$name, x = layout[, 1], y = layout[, 2],
                       candidate = igraph::V(graph)$name %in% candidate_names)
edgelist <- as.data.frame(igraph::as_edgelist(graph), stringsAsFactors = FALSE)
names(edgelist) <- c("from", "to")
edges <- merge(edgelist, vertices[c("name", "x", "y")], by.x = "from", by.y = "name")
names(edges)[names(edges) == "x"] <- "x"
names(edges)[names(edges) == "y"] <- "y"
edges <- merge(edges, vertices[c("name", "x", "y")], by.x = "to", by.y = "name",
               suffixes = c("", "_end"))
names(edges)[names(edges) == "x_end"] <- "xend"
names(edges)[names(edges) == "y_end"] <- "yend"
# Restore edge scores from input without relying on vertex ordering.
score_key <- paste(network$preferredName_A, network$preferredName_B, sep = "|")
reverse_key <- paste(network$preferredName_B, network$preferredName_A, sep = "|")
edge_key <- paste(edges$from, edges$to, sep = "|")
score_lookup <- c(stats::setNames(as.numeric(network$score), score_key),
                  stats::setNames(as.numeric(network$score), reverse_key))
edges$score <- unname(score_lookup[edge_key])
network_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = edges,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = score),
    colour = "#9AA6B2", alpha = .48) +
  ggplot2::geom_point(data = vertices,
    ggplot2::aes(x = x, y = y, fill = candidate), shape = 21, size = 3.4,
    colour = "white", stroke = .7) +
  ggplot2::geom_text(data = vertices[vertices$candidate, , drop = FALSE],
    ggplot2::aes(x = x, y = y, label = name), nudge_y = .07, size = 3,
    fontface = "bold", check_overlap = TRUE) +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "#A9B4BE", `TRUE` = "#D96B32"),
                             labels = c("Network neighbor", "Candidate protein"), name = NULL) +
  ggplot2::scale_linewidth(range = c(.25, 1.2), guide = "none") +
  ggplot2::labs(title = "STRING induced functional association network",
                subtitle = "Candidate proteins are highlighted; edges are not direct drug effects",
                x = NULL, y = NULL) +
  ggplot2::coord_equal(ylim=c(-1.15,1.15)) + theme_course() + ggplot2::theme_void(base_size=12) +
  ggplot2::theme(legend.position="bottom",plot.title=ggplot2::element_text(face="bold"),plot.margin=ggplot2::margin(12,12,12,12))
save_plot(network_plot, out, "string_induced_network", width = 8.5, height = 7)

# 5. Keep drug evidence classes separate: direct target, complex/family component,
# and candidate-to-drug-target STRING association.
direct <- read_tsv(input("drug_direct_targets.tsv"))
components <- read_tsv(input("drug_complex_components.tsv"))
links <- read_tsv(input("drug_network_associations.tsv"))
stopifnot(nrow(direct) == 0L, nrow(components) == 5L, nrow(links) == 53L)
class_counts <- list(
  direct = data.frame(uniprot = as.character(direct$uniprot), drug = as.character(direct$molecule_chembl_id)),
  component = data.frame(uniprot = as.character(components$uniprot), drug = as.character(components$molecule_chembl_id)),
  network = data.frame(uniprot = as.character(links$candidate_uniprot), drug = as.character(links$molecule_chembl_id))
)
association_rows <- do.call(rbind, lapply(names(class_counts), function(class_name) {
  x <- unique(class_counts[[class_name]])
  data.frame(uniprot = x$uniprot, drug = x$drug, evidence_class = rep(class_name,nrow(x)),
             stringsAsFactors = FALSE)
}))
association_rows$gene_symbol <- candidates$gene_symbol[match(association_rows$uniprot, candidates$uniprot)]
association_rows <- association_rows[!is.na(association_rows$gene_symbol), , drop = FALSE]
write_tsv(association_rows, file.path(out, "drug_association_evidence.tsv"))
counts <- as.data.frame(table(factor(association_rows$gene_symbol, levels = candidates$gene_symbol),
                              factor(association_rows$evidence_class,
                                     levels = c("direct", "component", "network"))),
                        stringsAsFactors = FALSE)
names(counts) <- c("gene_symbol", "evidence_class", "drug_count")
counts$drug_count <- as.integer(counts$drug_count)
write_tsv(counts, file.path(out, "drug_association_counts.tsv"))
drug_plot <- ggplot2::ggplot(counts,
  ggplot2::aes(x = gene_symbol, y = drug_count, fill = evidence_class)) +
  ggplot2::geom_col(position = "dodge", width = .72) +
  ggplot2::scale_fill_manual(values = c(direct = "#B44A4A", component = "#D99B31",
                                        network = "#4D8C78"),
    breaks = c("direct", "component", "network"),
    labels = c("Direct target", "Complex/family component", "STRING-associated target"),
    name = NULL) +
  ggplot2::labs(title = "Curated drug evidence by candidate",
    subtitle = "Unique molecule counts by evidence class; STRING links indicate association only",
    x = NULL, y = "Unique molecules") +
  theme_course() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
save_plot(drug_plot, out, "drug_association_classes", width = 8.5, height = 5.4)

cat("Candidates:", length(query), "of", length(background_genes), "unique GeneIDs from 2045 proteins\n")
print(comparison, row.names = FALSE, digits = 5)
cat("STRING induced edges:", nrow(network), "\n")
cat("Drug evidence rows: direct=", nrow(direct), ", complex/family components=", nrow(components),
    ", network associations=", nrow(links), "\n", sep = "")
cat("Results:", normalizePath(out), "\n")
