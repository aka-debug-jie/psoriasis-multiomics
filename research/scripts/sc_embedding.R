# Seurat pipeline on every cell retained by the paper's explicit QC thresholds.
sc_embedding <- function(root, config, output) {
  suppressPackageStartupMessages(library(Seurat))
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add=TRUE)
  future::plan(future::sequential)
  s <- config$sc
  approved <- jsonlite::fromJSON(file.path(root, config$approved_parameters))$single_cell
  for (key in names(approved)) {
    if (!isTRUE(all.equal(s[[key]], approved[[key]]))) stop("Approved single-cell parameter changed: ", key)
  }
  stopifnot(!s$batch_correction)
  counts <- Matrix::readMM(gzfile(file.path(root, s$counts)))
  counts <- as(counts, "CsparseMatrix")
  features <- read.delim(file.path(root, s$features), stringsAsFactors=FALSE, check.names=FALSE)
  metadata <- read.delim(file.path(root, s$cells), stringsAsFactors=FALSE, check.names=FALSE)
  stopifnot(nrow(counts) == nrow(features), ncol(counts) == nrow(metadata),
            !anyDuplicated(features$gene_id), !anyDuplicated(metadata$cell_id),
            all(Matrix::colSums(counts) == metadata$n_counts),
            all(Matrix::colSums(counts > 0) == metadata$n_features))
  rownames(counts) <- features$gene_id
  colnames(counts) <- metadata$cell_id
  rownames(metadata) <- metadata$cell_id
  object <- CreateSeuratObject(counts, meta.data=metadata, min.cells=0, min.features=0)
  rm(counts); gc()
  object <- NormalizeData(object, normalization.method=s$normalization, scale.factor=s$scale_factor, verbose=FALSE)
  object <- FindVariableFeatures(object, selection.method=s$variable_features_method, nfeatures=s$nfeatures, verbose=FALSE)
  object <- ScaleData(object, features=VariableFeatures(object), verbose=FALSE)
  object <- RunPCA(object, features=VariableFeatures(object), npcs=max(s$PCs), seed.use=s$seed, verbose=FALSE)
  object <- FindNeighbors(object, dims=s$PCs, k.param=s$k_param, verbose=FALSE)
  object <- FindClusters(object, algorithm=s$algorithm, resolution=s$resolution, random.seed=s$seed, verbose=FALSE)
  object <- RunUMAP(object, dims=s$PCs, umap.method=s$umap$implementation,
                    n.neighbors=s$umap$n_neighbors, min.dist=s$umap$min_dist,
                    metric=s$umap$metric, seed.use=s$seed, uwot.sgd=FALSE, verbose=FALSE)
  # Save a checkpoint before annotation; gene and cell identifiers stay unchanged.
  if (s$save_checkpoint) saveRDS(object, file.path(output, "seurat_checkpoint.rds"), compress=FALSE)
  write_tsv <- function(d, name) write.table(d, file.path(output, name), sep="\t", row.names=FALSE, quote=FALSE, na="NA")
  coords <- Embeddings(object, "umap")
  cells <- cbind(object[[]], cluster=as.character(Idents(object)), UMAP_1=coords[,1], UMAP_2=coords[,2])
  write_tsv(cells, "cell_embeddings.tsv")
  pca <- data.frame(cell_id=colnames(object), Embeddings(object, "pca"), check.names=FALSE)
  write_tsv(pca, "cell_pca.tsv")
  write_tsv(data.frame(gene_id=VariableFeatures(object)), "variable_features.tsv")
  group_counts <- as.data.frame(table(cluster=cells$cluster, donor_id=cells$donor_id, tissue_type=cells$tissue_type))
  write_tsv(group_counts, "cluster_donor_counts.tsv")
  panels <- stack(s$marker_panels); names(panels) <- c("symbol","proposed_cell_type")
  panels$gene_id <- vapply(panels$symbol, function(symbol) {
    hit <- features$gene_id[!is.na(features$symbol) & features$symbol == symbol]
    if (length(hit) == 1) hit else NA_character_
  }, character(1))
  write_tsv(panels, "marker_mapping.tsv")
  normalized <- LayerData(object, assay="RNA", layer="data")
  marker_results <- do.call(rbind, lapply(sort(unique(cells$cluster)), function(cluster) {
    idx <- cells$cluster == cluster
    d <- panels
    d$cluster <- cluster
    d$mean_log_expression <- NA_real_
    d$fraction_detected <- NA_real_
    mapped <- which(!is.na(d$gene_id))
    values <- normalized[d$gene_id[mapped], idx, drop=FALSE]
    d$mean_log_expression[mapped] <- Matrix::rowMeans(values)
    d$fraction_detected[mapped] <- Matrix::rowMeans(values > 0)
    d
  }))
  write_tsv(marker_results, "cluster_marker_expression.tsv")
  pdf(file.path(output, "umap_clusters.pdf"), width=9, height=7)
  print(DimPlot(object, reduction="umap", group.by="seurat_clusters", label=TRUE, raster=TRUE))
  dev.off()
  pdf(file.path(output, "umap_samples.pdf"), width=12, height=7)
  print(DimPlot(object, reduction="umap", group.by=c("tissue_type", "donor_id"), raster=TRUE))
  dev.off()
  observed <- panels[!is.na(panels$gene_id), ]
  pdf(file.path(output, "annotation_marker_dotplot.pdf"), width=12, height=8)
  print(DotPlot(object, features=unname(observed$gene_id)) +
          ggplot2::scale_x_discrete(labels=setNames(observed$symbol, observed$gene_id)) + RotatedAxis())
  dev.off()
  generated <- c("cell_embeddings.tsv", "cell_pca.tsv", "variable_features.tsv", "cluster_donor_counts.tsv",
                 "marker_mapping.tsv", "cluster_marker_expression.tsv", "umap_clusters.pdf", "umap_samples.pdf", "annotation_marker_dotplot.pdf")
  if (s$save_checkpoint) generated <- c(generated, "seurat_checkpoint.rds")
  list(status="embedding_completed_annotation_pending_review", cells=ncol(object), genes=nrow(object),
       clusters=length(unique(cells$cluster)), donors=length(unique(cells$donor_id)),
       Seurat=as.character(packageVersion("Seurat")), future_workers=future::nbrOfWorkers(), settings=s,
       generated_files=generated, formal_candidate_list_used=FALSE)
}
