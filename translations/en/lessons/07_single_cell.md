# Lesson 07: Cell Types and Candidate Genes in Single-Cell Data

## Questions and inputs

What major cell types are present in the skin single-cell data? Does expression of the 7 candidate genes differ across cell types? The input matrix contains 33,538 genes and 10,708 cells before QC, from 22 samples (8 normal and 14 psoriasis). To limit memory use, up to 500 cells per sample were sampled before QC; this teaching sample is not a full cell atlas.

## Core steps

The script first performs QC using the number of detected genes, total counts, and mitochondrial percentage for each cell; 7,938 cells are retained in this run. It then applies LogNormalize, selects 2,000 highly variable genes, scales only those genes, calculates 30 principal components, clusters cells, and plots UMAP. Mean marker expression and the fraction of expressing cells are summarized in a dot plot. Cell-type labels suggested by the rules should be checked against the dot plot and biological knowledge; cluster IDs from a different full-data analysis must not be carried over.

Cell-level Wilcoxon tests perform two separate comparisons: lesional versus normal cells within the same cell type, and one cell type versus all remaining cells. Each comparison is corrected independently. Both tests examine only the 7 preselected candidate genes, with BH correction within the family of candidate genes × evaluable cell types. Cells are not independent biological replicates, so these results describe only the current sample. The script separately analyzes donor-level pseudobulk: raw counts aggregated from the full dataset by sample and cell type are used with edgeR to compare disease and normal groups. Each sample is a replicate; this answers a different question from the Wilcoxon tests on sampled cells.

## How to read the result plots

Each point in the UMAP represents one cell, and its color indicates the suggested type inferred from marker expression; the two panels are separated by disease status. Neighboring points do not imply a direct relationship between cells, and separation does not prove a disease effect. In the dot plot, point size indicates the fraction of cells in that group expressing the marker, and color indicates mean log expression. First check whether a set of markers jointly supports a label, then inspect individual genes. The results table records 25 clusters and sample counts at each step.

## Exercises and answers

1. **What is expected if the mitochondrial QC threshold is made stricter than 10%?** More cells may be excluded, and cell composition and clustering may also change. Recheck the per-sample QC table; do not compare only the appearance of the UMAP.
2. **Is it sufficient to label a cluster as a cell type based on only one marker?** No. Return to the dot plot and check multiple markers, expression percentages, and means for that type; if the evidence conflicts, the cluster can be left unassigned.
3. **Why can cell-level Wilcoxon tests not be treated as donor-level disease tests?** Cells from the same donor are not independent replicates. Cell-level tests describe this sample; donor-level pseudobulk aggregates cells within samples and uses the sample as the statistical unit.

In RStudio, run `source("lessons/07_single_cell.R")`. See the [script](07_single_cell.R), [UMAP](../results/07_single_cell/cell_atlas.png), and [marker plot](../results/07_single_cell/marker_dotplot.png).
