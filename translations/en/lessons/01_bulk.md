# 01 | Bulk RNA-seq differential expression

## Inputs

This lesson uses the NCBI-reprocessed integer count matrix and sample metadata for GSE54456. The count matrix contains 171 samples; the metadata retains three samples that are absent from the matrix and marks them in `sample_inclusion.tsv`. Expression values are keyed by GeneID and mapped to current NCBI symbols. The annotation file contains only the GeneID, Symbol, and Synonyms needed for this lesson. The current symbol for GeneID 23586 is RIGI; `DDX58` is an alias.

NCBI's data documentation establishes that these values are raw counts. The integer check in the code helps detect input errors; integer values alone do not prove the data's origin.

The analysis compares lesional psoriatic skin with normal skin. A gene is called differentially expressed when its BH-adjusted P value is below 0.05 and its absolute log2 fold change is at least 1.

## Questions

1. How can the data source and value checks help distinguish raw counts from RPKM/TPM?
2. What problems do TMM and voom address, respectively?
3. Why is the BH correction based on all tested genes rather than only candidate genes?
4. Do the direction, effect size, and statistical significance agree for these seven candidate genes?

## Key code

```r
# Check counts and remove all-zero genes
stopifnot(all(counts >= 0), all(counts == floor(counts)))
keep <- rowSums(counts) > 0
y <- edgeR::DGEList(counts[keep, , drop = FALSE], group = group)

# TMM, voom, linear modeling, and empirical Bayes testing
y <- edgeR::calcNormFactors(y, method = "TMM")
design <- model.matrix(~ group)
v <- limma::voom(y, design)
fit <- limma::eBayes(limma::lmFit(v, design))
result <- limma::topTable(fit, coef = "grouppsoriasis", number = Inf,
                          adjust.method = "BH", sort.by = "none")
```

## Interpreting the results

The analysis tests 37,065 non-zero genes, with the contrast defined as psoriasis − normal. A positive `logFC` means higher expression in lesional skin; a negative value means lower expression. `adj.P.Val` is the BH-adjusted result across all tested genes. The volcano plot marks genes using both adjusted P value and effect size; the candidate-effect plot shows the seven candidates' log2 fold changes and empirical-Bayes moderated 95% intervals.

The candidates are IL12B, SGSH, TNFAIP3, DDX58, IL7R, PRSS53, and UBLCP1. DDX58 is mapped through the NCBI alias for GeneID 23586. The sample QC table reports library size and TMM factor for each sample. The sample QC plot shows the distribution of library sizes; it does not replace checks of sample identity or batch.

This lesson uses NCBI-reprocessed counts from the same study, not the authors' original historical count matrix. Expression differences describe group-level associations and, by themselves, do not establish causality or clinical effects.

## Exercises and answers

### Exercise 1: Concept check

If a gene has a count of zero in every sample, does it enter the model? How does all-zero filtering differ from low-expression filtering?

**Answer:** It does not enter the model. All-zero filtering removes only genes with no observed counts; low-expression filtering requires an additional threshold. This lesson does not add a low-expression filter.

### Exercise 2: Reading the code

Why is `normal` the reference group in `model.matrix(~ group)`? What does the `grouppsoriasis` coefficient represent?

**Answer:** The factor levels are ordered as normal, then psoriasis, so normal is the reference. `grouppsoriasis` represents the log2 expression change in lesional skin relative to normal skin.

### Exercise 3: Interpreting results

In `candidate_effects.tsv`, find candidates with `adj.P.Val < 0.05` and `abs(logFC) >= 1`, and state the coefficient direction.

**Answer:** IL12B, DDX58, IL7R, and PRSS53 meet both thresholds and have positive `logFC`, indicating higher expression in lesional psoriatic skin. The other candidates do not meet the combined criterion; passing the thresholds is not evidence of mechanism or causality.

Run `source("lessons/01_bulk.R")` and compare the [full script](01_bulk.R), [volcano plot](../results/01_bulk/bulk_volcano.png), and [candidate-effect plot](../results/01_bulk/candidate_effects.png).
