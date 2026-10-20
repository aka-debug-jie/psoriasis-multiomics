# 02 | A single genetic instrument: the Wald ratio

## Research question and inputs

This lesson asks whether genetically predicted IL7R protein levels are associated with psoriasis risk. The inputs are an IL7R cis-pQTL exposure record from UKB-PPP and Stuart psoriasis outcome summary statistics. The code selects the variant with the largest F statistic among the IL7R UKB-PPP records, then aligns the outcome effect to the same effect-allele direction.

On Windows, open `teaching.Rproj` and run `source('lessons/02_single_mr.R')` in the RStudio Console. See the [R script](02_single_mr.R).

## Statistical steps

The Wald ratio is the outcome effect divided by the exposure effect. Its standard error uses the first-order NOME formula; a delta-method standard error that includes exposure uncertainty is also calculated for comparison. The 95% interval is calculated on the log-OR scale and then exponentiated to obtain the OR. The vertical line at OR=1 in the plot marks no association.

In the results table, `beta_exposure` is the SNP–protein association, and `beta_outcome_aligned` is the disease effect oriented to the same effect allele. `orientation` is 1 when the directions agree and −1 when the outcome effect must be flipped. `se` is from the first-order formula; `se_delta` includes uncertainty in the exposure estimate, so their interval widths can be compared.

The table's P value comes from a two-sided normal test of the MR beta and its standard error. When viewing the OR, also inspect the exponentiated confidence interval; a beta interval is not an OR interval.

## Reading the results

The selected variant is rs11742270, with instrument F=16,788.6. The results table reports beta=0.0667, SE=0.0215, OR=1.0689 (95% CI 1.0249–1.1149), and P=1.91×10⁻³. An OR above 1 means that each analysis unit increase in genetically predicted IL7R is associated with higher outcome odds. This is a single-instrument estimate. See the [Wald result plot](../results/02_single_mr/wald_ratio.png) and [results table](../results/02_single_mr/wald_ratio.tsv).

## Exercises and answers

1. What happens to the Wald ratio if the aligned outcome beta changes sign? **Answer:** With the exposure beta unchanged, the MR beta changes sign and the OR moves from above 1 to below 1, or vice versa.
2. What would change if instrument selection used the smallest F instead? **Answer:** A different SNP might be selected, so the exposure and outcome effects and the estimate could change. Selection should still be based only on exposure data. A larger F usually indicates a stronger instrument association, but does not replace checking the direction.
3. What should you do if the outcome row is missing after `match()`, or the alleles cannot be aligned? **Answer:** Check whether the rsID is present and whether EA/OA are clear. Do not treat a missing effect as zero or divide without alignment. Also confirm that the outcome beta and standard error use the same effect allele.
