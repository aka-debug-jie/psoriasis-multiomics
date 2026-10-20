# 03 | Proteome-wide MR and external outcomes

## Research question and inputs

This lesson screens a fixed family of a fixed testing family of 2,709 proteins for associations with psoriasis, then checks the direction in the Stuart external outcome. Inputs include exposure data, within-protein LD matrices, and outcome summary statistics from MVP and Stuart. Exposure rows are deduplicated by protein and variant, then instruments are frozen by F statistic; the `|r|≤0.95` rule is used to retain instruments. After allele harmonization, the analysis uses generalized IVW with signed-LD covariance.

Each result represents a protein, not an SNP: variants retained for the same protein are combined into one generalized IVW estimate. `frozen_variants` and `observed_variants` give the number of frozen instruments and the number successfully harmonized in the outcome, respectively. `Q` and `se_mre` record heterogeneity and dispersion-adjustment information.

On Windows, open `teaching.Rproj` and run `source('lessons/03_proteome_mr.R')` in the RStudio Console. See the [R script](03_proteome_mr.R).

## Screening and results

The discovery criterion is MVP P<0.05/2,709; external nominal support then requires P<0.05 in Stuart. The results table contains all 2,709 proteins, and seven candidates pass both stages: DDX58, IL7R, TNFAIP3, IL12B, SGSH, PRSS53, and UBLCP1. All seven have the same direction in the externally nominally significant outcome(s) as in MVP. For example, IL12B has an MVP beta of −0.1726 and P=1.27×10⁻¹⁴. MVP standard errors were inferred from reported intervals, and the source P-value convention has not been confirmed, so the source `p_reported` may differ from the normal-test P value implied by beta/SE.

`candidates.tsv` includes only proteins that meet both discovery and external nominal-support criteria. The other 2,702 proteins remain in the full results table, allowing estimates near the thresholds to be checked.

See the [cross-outcome effect plot](../results/03_proteome_mr/discovery_and_validation.png), [all-proteins table](../results/03_proteome_mr/all_proteins.tsv), and [candidate table](../results/03_proteome_mr/candidates.tsv).

## Exercises and answers

1. Why is the denominator for the discovery threshold 2,709? **Answer:** This is the protein-family size fixed before testing; reducing it based on significance would change the testing family.
2. How could you change the stringency of the discovery criterion? **Answer:** Clarify the study design before changing the numerator in the Bonferroni threshold or the prespecified family size. Raising the P-value threshold would produce more candidates, but it would change the discovery rule.
3. What does NA mean when a protein has no instrument for an outcome? **Answer:** The current data do not provide a harmonizable effect estimate; NA does not mean a zero effect. First check the frozen instruments and harmonization table. Also distinguish `frozen_variants=0` from instruments that were frozen but have `observed_variants=0`.
