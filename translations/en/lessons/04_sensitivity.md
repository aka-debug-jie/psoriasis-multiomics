# 04 | Sensitivity to instrument selection, direction, and model

## Research question and inputs

This lesson examines whether genetic effect estimates change with instrument selection, assumed disease prevalence, or the MR method. Before this lesson, 30 MHC proteins were excluded from the 2,709-protein family based on the coding gene's location. This lesson separately flags whether an instrument SNP lies in the 25–34 Mb region; gene location and instrument location are distinct checks. The script calculates F statistics for 3,802 exposure records and flags the MHC region, then compares the primary analysis with a stricter LD rule. The multi-instrument section uses a fixed input table already clumped at P<5×10⁻⁶, r²<0.1, and 10 Mb; it does not select instruments again here. A separate PGS analysis uses P<5×10⁻⁸ and r²<0.01, so those thresholds should not be mixed with the MR settings.

On Windows, open `teaching.Rproj` and run `source('lessons/04_sensitivity.R')` in the RStudio Console. See the [R script](04_sensitivity.R).

## Analyses and results

The primary LD branch uses |r|≤0.95; the strict branch uses |r|≤√0.1. Each branch has results for 2,709 proteins, with 193 rows at P<0.05 in both. Beta could be estimated for 2,096 and 2,091 rows, respectively. The minimum input F is 21.54 and the median is 479.61; 4 of 3,802 rows lie in the MHC region flagged by the script. The Steiger grid has 72 scenario combinations, and `forward` is TRUE for K=1%, 2%, 3%, and 5% under both R² algorithms. Of 88 rows in the multi-instrument methods table, 55 have estimates and 33 were not run because they had fewer than three instruments, as required for the median, mode, or Egger methods.

`LD_and_dispersion.tsv` stores primary and strict LD results and heterogeneity-adjustment columns for each protein. The same number of nominal P<0.05 rows does not mean every protein has the same effect estimate. A `forward` flag in the Steiger table applies only to the listed prevalence and R² algorithm. Multi-instrument comparisons include IVW, median, mode, and Egger; single-instrument rows use the Wald ratio in the script.

See the [method-comparison plot](../results/04_sensitivity/method_comparison.png), [LD/dispersion table](../results/04_sensitivity/LD_and_dispersion.tsv), [Steiger grid](../results/04_sensitivity/steiger_grid.tsv), and [multi-instrument methods table](../results/04_sensitivity/multi_instrument_methods.tsv).

## Exercises and answers

1. Why does the strict LD branch use √0.1? **Answer:** The helper takes an upper bound for |r|, while the design threshold is r²=0.1, so the input is √0.1.
2. How can you change the Steiger assumptions? **Answer:** Change the K scenarios or R² method and compare the combinations again. The conclusion applies only to the listed input assumptions; all combinations being forward does not establish the result for other prevalences.
3. Why do some methods require at least three instruments? **Answer:** Under this workflow, the weighted median, mode, and Egger methods cannot be estimated with too few instruments. The status column records this limitation; do not replace missing results with zero. Use `cohort`, `gene`, and `mode` to inspect which input groups can be analyzed.
