# 11 | Genetic associations across related disease outcomes

This public teaching panel contains 22 indexed MVP European-ancestry phenotypes from GWAS Catalog. It differs from the research FinnGen panel; myocardial infarction, basal cell carcinoma, and hypothyroidism are additional examples, and broad phenotypes retain their actual names. The input table records exact definitions, accessions, sample sizes, and differences. These outcomes share a cohort and can overlap with the MVP discovery phenotype. They are not independent replications.

## Inputs and questions

This lesson uses the 17 proteins in `family_17.tsv` as a fixed analysis family, along with prespecified exposure data, within-protein LD, and 22 MVP EUR outcomes. The seven formal candidates come from `primary_candidates.tsv`. Outcome codes, definitions, case counts, and control counts are in `endpoints.tsv`.

The questions are: after instruments are frozen, what are the effect directions and intervals for the 17 proteins across the 22 outcomes? How do the seven formal candidates fare after adjustment across the full set of 374 tests? If we look only at the 154 candidate tests, how does exploratory BH adjustment differ?

## Methods

First, `freeze_instruments()` freezes instruments using the within-protein LD rule. Then `analyse_outcome()` is called for each outcome. This function returns allele-harmonization results and generalized IVW estimates. Each protein–outcome combination has one row; when no effect can be observed, the row is retained with NA.

The 17 × 22 design defines 374 prespecified tests. BH adjustment is first calculated for the P values across all 374 rows. The seven formal candidates are then selected, producing 154 rows that retain their original full-family BH values. A separate exploratory BH adjustment is also calculated for those 154 candidate tests. The full-family adjustment cannot be recreated from 154 P values by setting `n` to 374, because the other computable tests and the missing-value pattern and their ranks are missing.

`endpoints.tsv` retains the actual MVP PheCode definitions. Localized-and-systemic lupus includes both groups and should not be relabelled SLE. Phenotype names and case selection together determine what a result describes.

## Key code

```r
plan <- freeze_instruments(exposures, ld, max_abs_r = .95)
fit <- analyse_outcome(plan, y, ld)

# 先对完整 374 项校正
all_results$BH_original_family_374 <-
  p.adjust(all_results$p, method = "BH", n = 374)

# 再取 7 候选，并单独保留探索性 154 项校正
candidate_results <- all_results[all_results$uniprot %in% formal_ids, ]
candidate_results$BH_candidate_family_exploratory <-
  p.adjust(candidate_results$p, method = "BH", n = 154)
```

## Interpreting the results

`clinical_all_374.tsv` contains the complete analysis family of 374 protein–outcome combinations. Of these, 308 have an estimable P value and 66 are retained as NA because no instrument could be harmonized. In the full family, 37 results have BH<0.05. `clinical_candidates_154.tsv` contains 154 combinations for the seven formal candidates, all with estimable results. Of these, 21 have BH<0.05 under the full 374-test family, and 22 have BH<0.05 under exploratory adjustment within the 154 candidates. In the primary heatmap, color represents the log odds ratio; filled circles mark results passing the 374-test BH threshold of 0.05, and open circles mark results passing only the exploratory candidate-set BH threshold of 0.05. The point shapes distinguish adjustment families; they do not indicate independent replication.

For each row, inspect `phenocode`, `phenotype`, effect direction, confidence interval, and BH columns together. Different disease phenotypes are not the same clinical endpoint, and statistical association alone does not establish causality or clinical utility.

## Exercises and reference answers

### Exercise 1: Multiple testing

Why should we not select 154 rows first and then use `p.adjust(p, "BH", n = 374)` to obtain the primary adjusted values for the seven candidates?

**Answer:** BH ranks depend on the P-value ordering in the complete testing family. Setting `n` to 374 supplies only the family size; it does not restore the other computable tests and the missing-value pattern or their ranks. Calculate BH across all 374 rows first, then retain the results for the seven candidates.

### Exercise 2: Check phenotype definitions

Why should localized-and-systemic lupus not be relabelled SLE?

**Answer:** This phenotype includes localized lupus. Retain the source name and PheCode instead of assuming equivalence to an SLE-only phenotype.

### Exercise 3: Read the results

Choose a formal candidate and compare its `BH_original_family_374` and `BH_candidate_family_exploratory` values in `clinical_candidates_154.tsv`. If the conclusions differ, how should you report them?

**Answer:** Use BH across the complete family of 374 tests as the primary result, and clearly label adjustment within the 154-candidate set as exploratory sensitivity analysis. Do not describe the smaller candidate-set adjustment denominator as new confirmatory evidence.

