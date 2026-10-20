# 12 | Evidence scoring and coverage heatmaps

This chapter integrates a fixed snapshot of the research reconstruction. The public examples in Chapters 05, 06, and 11 use different inputs and do not automatically overwrite this evidence table. The example PGS and Cox datasets supply no UK Biobank validation evidence. Keep the snapshot source attached to every evidence column.

## Inputs

This lesson uses two groups of inputs. `author_score_input.tsv` retains the authors' 29-row main table, PWAS directions, nine raw evidence flags, and author-reported tier fields. `primary_candidates_input.tsv`, `candidate_evidence_input.tsv`, and `candidate_context_input.tsv` come from the study's final data and record selection information, evidence coverage, and context for seven reselected candidates.

The authors' scores are recalculated slot by slot using nine fields: `sensitivity`, two SMR types, two PGS types, bulk, single-cell, colocalization, and PPI. The total reported by the authors is used only to compare with the sum of the nine slots. The seven reselected candidates lack the UKB PGS entries, so this lesson shows only evidence coverage for them and does not calculate an author-comparable total score.

## Questions

1. How is each of the nine evidence slots scored based on the PWAS direction?
2. How does the heatmap distinguish blank raw flags, opposite directions, mixed directions, and support signals?
3. What is the difference between “not evaluated” and “evaluated but not supported”?
4. Why can't we fill in a total comparable to the authors' score for the seven reselected candidates?

## Key code

```r
if (slot %in% c("coloc", "ppi")) {
  return(as.integer(grepl("√", value, fixed = TRUE)))
}
if (reference == "+") {
  if (slot %in% c("bulk", "single_cell"))
    return(as.integer(grepl("↑", value, fixed = TRUE)))
  return(as.integer(grepl("+", value, fixed = TRUE)))
}
```

The remaining negative-direction checks, author-score parsing, and summation across nine slots follow the same rules in the study script. In the heatmap, “direction conflict” means the raw flag contains only a direction opposite to the PWAS direction. If a raw flag contains both the expected and opposite directions, it is scored according to the original rule and separately classified as “mixed directions.”

## Interpreting the results

For all 29 author rows, the reported score agrees with the recalculated sum across nine slots. The authors' heatmap retains the raw flag for each slot and distinguishes scored support, direction conflict, mixed directions, unscored, and unreported values. Whether a mixed-direction flag receives a score is still determined by the original rule; the plot classification helps with interpretation and does not change the contribution.

The seven reselected candidates have 63 evidence-component records: 21 show observed support, 20 were tested without support, 19 are not evaluable, and 3 have incomplete coverage with no observed support. The coverage heatmap uses the status definitions in the final evidence matrix. Not evaluated remains NA and is not treated as zero or as evidence against a candidate. The context table retains the original direction descriptions for cell-level and pseudobulk analyses.

UKB psoriasis and PsA PGS values are missing for all seven candidates, and the reconstructed total score in the context table is also NA. Therefore, the output presents evidence coverage and original context only; it does not calculate or fill in a total score for the seven candidates that could be compared with the authors' score.

## Exercises and reference answers

### Exercise 1: Read the rule

One candidate has a PWAS direction of `+` and a `bulk` flag of `↑`; another has a PWAS direction of `−` and a `bulk` flag of `↑`. How many points does each receive under this lesson's rule?

**Answer:** The first receives 1, and the second receives 0 and is shown as a direction conflict. Bulk uses up/down arrows; regular genetic-evidence slots use plus/minus signs.

### Exercise 2: Classify results

In the candidate evidence heatmap, what is the difference between `not_evaluable` and `not_supported_in_tested_data`?

**Answer:** `not_evaluable` means there is no observation available for evaluation, so the support value remains NA. `not_supported_in_tested_data` means the item was fully tested but did not meet the support criterion. They should not both be interpreted as negative results.

### Exercise 3: Check completeness

Why can we compare recalculated scores with reported scores for the authors' 29 rows but not calculate the same kind of score for the seven reselected candidates?

**Answer:** The authors' main table provides all nine scoring slots and their reported scores, so the arithmetic can be reproduced and checked slot by slot. The seven new candidates lack the UKB PGS entries, and the study's context table keeps the total score as NA. Missing entries cannot be filled with 0, so a score comparable to the authors' total cannot be claimed.

Run `source("lessons/12_evidence.R")`. See the [scoring script](12_evidence.R) and [candidate evidence heatmap](../results/12_evidence/candidate_evidence_heatmap.png).
