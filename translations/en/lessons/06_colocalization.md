# 06 | Regional Colocalization and Prior Sensitivity

## Research question and inputs

This lesson tests whether the regional signals for two molecular traits and psoriasis are more consistent with a shared causal variant. The traits are IL7R protein and ERAP2 mRNA; for each, we compare windows of index ±500 kb and ±1,000 kb. The inputs are shared, harmonizable SNP summary effects in the corresponding regions, with the outcome direction aligned by allele. The ±1 Mb request may extend beyond the cis coverage of the source study; see the actual boundaries in the [region table](../data/coloc/regions.tsv). The analysis retains all shared variants available in each region.

On Windows, open `teaching.Rproj` and run `source('lessons/06_colocalization.R')` in the RStudio Console. See the [R script](06_colocalization.R) for details.

## Method and results

Each region is analyzed with `coloc.abf`, keeping p1=p2=10⁻⁴ fixed and setting the shared causal prior p12 in turn to 10⁻⁶, 10⁻⁵, and 10⁻⁴. The five posterior probabilities H0–H4 represent: neither trait is associated; only trait one is associated; only trait two is associated; both are associated but have different causal variants; and both share a causal variant. The script labels PP.H4>0.8 as `strong_H4`.

Each result row corresponds to one gene, window, and p12 combination; `n_variants` is the number of shared SNPs supplied for that region. The prior-sensitivity table retains the full posterior for all five hypotheses rather than reporting only H4. The regional plot shows association positions and -log10(P) for both traits in the 500 kb window, to help inspect whether the signal positions overlap.

There are 12 scenario results in total. The number of shared SNPs in the IL7R region is 3,935 (500 kb) and 4,042 (1,000 kb); at p12=10⁻⁵, PP.H4 is 0.8467 and 0.8465, respectively, exceeding the script threshold. The ERAP2 region has 3,313 and 5,736 SNPs; under p12=10⁻⁵, H1 is about 0.6126/0.5615, H3 is about 0.3650/0.4180, and H4 is about 0.0224/0.0205. See the [colocalization summary table](../results/06_colocalization/coloc_summary.tsv) and [regional plot](../results/06_colocalization/regional_associations.png).

## Exercises and answers

1. How should PP.H4=0.8467 be interpreted? **Answer:** Under the current region, inputs, and priors, the posterior probability gives greater support to the shared causal variant model. It is not the frequency of colocalization in this region and does not prove that the molecular trait causes disease. Also inspect H3 and the other prior scenarios.
2. What is affected by changing p12? **Answer:** It changes the prior weight of the shared causal model. Compare the 12 result rows and check whether the conclusion is stable. Do not report only the row with the largest H4.
3. How should you troubleshoot SNPs that cannot be harmonized or have duplicate IDs? **Answer:** Check SNP uniqueness and EA/OA coding. The current script stops when it encounters duplicates or effects that cannot be harmonized; do not simply delete them and continue interpreting. Confirm that the coordinates and both types of effects refer to the same shared variant set.
