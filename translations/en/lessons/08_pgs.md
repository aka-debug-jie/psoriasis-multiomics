# Lesson 08: Build an IL7R Protein Score from Two Genetic Loci

## Questions and inputs

How can genetic variant effect weights be combined with individual genotypes to calculate a score for a target protein? This lesson uses 503 individuals of European ancestry from 1000 Genomes Phase 3 and IL7R protein cis loci from the same UK Biobank Pharma Proteomics Project (UKB-PPP) source study. Variants are selected using P<5×10⁻⁸, clear allele direction, and presence in the genotype reference; two variants are then retained under r²<0.01. These weights correspond to IL7R protein, not psoriasis effects.

## Core steps

First harmonize the effect-allele direction, then extract 0/1/2 dosages for the two variants. Each person's score is the sum across variants of “effect-allele dosage × protein effect size,” then standardized to inspect its distribution. Individuals missing a dosage at either variant are excluded from the primary score; PLINK results are used only to check the weighted sum calculated in R. In the current input summary, all 503 people have complete dosages at both variants, and the observed r² is about 0.0075.

A separate set of 5,000 genome-wide LD-pruned SNPs is used to calculate genetic principal components as a background check of population structure. These serve a different purpose from the two score variants. The score measures the genetically predicted IL7R protein signal carried by those two variants; it is not a clinical measurement or an individual's probability of psoriasis. It must not be joined by row to the lung cancer survival lesson data: the two datasets contain different populations, individuals, and outcomes.

## How to read the plots

The score histogram shows the distribution of this small score among 503 samples of European ancestry; the x-axis is the standardized score and the y-axis is the number of individuals. It does not show disease risk. The principal-component plot shows PC1 and PC2 from 5,000 genome-wide variants, with color indicating the reference population. It is used to inspect genetic structure and does not represent score performance. If the plot or score table has not yet been generated in the current run, rerun the script before reading it; do not infer plot results from the method description.

## Exercises and answers

1. **What might happen if the selection threshold is changed to P<1×10⁻⁶?** The selected variant set may change. Recheck effect alleles, reference variants, LD, and the number of people with scores; do not simply relax the threshold and interpret the result.
2. **Why doesn't the primary analysis automatically impute the mean when one variant is missing?** This lesson is fixed to use individuals with complete dosages at both variants, avoiding treating missing individuals as having an average genetic dosage. If missing-data handling is changed, state it clearly and perform a sensitivity analysis.
3. **Does a higher score mean a higher chance of psoriasis?** No. The weights come from IL7R protein associations; the score is not a disease-risk model and has not been calibrated against disease outcomes here.

Run `source("lessons/08_pgs.R")`. See the [script](08_pgs.R), [score distribution](../results/08_pgs/protein_score.png), and [genetic principal components](../results/08_pgs/ancestry_pcs.png). For these 503 people, both score variants are complete, and the maximum absolute difference between the R and PLINK sums is 2.22×10⁻¹⁶.
