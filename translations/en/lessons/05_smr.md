# 05 | SMR Association and HEIDI Heterogeneity

The public-data edition uses reprocessed GTEx whole blood from eQTL Catalogue QTD000356 (670 participants; variant-level an/2 is 636–670). ERAP2 and LIPA were chosen as teaching examples using exposure eQTL strength, before examining disease effects. Their inputs and results are separate from the original research eQTLGen analysis.

Is the expression of ERAP2 and LIPA in blood associated with psoriasis? If so, do the regional signals for expression and disease support a shared variant? SMR and HEIDI provide information about these two questions, respectively.

Run the code in order in the [cell-by-cell lesson](../notebooks/05_smr.ipynb). The [script](05_smr.R) uses the complete public cis statistics for two genes in GTEx / eQTL Catalogue blood: 7,131 associations for ERAP2 and 7,199 for LIPA. It also provides the corresponding MVP GWAS and a 1000 Genomes EUR LD subset of 503 people. See the [data documentation](../data/README.md) for sources.

The script first converts readable ESD tables to BESD to demonstrate the tool's input format. The primary analysis uses `public_eqtl`, a fixed-order BESD generated from the same ESD inputs. ESD conversion sorts by genomic position, while HEIDI's LD pruning in tied cases can depend on input order; reordering the same summary statistics may select different variants for testing.

The SMR instrument-selection threshold is eQTL P<5×10⁻⁸; the example association rule is SMR P<0.05 and HEIDI P>0.01. HEIDI also uses eQTL P<0.0015654 and an LD r² range of 0.05–0.9. These three P thresholds serve different purposes. All parameters are listed at the start of the script.

The observed results with Linux SMR 1.4.3 are:

| Gene | SMR beta | SMR P | HEIDI P |
|---|---:|---:|---:|
| ERAP2 | 0.0254694 | 0.06074163 | 0.2480520 |
| LIPA | 0.0106345 | 0.7502885 | 0.2011178 |

Neither gene meets the combined rule: neither SMR test reaches the association threshold, and neither HEIDI test detects heterogeneity. In the [two-dimensional plot](../results/05_smr/SMR_HEIDI.png), a larger x value indicates a stronger SMR association, and a larger y value indicates stronger evidence of HEIDI heterogeneity; “larger is better” does not apply to both axes. Missing HEIDI values remain NA.

This lesson uses the bundled SMR 1.4.3. The program, example inputs, and reference results are included in the course resources, and the run entry point selects the version matching the current system.

Exercises and answers:

1. How should you interpret significant SMR P and very small HEIDI P? **Answer:** An expression–disease association is detected, but the regional effects are heterogeneous; under the rule used in this example, this does not count as support for a shared signal.
2. What might change if `heidi_ld_upper` is changed from 0.9 to 0.8? **Answer:** The set of variants included in HEIDI may change. Compare the variant count and P values; do not choose parameters based on significance.
3. Can `p_HEIDI=NA` be counted as P>0.01? **Answer:** No. First check the number of instruments, regional coverage, alleles, and program output.
