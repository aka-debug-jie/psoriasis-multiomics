# TwoSampleMR method sources

`mr.R`, `mr_mode.R`, and `add_rsq.R` retain the fixed [MRCIEU/TwoSampleMR source](https://github.com/MRCIEU/TwoSampleMR/tree/d8b07182da23fb984157d9335fefcb83e894a944/R) at commit `d8b07182da23fb984157d9335fefcb83e894a944`. Chapter 04 calls IVW, weighted median, weighted mode, MR-Egger, and binary-trait correlation conversion functions. Resampling uses the recorded seed and 1,000 bootstrap replicates.

Authors: Gibran Hemani and contributors. The package declares MIT + LICENSE. The original copyright notice is retained in `TwoSampleMR_LICENSE`, with the complete MIT text in `MIT.txt`. Method reference: Hemani et al. *The MR-Base platform supports systematic causal inference across the human phenome*. eLife 2018;7:e34408. https://doi.org/10.7554/eLife.34408.
