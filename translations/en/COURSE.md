# Statistical functions and references

<a id="appendix-genetics"></a>

## Appendix A | Complete implementation of shared genetic-statistics functions

Chapters 2, 3, 4, 6, 8, and 11 load the functions below with `source("R/genetics.R")`. The [source code](R/genetics.R) is shown here so you can trace a function call through to its calculation. Read in this order: first see how alleles change the sign, then how duplicate records for the same variant are handled, and finally the matrix formula for generalized IVW. The main lesson code already loads these functions; you do not need to copy this appendix to run it.

| Function | Inputs and purpose |
|---|---|
| `allele_sign()` | Returns 1, −1, or NA for a pair of allele pairs; under the current rule, palindromic sites remain unresolved when the direction cannot be uniquely harmonized. |
| `load_ld()` / `ld_for()` | Reads within-protein correlations and converts their signs according to the exposure effect alleles. |
| `freeze_instruments()` | Deduplicates first by exposure F and a fixed order, then selects instruments by LD; the result does not depend on disease P values. |
| `givw()` | Constructs the disease-effect covariance matrix S and calculates the common effect, SE, Q, and dispersion sensitivity. |
| `analyse_outcome()` | Aligns one outcome, records available instruments, then calls the estimator for each protein. |

In `givw()`, `solve(S, x)` directly solves a linear system. The denominator `sum(x * solve(S, x))` represents the information contributed by the current instruments under the covariance structure. With one instrument, generalized IVW reduces to the Wald ratio. `se_mre` inflates the standard error by `max(1, Q/df)`, so under this rule heterogeneity cannot make the interval narrower; it is reported separately from the primary fixed-effect SE. This first-order SE does not include uncertainty in the exposure estimate; the delta formula in Chapter 2 shows another approximation for a single instrument.

```r
allele_sign <- function(ea, oa, target_ea, target_oa) {
    a <- toupper(c(ea, oa, target_ea, target_oa))
    if (anyNA(a) || any(!a %in% c("A", "C", "G", "T")) || a[1] == a[2] || a[3] == a[4]) 
        return(NA_real_)
    if (paste(sort(a[1:2]), collapse = "") %in% c("AT", "CG")) 
        return(NA_real_)
    if (identical(a[1:2], a[3:4])) 
        return(1)
    if (identical(a[1:2], rev(a[3:4]))) 
        return(-1)
    complement <- c(A = "T", T = "A", C = "G", G = "C")[a[3:4]]
    if (identical(a[1:2], unname(complement))) 
        return(1)
    if (identical(a[1:2], unname(rev(complement)))) 
        return(-1)
    NA_real_
}
load_ld <- function() {
    p <- read_tsv("data/genetics/within_protein_ld.tsv")
    setNames(p$r, paste(p$i, p$j, sep = ":"))
}
ld_for <- function(d, ld) {
    if (nrow(d) == 1) 
        return(matrix(1, 1, 1))
    ix <- outer(d$LD_index, d$LD_index, paste, sep = ":")
    matrix(unname(ld[ix]), nrow(d)) * outer(d$LD_sign, d$LD_sign)
}
freeze_instruments <- function(x, ld, max_abs_r = 0.95) {
    x <- x[order(-x$F_statistic, x$cohort, x$assay, x$rsid), ]
    x <- x[!duplicated(x[c("uniprot", "variant_key")]), ]
    out <- lapply(split(x, x$uniprot), function(d) {
        d <- d[order(-d$F_statistic, d$cohort, d$assay, d$rsid), ]
        keep <- 1L
        reason <- "single_unique_variant"
        if (nrow(d) > 1) {
            r <- ld_for(d, ld)
            if (any(d$LD_index < 0) | any(!is.finite(r))) 
                reason <- "single_IV_missing_or_unresolved_LD"
            else {
                keep <- integer()
                for (j in seq_len(nrow(d))) if (!length(keep) || all(abs(r[j, keep]) <= max_abs_r)) 
                  keep <- c(keep, j)
                z <- r[keep, keep, drop = FALSE]
                if (min(eigen(z, symmetric = TRUE, only.values = TRUE)$values) <= 0 || kappa(z, 
                  exact = TRUE) > 10000) {
                  keep <- 1L
                  reason <- "single_IV_unstable_LD"
                }
                else reason <- if (length(keep) < nrow(d)) 
                  "LD_pruned"
                else "all_unique_variants"
            }
        }
        d$selected <- seq_len(nrow(d)) %in% keep
        d$selection_reason <- reason
        d
    })
    do.call(rbind, out)
}
givw <- function(x, y, sy, r) {
    S <- outer(sy, sy) * r
    wx <- solve(S, x)
    denominator <- sum(x * wx)
    beta <- sum(wx * y)/denominator
    se <- sqrt(1/denominator)
    residual <- y - beta * x
    Q <- as.numeric(t(residual) %*% solve(S, residual))
    df <- length(x) - 1L
    se_mre <- se * sqrt(if (df > 0) max(1, Q/df) else 1)
    data.frame(beta = beta, se = se, p = 2 * pnorm(-abs(beta/se)), CI_log_low = beta - qnorm(0.975) * 
        se, CI_log_high = beta + qnorm(0.975) * se, Q = Q, Q_df = df, Q_p = if (df > 0) 
        pchisq(Q, df, lower.tail = FALSE)
    else NA_real_, se_mre = se_mre, p_mre = 2 * pnorm(-abs(beta/se_mre)), n_variants = length(x))
}
analyse_outcome <- function(plan, y, ld) {
    stopifnot(!anyDuplicated(y$rsid))
    ans <- list()
    audits <- list()
    for (uid in unique(plan$uniprot)) {
        d <- plan[plan$uniprot == uid & plan$selected, ]
        a <- y[match(d$rsid, y$rsid), ]
        sign <- mapply(allele_sign, d$EA, d$OA, a$EA, a$OA)
        d$y <- a$beta * sign
        d$sy <- a$se
        d$included <- is.finite(d$y) & is.finite(d$sy) & d$sy > 0 & !grepl("^ambiguous_aliases", 
            d$selection_status)
        d$orientation <- sign
        audits[[uid]] <- d
        v <- d[d$included, ]
        row <- data.frame(uniprot = uid, gene = d$gene[1], frozen_variants = nrow(d), observed_variants = nrow(v))
        if (nrow(v)) 
            row <- cbind(row, givw(v$beta, v$y, v$sy, ld_for(v, ld)))
        ans[[uid]] <- row
    }
    list(results = as.data.frame(data.table::rbindlist(ans, fill = TRUE)), harmonized = do.call(rbind, 
        audits))
}
```

<a id="appendix-common"></a>

## Appendix B | Reading tables, plotting, and external programs

The run entry point has already selected R and the course package library. Shared functions read tables, save results, and create plots; statistical modeling remains in each lesson's code. The implementation below matches the source files shipped with the course.

```r
source('R/load_environment.R')
# 在 teaching.Rproj 打开的项目根目录中运行各课脚本。
options(stringsAsFactors = FALSE, scipen = 6)
Sys.setenv(OMP_NUM_THREADS = "2", OPENBLAS_NUM_THREADS = "2", MKL_NUM_THREADS = "2")
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(2)
read_tsv <- function(path) read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path,
  check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
  quote = FALSE, na = "NA", fileEncoding = "UTF-8")
lesson_output <- function(id) {
  path <- file.path(Sys.getenv("TEACHING_OUTPUT", "results"), id)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
theme_course <- function() ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(colour = "#566573"),
    axis.line = ggplot2::element_line(linewidth = .35),
    legend.position = "bottom", plot.margin = ggplot2::margin(12, 16, 12, 12))
save_plot <- function(p, output, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(output, paste0(name, ".pdf")), p, width = width, height = height)
  ggplot2::ggsave(file.path(output, paste0(name, ".png")), p, width = width, height = height,
                 dpi = 320, bg = "white")
  if(isTRUE(getOption("teaching.inline",FALSE))) {
    options(repr.plot.width=width,repr.plot.height=height)
    print(p)
  }
}
source('R/tool_support.R')
run_tool <- function(program, args, log) {
  code <- system2(program, vapply(args, shQuote, ""), stdout = log, stderr = log)
  if (code != 0) stop("外部程序运行失败；请查看 ", log)
  invisible(code)
}

```

`tool_path()` selects the bundled PLINK or SMR program based on the operating system and processor; see [R/tool_support.R](R/tool_support.R) for program-location code. Results are written to `resources/work/results`, with logs and outputs in the folder for the corresponding lesson.

<a id="appendix-methods"></a>

## Appendix C | Third-party methods, data citations, and further reading

Chapter 4 calls fixed-version TwoSampleMR method source code stored under `R/vendor/`. The main lesson shows the function arguments in full; see [method sources](R/vendor/README.md) for the original implementations, copyright, and version notes.

| Call | Use in this course | What to look for |
|---|---|---|
| `mr_ivw()` | Multi-instrument IVW estimation | Weights, heterogeneity, and standard-error definitions |
| `mr_weighted_median()` | Weighted-median estimation | Proportion of valid instrument weight and bootstrap |
| `mr_weighted_mode()` | Weighted-mode estimation | Ratio distribution, bandwidth, and largest consistent group |
| `mr_egger_regression()` | MR regression with an intercept | Intercept, instrument direction, and additional identification assumptions |
| `get_r_from_lor()` | Converts a binary-outcome effect to a correlation | Population prevalence, case/control counts, and allele frequency |
| `coloc::coloc.abf()` | Posterior comparison across five regional hypotheses | Regional completeness, effect variance, and priors |
| `survival::coxph()` / `cox.zph()` | Cox estimation and PH test | Time scale, event coding, covariates, and risk sets |

Look up full package-function arguments in R help pages, for example `?limma::voom`, `?coloc::coloc.abf`, and `?survival::coxph`. Method references include Hemani et al.'s [MR-Base paper](https://doi.org/10.7554/eLife.34408), the [official SMR documentation](https://yanglab.westlake.edu.cn/software/smr/), and the [PLINK documentation](https://www.cog-genomics.org/plink/1.9/). Specific versions and access points for all data are recorded in the [source manifest](data/sources.tsv); use remains subject to each provider's terms.

Choose a chapter from the [course start page](notebooks/00_start.nb.html), run the code in the order shown for each lesson, and complete the exercises.
