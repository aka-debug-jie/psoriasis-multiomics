# Small statistical primitives; no data selection or discovery policy is hidden here.
wald_ratio <- function(beta_exposure, se_exposure, beta_outcome, se_outcome) {
  stopifnot(all(is.finite(beta_exposure)), all(beta_exposure != 0),
            all(is.finite(beta_outcome)), all(se_exposure > 0), all(se_outcome > 0))
  beta <- beta_outcome / beta_exposure
  se <- se_outcome / abs(beta_exposure)
  z <- beta / se
  log_p <- log(2) + pnorm(-abs(z), log.p = TRUE)
  data.frame(beta = beta, se = se, OR = exp(beta),
             CI_low = exp(beta - qnorm(.975) * se),
             CI_high = exp(beta + qnorm(.975) * se),
             p = exp(log_p), neg_log10_p = -log_p / log(10),
             F_statistic = (beta_exposure / se_exposure)^2)
}

fixed_effect <- function(beta, se) {
  stopifnot(length(beta) == length(se), all(is.finite(beta)), all(se > 0))
  weights <- 1 / se^2
  estimate <- sum(weights * beta) / sum(weights)
  standard_error <- sqrt(1 / sum(weights))
  log_p <- log(2) + pnorm(-abs(estimate / standard_error), log.p = TRUE)
  data.frame(beta = estimate, se = standard_error, OR = exp(estimate),
             CI_low = exp(estimate - qnorm(.975) * standard_error),
             CI_high = exp(estimate + qnorm(.975) * standard_error),
             p = exp(log_p), neg_log10_p = -log_p / log(10),
             n_cohorts = length(beta), Q = sum(weights * (beta - estimate)^2))
}

# Orientation only. The caller must first establish a unique, build-matched variant.
# Palindromic alleles require additional evidence and are deliberately left unresolved.
orient_alleles <- function(exposure_ea, exposure_oa, outcome_ea, outcome_oa, beta, eaf = NA_real_) {
  a <- toupper(c(exposure_ea, exposure_oa, outcome_ea, outcome_oa))
  if (anyNA(a) || any(!a %in% c("A", "C", "G", "T")))
    return(list(status = "requires_normalized_variant", beta = NA_real_, eaf = NA_real_))
  if (a[1] == a[2] || a[3] == a[4])
    return(list(status = "invalid_alleles", beta = NA_real_, eaf = NA_real_))
  if (paste(sort(a[1:2]), collapse = "") %in% c("AT", "CG"))
    return(list(status = "palindromic_unresolved", beta = NA_real_, eaf = NA_real_))
  complement <- c(A = "T", T = "A", C = "G", G = "C")
  sign <- NA_real_
  status <- "allele_mismatch"
  if (all(a[1:2] == a[3:4])) { sign <- 1; status <- "aligned" }
  else if (all(a[1:2] == a[4:3])) { sign <- -1; status <- "swapped" }
  else if (all(a[1:2] == unname(complement[a[3:4]]))) { sign <- 1; status <- "strand_aligned" }
  else if (all(a[1:2] == unname(complement[a[4:3]]))) { sign <- -1; status <- "strand_swapped" }
  list(status = status, beta = beta * sign,
       eaf = if (is.na(sign)) NA_real_ else if (sign == -1) 1 - eaf else eaf)
}
