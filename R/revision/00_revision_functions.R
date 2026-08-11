# R/revision/00_revision_functions.R --------------------------------------------
# Shared functions for the revision sensitivity analyses.
#
# Sourced by 01-05. Defines nothing that touches the submitted analysis:
# `S2_v2.R` is never sourced, edited or run by anything in R/revision/.
#
# Two deliberate differences from R/00_setup.R, both documented here because they
# change reported numbers:
#
#   1. Type M is computed in CLOSED FORM, not by Monte Carlo. `S2_v2.R` and
#      R/00_setup.R draw 10,000 normal deviates per unit, which makes results
#      depend on the RNG state and quantises Type M to 3 decimals. The closed form
#      is exact, deterministic, and agrees with the simulation to 0.34% at
#      N = 200,000. Every Type M in results/revision/ is closed-form.
#   2. Data loading and the optimizer table are REUSED from R/00_setup.R rather
#      than duplicated, so the two cannot drift apart.

suppressMessages({
  library(here); library(dplyr); library(tibble); library(readr)
  library(metafor); library(lme4)
})

source(here::here("R", "00_setup.R"))   # load_datasets(), dataset_index(), fit_rma(), optimizer_for()

REV_OUT <- here::here("results", "revision")
# Chaining artifacts (fitted objects passed from one script to the next) live in a
# separate, gitignored subdirectory so that results/revision/ contains only the
# canonical result tables and the README.
REV_TMP <- file.path(REV_OUT, "intermediate")
dir.create(REV_TMP, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
CRIT  <- stats::qnorm(1 - ALPHA / 2)     # 1.959964

# --- design analysis metrics, closed form ------------------------------------
# All three are functions of t = mu/se alone. Verified: holding t fixed and
# varying se over {0.01, 0.1, 1, 10} leaves all three unchanged to 10 decimals.

power_two_tailed_cf <- function(mu, se) {
  t <- abs(mu) / se
  2 - stats::pnorm(CRIT - t) - stats::pnorm(CRIT + t)
}

type_S_cf <- function(mu, se) {
  t  <- abs(mu) / se
  pu <- 1 - stats::pnorm(CRIT - t)
  pl <- stats::pnorm(-CRIT - t)
  pl / (pu + pl)
}

# E[|X| | |X| > CRIT*se] / |mu| for X ~ N(mu, se). As t -> 0 this tends to
# dnorm(CRIT)/pnorm(-CRIT) / |t| = 2.3378 / |t|, i.e. it diverges. That divergence
# is why a corrected mean close to zero produces an arbitrarily large Type M, and
# it holds regardless of how the corrected mean was obtained.
type_M_cf <- function(mu, se) {
  t   <- mu / se
  num <- t * stats::pnorm(t - CRIT) + stats::dnorm(CRIT - t) -
         t * stats::pnorm(-CRIT - t) + stats::dnorm(CRIT + t)
  den <- abs(t) * (stats::pnorm(t - CRIT) + stats::pnorm(-CRIT - t))
  num / den
}

# --- Yang et al. (2024) step one: bias-robust point estimation ---------------
# Yang Y, Lagisz M, Williams C, Noble DWA, Pan J, Nakagawa S (2024). Robust point
# and variance estimation for meta-analyses with selective reporting and dependent
# effect sizes. Methods in Ecology and Evolution 15(9). doi:10.1111/2041-210X.14377
#
# The paper's supported claim is that the two-step approach "does not rely on
# extrapolation" (Section 5.1, p.1605), in contrast to PET-PEESE, whose intercept
# is a marginalised mean "assuming an infinitely large sample size". It is NOT
# claimed anywhere in the paper that the estimator cannot cross zero, and
# empirically it does: see results/revision/reversal_counts.csv.
#
# VCV construction (paper Equation 4). Block-diagonal by CLUSTERING UNIT, which
# here is `study_ID` (the same grouping the submitted analysis uses as the outer
# random effect in ~1|study_ID/obs_ID):
#
#   diagonal    V[i,i] = v_i                       (the effect size's sampling variance)
#   off-diagonal V[i,j] = rho * sqrt(v_i * v_j)    for i != j in the SAME study
#   between studies      0
#
# rho is the assumed constant within-study sampling correlation. The tutorial uses
# rho = 0.5; RHO_DEFAULT below matches it, and 03_yang2024_bias_robust.R reports a
# sensitivity analysis over rho in {0, 0.25, 0.5, 0.75}.
RHO_DEFAULT <- 0.5

build_vcv <- function(vi, cluster, rho = RHO_DEFAULT) {
  n <- length(vi)
  V <- diag(vi, nrow = n)
  s <- sqrt(vi)
  for (lev in unique(cluster)) {
    ix <- which(cluster == lev)
    if (length(ix) > 1L) {
      for (a in ix) for (b in ix) if (a != b) V[a, b] <- rho * s[a] * s[b]
    }
  }
  V
}

# FE + VCV (paper Equation 3): a fixed-effect (no random effects) GLS intercept.
#   beta = (1' V^-1 1)^-1 1' V^-1 y, i.e. a weighted mean with weights a = V^-1 1.
# Note the weights are NOT constrained to be positive. The paper states this
# directly (Section 2.2.1, p.1596): "when turning the VCV matrix into a weighting
# scheme, the off-diagonal elements become negative values, which can de-emphasize
# the studies reporting more effect size estimates". So FE + VCV carries no
# convex-combination guarantee; `n_negative_weight` is recorded per meta-analysis.
#
# UWLS (Stanley & Doucouliagos): the same estimator with a DIAGONAL VCV, i.e.
# rho = 0, so a = 1/v_i and the estimate is strictly a convex combination of the
# observed effect sizes. The paper notes (Section 2.2.1) that UWLS, Henmi-Copas and
# IVhet are special cases of its framework with zero sampling correlation.
#
# Implemented in closed form rather than via rma.mv so the weights are inspectable
# and the result cannot depend on an optimizer. Agreement with a metafor fit was
# checked at 1.1e-04, the residual being a Hedges' g reimplementation difference in
# the five SMD/des_stat datasets.
bias_robust_fit <- function(y, vi, cluster, rho = RHO_DEFAULT, diagonal = FALSE) {
  n <- length(y)
  if (diagonal) {
    a <- 1 / vi                                    # UWLS
  } else {
    V  <- build_vcv(vi, cluster, rho)
    Vi <- tryCatch(solve(V), error = function(e) NULL)
    if (is.null(Vi)) return(NULL)                  # non-invertible VCV
    a <- as.vector(Vi %*% rep(1, n))               # FE + VCV
  }
  beta <- sum(a * y) / sum(a)
  e    <- y - beta
  J    <- dplyr::n_distinct(cluster)

  # Step two: cluster-robust variance estimation (paper Equation 5), clustered by
  # `study_ID`. Sandwich estimator:
  #   Var(beta) = (sum a)^-2 * sum_over_clusters ( sum_{i in cluster} a_i e_i )^2
  # CR0 is the unadjusted sandwich; CR1 applies the small-sample factor J/(J-1).
  # Intervals below use t on J-1 degrees of freedom.
  #
  # NOT IMPLEMENTED HERE: the CR2 / Satterthwaite correction (clubSandwich, and
  # metafor::robust(..., clubSandwich = TRUE)). An earlier pass used CR2 and its
  # standard errors differ from CR1 by up to 0.036 on standard errors of order
  # 0.1-0.2. That variant has not been independently replicated, so it is excluded
  # from the canonical tables and the CRVE columns are marked
  # verification_status = "single_derivation".
  meat <- sum(tapply(a * e, cluster, sum)^2)
  v0   <- meat / (sum(a))^2
  v1   <- v0 * J / (J - 1)

  list(beta = beta, se_cr0 = sqrt(v0), se_cr1 = sqrt(v1),
       n_cluster = J, df = J - 1L,
       n_negative_weight = sum(a < 0), prop_negative_weight = mean(a < 0),
       y_min = min(y), y_max = max(y))
}

# --- aggregation ------------------------------------------------------------
# Two aggregation schemes, matching the submitted analysis so the revision rows
# are comparable with the published ones. Both are reported with explicit
# provenance in the canonical tables (`aggregation`, `weighting` columns).
#
#   meta-analysis level : lm(log(metric) ~ 1, weights = k), k = number of effect
#                         sizes in the meta-analysis. 48 units.
#   primary-study level : lmer(log(metric) ~ 1 + (1 | study_ID)), unweighted.
#                         5,740 units.
#
# The reported summary is the back-transformed intercept. Note this is a weighted
# GEOMETRIC mean, which equals the median only under log-symmetry; the submitted
# manuscript calls it a median. Column name below says geometric_mean to avoid
# repeating that. The offset for Type S (which can be exactly zero) is 0.025,
# subtracted after back-transformation, as in the submitted analysis; the resulting
# lower confidence limit can be negative and is reported unfloored.
TYPE_S_OFFSET <- 0.025
offset_for    <- function(metric) if (metric == "type_S") TYPE_S_OFFSET else 0

aggregate_ma <- function(values, k, metric) {
  off <- offset_for(metric)
  fit <- stats::lm(log(values + off) ~ 1, weights = k)
  ci  <- stats::confint(fit)
  tibble::tibble(
    geometric_mean = exp(stats::coef(fit)[[1]]) - off,
    ci_lower = exp(ci[1]) - off, ci_upper = exp(ci[2]) - off,
    raw_median = stats::median(values), raw_min = min(values), raw_max = max(values),
    arithmetic_mean_unweighted = mean(values),
    arithmetic_mean_kweighted  = sum(values * k) / sum(k),
    n_unit = length(values)
  )
}

aggregate_primary <- function(values, study_id, metric) {
  off <- offset_for(metric)
  d   <- data.frame(y = log(values + off), study_ID = study_id)
  fit <- lme4::lmer(y ~ 1 + (1 | study_ID), data = d)
  ci  <- stats::confint(fit, method = "Wald")
  ci  <- ci[rownames(ci) == "(Intercept)", ]
  tibble::tibble(
    geometric_mean = exp(lme4::fixef(fit)[[1]]) - off,
    ci_lower = exp(ci[[1]]) - off, ci_upper = exp(ci[[2]]) - off,
    raw_median = stats::median(values), raw_min = min(values), raw_max = max(values),
    arithmetic_mean_unweighted = mean(values),
    arithmetic_mean_kweighted  = NA_real_,   # the primary-level model is unweighted
    n_unit = length(values)
  )
}

write_revision <- function(x, filename) {
  p <- file.path(REV_OUT, filename)
  readr::write_csv(x, p)
  message("wrote ", filename, "  (", nrow(x), " rows)")
  invisible(p)
}
