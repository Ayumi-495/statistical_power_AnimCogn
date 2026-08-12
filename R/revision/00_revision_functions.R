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

# clubSandwich is a HARD requirement, not an optional enhancement. The Yang et al.
# (2024) specification is CR2 with Satterthwaite degrees of freedom (verified from
# the primary sources - see the block above `fit_fe_vcv_cr2()`), and CR2 is not
# implemented anywhere else in this repository. Falling back to the CR1 sandwich
# would silently change the reported specification, so this stops instead.
if (!requireNamespace("clubSandwich", quietly = TRUE)) {
  stop("clubSandwich is required for the Yang-2024 CR2 specification.\n",
       "  install.packages(\"clubSandwich\")\n",
       "  Do NOT substitute the CR1 sandwich: it is a different specification.")
}
PKG_VERSIONS <- c(
  R            = paste(R.version$major, R.version$minor, sep = "."),
  metafor      = as.character(utils::packageVersion("metafor")),
  clubSandwich = as.character(utils::packageVersion("clubSandwich")),
  lme4         = as.character(utils::packageVersion("lme4"))
)

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
#
# THE CRITICAL VALUE IS AN EXPLICIT ARGUMENT, and every caller records which one it
# used in a `crit_value_method` column. The default `CRIT` = qnorm(0.975) = 1.96 is
# the design-analysis convention of Gelman & Carlin (2014) and is what the submitted
# analysis uses. Adopting the Yang-2024 estimator introduces a second candidate -
# the Satterthwaite t critical value that clubSandwich reports for the CR2 test -
# and the two are not interchangeable: for a meta-analysis with ~10 study clusters
# qt(0.975, df_Satt) can exceed 2.4, which lowers power and raises Type M.
#
# CANONICAL CHOICE: z = 1.96 everywhere, so that a sensitivity analysis varies the
# assumed effect and its standard error and NOT the test used to define the metric.
# The Satterthwaite alternative is reported as an explicitly labelled diagnostic at
# the meta-analysis level, never mixed into the canonical rows.

power_two_tailed_cf <- function(mu, se, crit = CRIT) {
  t <- abs(mu) / se
  2 - stats::pnorm(crit - t) - stats::pnorm(crit + t)
}

type_S_cf <- function(mu, se, crit = CRIT) {
  t  <- abs(mu) / se
  pu <- 1 - stats::pnorm(crit - t)
  pl <- stats::pnorm(-crit - t)
  pl / (pu + pl)
}

# E[|X| | |X| > CRIT*se] / |mu| for X ~ N(mu, se). As t -> 0 this tends to
# dnorm(CRIT)/pnorm(-CRIT) / |t| = 2.3378 / |t|, i.e. it diverges. That divergence
# is why a corrected mean close to zero produces an arbitrarily large Type M, and
# it holds regardless of how the corrected mean was obtained.
type_M_cf <- function(mu, se, crit = CRIT) {
  t   <- mu / se
  num <- t * stats::pnorm(t - crit) + stats::dnorm(crit - t) -
         t * stats::pnorm(-crit - t) + stats::dnorm(crit + t)
  den <- abs(t) * (stats::pnorm(t - crit) + stats::pnorm(-crit - t))
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
  # THIS IS A DIAGNOSTIC, NOT THE REPORTED SPECIFICATION. Yang et al. (2024) use CR2
  # with Satterthwaite degrees of freedom; that is implemented in `fit_fe_vcv_cr2()`
  # below and is what the canonical tables report. CR0 / CR1 are retained only so the
  # effect of the small-sample correction is visible rather than assumed, and every
  # row derived from them carries verification_status = "single_derivation" and a
  # `role` beginning "diagnostic_". The unweighted median CR2/CR1 standard error ratio
  # is 1.008 across the 48 models, but it reaches 1.344, so the two are not
  # interchangeable in an aggregate.
  meat <- sum(tapply(a * e, cluster, sum)^2)
  v0   <- meat / (sum(a))^2
  v1   <- v0 * J / (J - 1)

  list(beta = beta, se_cr0 = sqrt(v0), se_cr1 = sqrt(v1),
       n_cluster = J, df = J - 1L,
       n_negative_weight = sum(a < 0), prop_negative_weight = mean(a < 0),
       y_min = min(y), y_max = max(y))
}

# --- Yang et al. (2024) as the source actually implements it ------------------
# VERIFIED FROM PRIMARY SOURCES on 2026-08-12. Do not re-derive from notes.
#
# Provenance of the sources. The paper's data-availability statement names two
# locations: the analysis repository `github.com/Yefeng0920/WLS_RVE` and the tutorial
# `yefeng0920.github.io/BiasRobustMA_tutorial/`. Both were fetched fresh from GitHub:
#   BiasRobustMA_tutorial/R/hands_on_R.Rmd  md5 24f70634e4157e63321faac86e38f3e7
#   WLS_RVE/R/final_GLS_V5.Rmd             md5 ddfdf14a1840070e0956ae7c794c283e
#
# The tutorial's step one and step two, verbatim:
#   VCV <- vcalc(vi = var.eff.size, cluster = study, rho = 0.5, obs = obs, data = dat)
#   mod_MLFE     <- rma.mv(yi = eff.size, V = VCV, method = "REML", test = "t",
#                          dfs = "contain", data = dat)
#   mod_MLFE_RVE <- robust(mod_MLFE, cluster = study, adjust = TRUE, clubSandwich = TRUE)
# The paper's own re-analysis of the 448 meta-analyses uses the same two calls
# (final_GLS_V5.Rmd:271 and :281; the CRVE call there omits `adjust`, which
# `robust()` ignores on the clubSandwich path).
#
# WHAT `clubSandwich = TRUE` MEANS. Read directly out of metafor 5.0.1's
# `robust.rma.mv` source rather than inferred from documentation: the defaults are
#   vcov      = "CR2"
#   coef_test = "Satterthwaite"      (and conf_test inherits coef_test)
# so the specification is CR2 with Satterthwaite degrees of freedom. The tutorial
# states the reason in a technical note: "CR2 correction performs better than CR1.
# However, CR2 is not applicable to models with non-nested random effects ... In our
# case, the model in the first step does not include random effects."
#
# EXTERNAL ANCHOR. `06_validate_yang2024_reference.R` runs this same code path on the
# tutorial's own example data and reproduces every digit of the tutorial's rendered
# output (step one beta 0.074, SE 0.018, CI [0.039, 0.108]; step two SE 0.053,
# p 0.168, CI [-0.032, 0.18], df 1632 -> 52). The prose of the paper quotes slightly
# different figures for the same example (beta 0.075, SE 0.054, t 1.375, p 0.175,
# CI [-0.034, 0.184]); those come from the 448-model corpus pipeline, which applies
# its own exclusions, not from the tutorial's data preparation. Both are CR2 with
# Satterthwaite df and they agree to two decimals; the tutorial is the reproducible
# one and is therefore the anchor used here.
#
# Returns NULL on failure so the caller can record a per-model status rather than
# silently dropping a meta-analysis from the 48.
fit_fe_vcv_cr2 <- function(y, vi, cluster, rho = RHO_DEFAULT) {
  d <- data.frame(es = y, var = vi, study_ID = as.character(cluster),
                  stringsAsFactors = FALSE)
  d$obs <- seq_len(nrow(d))

  V <- try(metafor::vcalc(vi = var, cluster = study_ID, obs = obs, rho = rho,
                          data = d), silent = TRUE)
  if (inherits(V, "try-error")) return(NULL)

  # Third check on our own VCV builder: metafor's vcalc and build_vcv() must agree.
  vcv_diff <- max(abs(as.matrix(V) - build_vcv(d$var, d$study_ID, rho)))

  fit <- try(metafor::rma.mv(yi = es, V = V, method = "REML", test = "t",
                             dfs = "contain", data = d, sparse = TRUE),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  # Canonical path: metafor's wrapper, exactly as the tutorial calls it.
  rve <- try(metafor::robust(fit, cluster = d$study_ID, adjust = TRUE,
                             clubSandwich = TRUE), silent = TRUE)
  if (inherits(rve, "try-error")) return(NULL)

  # Cross-check: clubSandwich called directly, bypassing metafor's wrapper. Same
  # package, different call path - this catches argument-passing errors in the
  # wrapper, which is the realistic failure mode. It is NOT a second derivation of
  # the CR2 algebra itself; that is what the external anchor above is for.
  ct <- try(clubSandwich::coef_test(fit, vcov = "CR2", cluster = d$study_ID,
                                    test = "Satterthwaite"), silent = TRUE)
  ci <- try(clubSandwich::conf_int(fit, vcov = "CR2", cluster = d$study_ID,
                                   test = "Satterthwaite"), silent = TRUE)
  if (inherits(ct, "try-error") || inherits(ci, "try-error")) return(NULL)

  list(
    beta = as.numeric(fit$b[1]), se_cr2 = as.numeric(rve$se[1]),
    df_satt = as.numeric(rve$ddf[1]), pval = as.numeric(rve$pval[1]),
    ci_lb = as.numeric(rve$ci.lb[1]), ci_ub = as.numeric(rve$ci.ub[1]),
    se_working = as.numeric(fit$se[1]),
    vcv_max_abs_diff = vcv_diff,
    wrapper_max_abs_diff = max(abs(c(rve$se[1] - ct$SE, rve$ddf[1] - ct$df_Satt,
                                     rve$ci.lb[1] - ci$CI_L, rve$ci.ub[1] - ci$CI_U)))
  )
}

# UWLS, matching the source's OWN UWLS specification, which differs from its
# FE + VCV specification and must not be copied across from it. In
# final_GLS_V5.Rmd:402 and :451 the carried-forward objects are
#   lm(eff.size ~ 1, weights = 1/var.eff.size)
#   coef_test(mod, vcov = "CR2", cluster = study, test = "naive-t",
#             target = var.eff.size)
# i.e. CR2 with naive-t degrees of freedom and an explicit `target` working variance,
# not Satterthwaite. (Satterthwaite variants appear at :491 and :497 as exploratory
# code, one annotated "too conservative", and are not the objects used downstream.)
# UWLS is supplementary here, so it follows its source specification rather than
# being forced into the FE + VCV one.
fit_uwls_cr2 <- function(y, vi, cluster) {
  d <- data.frame(es = y, var = vi, study_ID = as.character(cluster),
                  stringsAsFactors = FALSE)
  m <- try(stats::lm(es ~ 1, weights = 1 / var, data = d), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  ct <- try(clubSandwich::coef_test(m, vcov = "CR2", cluster = d$study_ID,
                                    test = "naive-t", target = d$var), silent = TRUE)
  ci <- try(clubSandwich::conf_int(m, vcov = "CR2", cluster = d$study_ID,
                                   test = "naive-t", target = d$var), silent = TRUE)
  if (inherits(ct, "try-error") || inherits(ci, "try-error")) return(NULL)
  list(beta = as.numeric(stats::coef(m)[[1]]), se_cr2 = as.numeric(ct$SE),
       df = as.numeric(ct$df_t), pval = as.numeric(ct$p_t),
       ci_lb = as.numeric(ci$CI_L), ci_ub = as.numeric(ci$CI_U))
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
