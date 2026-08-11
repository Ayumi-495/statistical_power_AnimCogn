# 01_estimates.R -------------------------------------------------------------
# Per-meta-analysis effect estimates: the uncorrected mean (beta0), the
# publication-bias detection model, the four-scenario bias correction, and the
# selected corrected mean (beta0_c3).
#
# This re-implements S2_v2.R:284-1111 as one loop over the 48 datasets instead
# of twelve hand-written per-metric/per-scenario blocks. The specifications are
# unchanged; only the control flow is. Three things are worth stating because
# they were not documented in the manuscript:
#
#   Scenario assignment uses the SIGN of beta0*beta1 and beta0*beta2 only --
#   NOT significance. The p < 0.05 tests in S2_v2.R:429-443 are used to report
#   which meta-analyses "have" a small-study or decline effect; they do not
#   affect which correction model is fitted.
#
#   beta0_c comes from the sampling-ERROR moderator model, beta0_c2 from the
#   sampling-VARIANCE moderator model, and the reported corrected mean beta0_c3
#   is beta0_c2 only when |beta0| > |beta0_c2|, otherwise it reverts to beta0.
#   That gate is one-directional: correction can shrink a magnitude but never
#   grow one, and "no bias detected" and "correction moved the estimate away
#   from zero" both collapse to beta0.
#
#   S2_v2.R:544 and :554 write `mod = sei + year_pub.l` and
#   `mod = var + year_pub.l`, omitting the `~`. metafor partially matches `mod`
#   to `mods` but the value is then an EXPRESSION, not a formula, so a single
#   composite moderator equal to the arithmetic sum of the two variables is
#   fitted instead of two separate moderators. The typo is confined to the lnRR
#   scenario-1 block, whose only member is MA09.csv -- the largest lnRR dataset
#   (1,297 of the 5,740 effect-size rows, 126 studies). Its corrected mean is
#   0.1060 under the composite moderator and 0.0681 under the intended
#   two-moderator model.
#
#   Both are produced here. `legacy = TRUE` reproduces the submitted
#   manuscript; `legacy = FALSE` (the default, used for all revision analyses)
#   fits the intended specification. S2_v2.R itself is never modified.

source(here::here("R", "00_setup.R"))

err_mod <- function(uses_ess) if (uses_ess) "ess.sei" else "sei"
var_mod <- function(uses_ess) if (uses_ess) "ess.var" else "var"

# Scenario-specific moderator structure.
#   1: error/variance moderator + year        2: year only
#   3: error/variance moderator only          4: intercept only
mods_for <- function(scenario, moderator) {
  switch(as.character(scenario),
    "1" = stats::as.formula(paste("~", moderator, "+ year_pub.l")),
    "2" = ~ year_pub.l,
    "3" = stats::as.formula(paste("~", moderator)),
    "4" = NULL)
}

# The baseline's composite-moderator variant, reproduced by evaluating
# `moderator + year_pub.l` as a numeric vector rather than a formula.
fit_composite <- function(data, moderator) {
  m <- data[[moderator]] + data$year_pub.l
  d <- data
  d$.composite <- m
  fit_rma(d, mods = ~ .composite, optimizer = "optim")
}

build_estimates <- function(dat, idx, legacy = FALSE) {
  L <- all_datasets(dat)
  stopifnot(identical(names(L), idx$case))

  # 1. uncorrected mean ------------------------------------------------------
  message("  intercept-only models (48)")
  f_null <- lapply(seq_along(L), function(i)
    fit_rma(L[[i]], optimizer = optimizer_for(idx$es_type[i], "null")))
  stopifnot("intercept-only fit failed" = !any(vapply(f_null, is.null, logical(1))))

  est <- idx |>
    dplyr::mutate(
      beta0      = vapply(f_null, function(m) as.numeric(m$beta), numeric(1)),
      se_beta0   = vapply(f_null, function(m) m$se, numeric(1)),
      pval_beta0 = vapply(f_null, function(m) m$pval, numeric(1)),
      # t-based CI read from the fitted object. NOT beta0 +/- 1.96*se: the
      # models use test = "t", ddf ranges 3-1296, and a normal-quantile
      # reconstruction differs by up to 0.825 on the effect-size scale.
      ci_lb = vapply(f_null, function(m) m$ci.lb, numeric(1)),
      ci_ub = vapply(f_null, function(m) m$ci.ub, numeric(1)),
      ddf   = vapply(f_null, function(m) as.numeric(m$ddf)[1], numeric(1)),
      k = vapply(L, nrow, integer(1)),
      n_study_id = vapply(L, function(d) dplyr::n_distinct(d$study_ID), integer(1))
    )

  # 2. bias detection -------------------------------------------------------
  # For the five SMD datasets built from raw descriptive statistics the
  # moderator is the effective-sample-size analogue ess.sei, which breaks the
  # artefactual correlation between an SMD estimate and its own sampling error.
  # lnRR has the same artefactual correlation but its datasets carry no group
  # sample sizes, so ess cannot be computed -- a data-limited deviation from
  # Yang et al. that applies to 5 of the 48 models.
  message("  bias-detection models (48)")
  f_det <- lapply(seq_along(L), function(i) {
    f <- stats::as.formula(paste("~", err_mod(idx$uses_ess[i]), "+ year_pub.l"))
    fit_rma(L[[i]], mods = f,
            optimizer = optimizer_for(idx$es_type[i], "detect", idx$uses_ess[i]))
  })
  stopifnot("detection fit failed" = !any(vapply(f_det, is.null, logical(1))))

  est <- est |>
    dplyr::mutate(
      beta1 = vapply(f_det, function(m) as.numeric(m$beta[2]), numeric(1)),
      pval_beta1 = vapply(f_det, function(m) m$pval[2], numeric(1)),
      beta2 = vapply(f_det, function(m) as.numeric(m$beta[3]), numeric(1)),
      pval_beta2 = vapply(f_det, function(m) m$pval[3], numeric(1)),
      beta0_x_beta1 = beta0 * beta1,
      beta0_x_beta2 = beta0 * beta2,
      has_small_study_effect = pval_beta1 < 0.05 & beta0_x_beta1 > 0,
      has_decline_effect     = pval_beta2 < 0.05 & beta0_x_beta2 < 0,
      scenario = dplyr::case_when(
        beta0_x_beta1 > 0 & beta0_x_beta2 < 0 ~ 1L,
        beta0_x_beta1 < 0 & beta0_x_beta2 < 0 ~ 2L,
        beta0_x_beta1 > 0 & beta0_x_beta2 > 0 ~ 3L,
        TRUE                                  ~ 4L
      )
    )

  # 3. scenario-specific correction models ----------------------------------
  message("  correction models (48 x 2)")
  use_composite <- legacy & est$es_type == "lnRR" & est$scenario == 1L

  corr <- lapply(seq_along(L), function(i) {
    sc <- est$scenario[i]
    if (use_composite[i]) {
      list(err = fit_composite(L[[i]], err_mod(est$uses_ess[i])),
           var = fit_composite(L[[i]], var_mod(est$uses_ess[i])))
    } else {
      m_err <- fit_rma(L[[i]], mods = mods_for(sc, err_mod(est$uses_ess[i])),
                       optimizer = "optim")
      m_var <- if (sc %in% c(1L, 3L)) {
        fit_rma(L[[i]], mods = mods_for(sc, var_mod(est$uses_ess[i])), optimizer = "optim")
      } else m_err
      list(err = m_err, var = m_var)
    }
  })

  pull <- function(which, field) {
    vapply(corr, function(z) {
      m <- z[[which]]
      if (is.null(m)) NA_real_
      else as.numeric(if (field == "beta") m$beta[1] else m[[field]][1])
    }, numeric(1))
  }

  est <- est |>
    dplyr::mutate(
      spec = if (legacy) "legacy (S2_v2.R as written)" else "corrected specification",
      composite_moderator_used = use_composite,
      beta0_c     = pull("err", "beta"),
      se_beta0_c  = pull("err", "se"),
      beta0_c2    = pull("var", "beta"),
      se_beta0_c2 = pull("var", "se"),
      ci_lb_c2    = pull("var", "ci.lb"),
      ci_ub_c2    = pull("var", "ci.ub"),
      gate_selects_c2 = abs(beta0) - abs(beta0_c2) > 0,
      beta0_c3 = dplyr::if_else(gate_selects_c2, beta0_c2, beta0),
      shrinkage = abs(beta0) - abs(beta0_c3),
      # Optimistic assumed effect: the t-based confidence limit FARTHER from
      # zero -- the upper limit for a positive beta0, the lower limit for a
      # negative one. Not "the upper 95% CI". Verified to preserve sign in all
      # 48 meta-analyses.
      mu_optimistic = dplyr::if_else(beta0 >= 0, ci_ub, ci_lb),
      ci_includes_zero = ci_lb < 0 & ci_ub > 0
    )

  n_failed <- sum(is.na(est$beta0_c2))
  if (n_failed > 0) {
    warning(sprintf("%d variance-moderator model(s) did not converge: %s",
                    n_failed, paste(est$case[is.na(est$beta0_c2)], collapse = ", ")))
  }
  est
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  dat <- load_datasets()
  idx <- dataset_index(dat)
  check_hierarchy(dat, idx)

  message("building estimates: corrected specification")
  est <- build_estimates(dat, idx, legacy = FALSE)
  message("building estimates: legacy specification (manuscript reproduction)")
  est_legacy <- build_estimates(dat, idx, legacy = TRUE)

  readr::write_csv(est, out_dir("estimates_per_meta_analysis.csv"))
  readr::write_csv(est_legacy, out_dir("estimates_per_meta_analysis_legacy.csv"))
  saveRDS(list(dat = dat, est = est, est_legacy = est_legacy), out_dir("estimates.rds"))

  # the one dataset the composite-moderator typo affects
  cmp <- dplyr::inner_join(
    dplyr::select(est_legacy, case, beta0, beta0_c2_legacy = beta0_c2, beta0_c3_legacy = beta0_c3),
    dplyr::select(est, case, beta0_c2, beta0_c3, k, n_study_id),
    by = "case") |>
    dplyr::filter(abs(beta0_c3_legacy - beta0_c3) > 1e-8)
  readr::write_csv(cmp, out_dir("composite_moderator_effect.csv"))

  message(sprintf("scenarios: %s",
          paste(sprintf("s%d=%d", 1:4, tabulate(est$scenario, 4)), collapse = " ")))
  message(sprintf("gate selects beta0_c2 for %d of %d; reverts to beta0 for %d",
          sum(est$gate_selects_c2), nrow(est), sum(!est$gate_selects_c2)))
  message(sprintf("CIs including zero: %d of %d", sum(est$ci_includes_zero), nrow(est)))
  message(sprintf("datasets affected by the composite-moderator typo: %d (%s)",
          nrow(cmp), paste(cmp$case, collapse = ", ")))
}
