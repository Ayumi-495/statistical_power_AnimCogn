# revision/R/analyse/03_yang2024_bias_robust.R -----------------------------------------
# Step 3: the Yang et al. (2024) two-step bias-robust estimators, fitted to the same
# 48 meta-analysis models, as a SENSITIVITY ANALYSIS. Not adopted as a replacement
# for the Yang-2023-style primary correction.
#
# Yang Y, Lagisz M, Williams C, Noble DWA, Pan J, Nakagawa S (2024).
# Methods in Ecology and Evolution 15(9). doi:10.1111/2041-210X.14377
#
# Implementation detail is in 00_revision_functions.R (`build_vcv`,
# `bias_robust_fit`) and summarised in revision/results/README.md. In brief:
#   step one : FE + VCV, a fixed-effect GLS intercept with a within-study sampling
#              variance-covariance matrix, and UWLS, the same estimator with a
#              diagonal VCV (rho = 0)
#   clustering unit : study_ID
#   rho      : 0.5 by default, with a sensitivity analysis over {0, 0.25, 0.5, 0.75}
#   step two : cluster-robust variance estimation
#
# STEP TWO USES THE SOURCE'S OWN SPECIFICATION: CR2 with Satterthwaite degrees of
# freedom, via `metafor::robust(..., clubSandwich = TRUE)`, exactly as the authors'
# tutorial and their 448-model re-analysis call it. Verified from the primary sources
# and from metafor's own source; see the long comment block in
# 00_revision_functions.R above `fit_fe_vcv_cr2()`, and the external validation in
# 06_validate_yang2024_reference.R.
#
# The hand-written CR0 / CR1 sandwich is RETAINED AS A DIAGNOSTIC ONLY, so that the
# effect of the small-sample correction is visible rather than assumed. It is not the
# reported specification and its columns are labelled accordingly.
#
# The paper's supported claim is that the approach "does not rely on extrapolation"
# (Section 5.1). It does NOT claim the estimator cannot cross zero, and it can:
# reversal counts are written to reversal_counts.csv.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))

message("== 03: Yang et al. (2024) bias-robust estimators ==")
S <- readRDS(file.path(REV_TMP, "original_estimates.rds"))
dat <- S$dat; o <- S$original
L <- all_datasets(dat)
stopifnot(identical(names(L), o$MA_model))

fit_all <- function(rho) {
  purrr::list_rbind(lapply(seq_along(L), function(i) {
    d <- L[[i]]; cl <- as.character(d$study_ID)
    f <- bias_robust_fit(d$es, d$var, cl, rho = rho, diagonal = FALSE)
    u <- bias_robust_fit(d$es, d$var, cl, rho = rho, diagonal = TRUE)
    # the reported specification
    g  <- fit_fe_vcv_cr2(d$es, d$var, cl, rho = rho)
    uc <- fit_uwls_cr2(d$es, d$var, cl)
    tibble::tibble(
      MA_model = o$MA_model[i], rho = rho,
      # --- point estimates. Canonical value is the closed form; the rma.mv fit is
      # the second, independent derivation and the two must agree to ~1e-10.
      FE_VCV_estimate = if (is.null(f)) NA_real_ else f$beta,
      FE_VCV_estimate_rma = if (is.null(g)) NA_real_ else g$beta,
      # --- reported CRVE: CR2, Satterthwaite df
      FE_VCV_CRVE_SE_CR2 = if (is.null(g)) NA_real_ else g$se_cr2,
      FE_VCV_CRVE_df_Satterthwaite = if (is.null(g)) NA_real_ else g$df_satt,
      FE_VCV_CR2_ci_lb = if (is.null(g)) NA_real_ else g$ci_lb,
      FE_VCV_CR2_ci_ub = if (is.null(g)) NA_real_ else g$ci_ub,
      FE_VCV_CR2_pval  = if (is.null(g)) NA_real_ else g$pval,
      FE_VCV_working_SE = if (is.null(g)) NA_real_ else g$se_working,
      # --- diagnostic CRVE: hand-written sandwich, t on J-1 df. NOT reported.
      FE_VCV_CRVE_SE_CR1 = if (is.null(f)) NA_real_ else f$se_cr1,
      FE_VCV_CRVE_SE_CR0 = if (is.null(f)) NA_real_ else f$se_cr0,
      # --- UWLS, supplementary, following the source's own UWLS call
      # (CR2 with naive-t df and an explicit `target`), not the FE + VCV call.
      UWLS_estimate = u$beta,
      UWLS_estimate_lm = if (is.null(uc)) NA_real_ else uc$beta,
      UWLS_CRVE_SE_CR2_naive_t = if (is.null(uc)) NA_real_ else uc$se_cr2,
      UWLS_CR2_df = if (is.null(uc)) NA_real_ else uc$df,
      UWLS_CR2_ci_lb = if (is.null(uc)) NA_real_ else uc$ci_lb,
      UWLS_CR2_ci_ub = if (is.null(uc)) NA_real_ else uc$ci_ub,
      UWLS_CRVE_SE_CR1 = u$se_cr1, UWLS_CRVE_SE_CR0 = u$se_cr0,
      n_cluster = u$n_cluster, crve_df = u$df,
      n_negative_weight = if (is.null(f)) NA_integer_ else f$n_negative_weight,
      prop_negative_weight = if (is.null(f)) NA_real_ else f$prop_negative_weight,
      observed_effect_min = u$y_min, observed_effect_max = u$y_max,
      # --- per-model status. No silent drops: if CR2 fails for any meta-analysis the
      # 48-model aggregation denominator would change, so it is recorded and counted.
      cr2_status = if (is.null(g)) "failed" else "ok",
      uwls_cr2_status = if (is.null(uc)) "failed" else "ok",
      vcv_max_abs_diff = if (is.null(g)) NA_real_ else g$vcv_max_abs_diff,
      wrapper_max_abs_diff = if (is.null(g)) NA_real_ else g$wrapper_max_abs_diff
    )
  }))
}

main <- fit_all(RHO_DEFAULT) |>
  dplyr::left_join(dplyr::select(o, MA_model, k, beta0), by = "MA_model") |>
  dplyr::mutate(
    reversal_FE_VCV = sign(FE_VCV_estimate) != sign(beta0),
    reversal_UWLS   = sign(UWLS_estimate)   != sign(beta0),
    # a weighted mean cannot lie outside the observed effect sizes; recorded because
    # it is the property that "does not rely on extrapolation" actually delivers.
    FE_VCV_within_observed_range = FE_VCV_estimate >= observed_effect_min & FE_VCV_estimate <= observed_effect_max,
    UWLS_within_observed_range   = UWLS_estimate   >= observed_effect_min & UWLS_estimate   <= observed_effect_max,
    observed_effects_straddle_zero = observed_effect_min < 0 & observed_effect_max > 0,
    # reported interval: CR2 with the Satterthwaite critical value, taken straight
    # from clubSandwich rather than reconstructed
    FE_VCV_ci_includes_zero = FE_VCV_CR2_ci_lb < 0 & FE_VCV_CR2_ci_ub > 0,
    UWLS_ci_includes_zero   = UWLS_CR2_ci_lb   < 0 & UWLS_CR2_ci_ub   > 0,
    # diagnostic interval: hand-written CR1 sandwich, t on J-1 df
    t_crit = stats::qt(0.975, crve_df),
    FE_VCV_ci_includes_zero_CR1 = (FE_VCV_estimate - t_crit * FE_VCV_CRVE_SE_CR1 < 0) &
                                  (FE_VCV_estimate + t_crit * FE_VCV_CRVE_SE_CR1 > 0),
    UWLS_ci_includes_zero_CR1   = (UWLS_estimate   - t_crit * UWLS_CRVE_SE_CR1   < 0) &
                                  (UWLS_estimate   + t_crit * UWLS_CRVE_SE_CR1   > 0)
  )

# --- verification gates. These stop the pipeline rather than warn. ------------
n_cr2_failed <- sum(main$cr2_status != "ok") + sum(main$uwls_cr2_status != "ok")
if (n_cr2_failed > 0L)
  stop(sprintf("CR2 failed for %d model(s): %s. The 48-model denominator would change; resolve before reporting.",
               n_cr2_failed,
               paste(main$MA_model[main$cr2_status != "ok" | main$uwls_cr2_status != "ok"],
                     collapse = ", ")))

# point estimate, two independent derivations: closed-form V^-1 1 vs metafor::rma.mv
pt_diff_fe   <- max(abs(main$FE_VCV_estimate - main$FE_VCV_estimate_rma))
pt_diff_uwls <- max(abs(main$UWLS_estimate   - main$UWLS_estimate_lm))
message(sprintf("point estimate agreement: FE+VCV closed form vs rma.mv max|diff| = %.3g | UWLS closed form vs lm max|diff| = %.3g",
                pt_diff_fe, pt_diff_uwls))
if (pt_diff_fe > 1e-8 || pt_diff_uwls > 1e-8)
  stop("The two point-estimate derivations disagree beyond 1e-8. ",
       "Primary-study-level results depend on this estimate, so stop here.")

message(sprintf("VCV construction: metafor::vcalc vs build_vcv max|diff| = %.3g | metafor::robust vs clubSandwich direct max|diff| = %.3g",
                max(main$vcv_max_abs_diff), max(main$wrapper_max_abs_diff)))
if (max(main$vcv_max_abs_diff) > 1e-10 || max(main$wrapper_max_abs_diff) > 1e-10)
  stop("A CR2 cross-check exceeded tolerance.")

saveRDS(main, file.path(REV_TMP, "bias_robust.rds"))

# rho sensitivity, retained because rho is an assumption with no external anchor
rho_grid <- c(0, 0.25, 0.5, 0.75)
rs <- purrr::list_rbind(lapply(rho_grid, function(r) {
  f <- fit_all(r) |> dplyr::left_join(dplyr::select(o, MA_model, beta0), by = "MA_model") |>
    dplyr::mutate(t_crit = stats::qt(0.975, crve_df))
  tibble::tibble(
    rho = r,
    n_meta_analysis = nrow(f),
    reversal_FE_VCV = sum(sign(f$FE_VCV_estimate) != sign(f$beta0), na.rm = TRUE),
    reversal_UWLS   = sum(sign(f$UWLS_estimate)   != sign(f$beta0), na.rm = TRUE),
    FE_VCV_ci_includes_zero = sum(f$FE_VCV_CR2_ci_lb < 0 & f$FE_VCV_CR2_ci_ub > 0, na.rm = TRUE),
    FE_VCV_ci_includes_zero_CR1_diagnostic =
      sum((f$FE_VCV_estimate - f$t_crit * f$FE_VCV_CRVE_SE_CR1 < 0) &
          (f$FE_VCV_estimate + f$t_crit * f$FE_VCV_CRVE_SE_CR1 > 0), na.rm = TRUE),
    median_abs_ratio_to_beta0 = stats::median(abs(f$FE_VCV_estimate) / abs(f$beta0), na.rm = TRUE),
    n_with_negative_weight = sum(f$n_negative_weight > 0, na.rm = TRUE),
    n_cr2_failed = sum(f$cr2_status != "ok"),
    crve_variant = "CR2_Satterthwaite", verification_status = "two_derivations"
  )
}))
write_revision(rs, "rho_sensitivity.csv")

message(sprintf("FE+VCV reversals %d of 48 | UWLS %d of 48",
        sum(main$reversal_FE_VCV), sum(main$reversal_UWLS)))
message(sprintf("within observed effect range: FE+VCV %d/48, UWLS %d/48; observed effects straddle zero in %d/48",
        sum(main$FE_VCV_within_observed_range), sum(main$UWLS_within_observed_range),
        sum(main$observed_effects_straddle_zero)))
message(sprintf("REPORTED  CRVE (CR2, Satterthwaite) CI includes zero: FE+VCV %d of 48 | UWLS %d of 48",
        sum(main$FE_VCV_ci_includes_zero), sum(main$UWLS_ci_includes_zero)))
message(sprintf("diagnostic CRVE (CR1, t on J-1)     CI includes zero: FE+VCV %d of 48 | UWLS %d of 48",
        sum(main$FE_VCV_ci_includes_zero_CR1), sum(main$UWLS_ci_includes_zero_CR1)))
message(sprintf("CR2 / CR1 standard error ratio, FE+VCV: median %.3f, range %.3f-%.3f | Satterthwaite df: median %.1f, range %.1f-%.1f",
        stats::median(main$FE_VCV_CRVE_SE_CR2 / main$FE_VCV_CRVE_SE_CR1),
        min(main$FE_VCV_CRVE_SE_CR2 / main$FE_VCV_CRVE_SE_CR1),
        max(main$FE_VCV_CRVE_SE_CR2 / main$FE_VCV_CRVE_SE_CR1),
        stats::median(main$FE_VCV_CRVE_df_Satterthwaite),
        min(main$FE_VCV_CRVE_df_Satterthwaite), max(main$FE_VCV_CRVE_df_Satterthwaite)))
