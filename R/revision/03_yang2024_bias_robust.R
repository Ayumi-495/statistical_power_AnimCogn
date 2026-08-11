# R/revision/03_yang2024_bias_robust.R -----------------------------------------
# Step 3: the Yang et al. (2024) two-step bias-robust estimators, fitted to the same
# 48 meta-analysis models, as a SENSITIVITY ANALYSIS. Not adopted as a replacement
# for the Yang-2023-style primary correction.
#
# Yang Y, Lagisz M, Williams C, Noble DWA, Pan J, Nakagawa S (2024).
# Methods in Ecology and Evolution 15(9). doi:10.1111/2041-210X.14377
#
# Implementation detail is in 00_revision_functions.R (`build_vcv`,
# `bias_robust_fit`) and summarised in results/revision/README.md. In brief:
#   step one : FE + VCV, a fixed-effect GLS intercept with a within-study sampling
#              variance-covariance matrix, and UWLS, the same estimator with a
#              diagonal VCV (rho = 0)
#   clustering unit : study_ID
#   rho      : 0.5 by default, with a sensitivity analysis over {0, 0.25, 0.5, 0.75}
#   step two : cluster-robust variance estimation, sandwich, CR0 and CR1, t on J-1 df
#
# The paper's supported claim is that the approach "does not rely on extrapolation"
# (Section 5.1). It does NOT claim the estimator cannot cross zero, and it can:
# reversal counts are written to reversal_counts.csv.

source(here::here("R", "revision", "00_revision_functions.R"))

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
    tibble::tibble(
      MA_model = o$MA_model[i], rho = rho,
      FE_VCV_estimate = if (is.null(f)) NA_real_ else f$beta,
      FE_VCV_CRVE_SE_CR1 = if (is.null(f)) NA_real_ else f$se_cr1,
      FE_VCV_CRVE_SE_CR0 = if (is.null(f)) NA_real_ else f$se_cr0,
      UWLS_estimate = u$beta,
      UWLS_CRVE_SE_CR1 = u$se_cr1, UWLS_CRVE_SE_CR0 = u$se_cr0,
      n_cluster = u$n_cluster, crve_df = u$df,
      n_negative_weight = if (is.null(f)) NA_integer_ else f$n_negative_weight,
      prop_negative_weight = if (is.null(f)) NA_real_ else f$prop_negative_weight,
      observed_effect_min = u$y_min, observed_effect_max = u$y_max
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
    t_crit = stats::qt(0.975, crve_df),
    FE_VCV_ci_includes_zero = (FE_VCV_estimate - t_crit * FE_VCV_CRVE_SE_CR1 < 0) &
                              (FE_VCV_estimate + t_crit * FE_VCV_CRVE_SE_CR1 > 0),
    UWLS_ci_includes_zero   = (UWLS_estimate   - t_crit * UWLS_CRVE_SE_CR1   < 0) &
                              (UWLS_estimate   + t_crit * UWLS_CRVE_SE_CR1   > 0)
  )

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
    FE_VCV_ci_includes_zero = sum((f$FE_VCV_estimate - f$t_crit * f$FE_VCV_CRVE_SE_CR1 < 0) &
                                  (f$FE_VCV_estimate + f$t_crit * f$FE_VCV_CRVE_SE_CR1 > 0), na.rm = TRUE),
    median_abs_ratio_to_beta0 = stats::median(abs(f$FE_VCV_estimate) / abs(f$beta0), na.rm = TRUE),
    n_with_negative_weight = sum(f$n_negative_weight > 0, na.rm = TRUE),
    crve_variant = "CR1", verification_status = "two_derivations"
  )
}))
write_revision(rs, "rho_sensitivity.csv")

message(sprintf("FE+VCV reversals %d of 48 | UWLS %d of 48",
        sum(main$reversal_FE_VCV), sum(main$reversal_UWLS)))
message(sprintf("within observed effect range: FE+VCV %d/48, UWLS %d/48; observed effects straddle zero in %d/48",
        sum(main$FE_VCV_within_observed_range), sum(main$UWLS_within_observed_range),
        sum(main$observed_effects_straddle_zero)))
message(sprintf("CRVE (CR1) CI includes zero: FE+VCV %d of 48 | UWLS %d of 48",
        sum(main$FE_VCV_ci_includes_zero), sum(main$UWLS_ci_includes_zero)))
