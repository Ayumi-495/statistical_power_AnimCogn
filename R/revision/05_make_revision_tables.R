# R/revision/05_make_revision_tables.R -----------------------------------------
# Step 5: assemble the canonical, machine-readable result tables.
#
# Writes to results/revision/:
#   per_meta_analysis_estimates.csv    provenance table, one row per meta-analysis
#   reversal_counts.csv                sign-reversal counts by assumed effect
#   primary_level_sensitivity.csv      primary-study-level summaries
#   meta_analysis_level_sensitivity.csv  meta-analysis-level summaries
#   rho_sensitivity.csv                written by 03_yang2024_bias_robust.R
#
# Every table carries `verification_status`. "two_derivations" means the value was
# obtained twice by different routes (typically a metafor fit and an independent
# closed-form implementation, or two separate reviewers). "single_derivation" means
# it has been computed once and is provisional.

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 05: assembling canonical tables ==")
o  <- readRDS(file.path(REV_TMP, "original_estimates.rds"))$original
d  <- readRDS(file.path(REV_TMP, "overshoot_diagnostics.rds"))
BR <- readRDS(file.path(REV_TMP, "bias_robust.rds"))
S  <- readRDS(file.path(REV_TMP, "summaries.rds"))
stopifnot(identical(o$MA_model, d$MA_model), identical(o$MA_model, BR$MA_model))

# --- 1. provenance table -----------------------------------------------------
per_ma <- o |>
  dplyr::select(MA_model, effect_size_type, source_paper, uses_effective_n,
                k, n_study_cluster, scenario,
                beta0, se_beta0, ci_lb_beta0, ci_ub_beta0,
                beta0_se_model, beta0_var_model, se_beta0_var_model,
                ci_lb_beta0_var_model, ci_ub_beta0_var_model,
                beta0_c3, gate_selects_var_model,
                has_small_study_effect, has_decline_effect) |>
  dplyr::left_join(
    BR |> dplyr::select(MA_model, rho,
            FE_VCV_estimate, FE_VCV_CRVE_SE_CR1, FE_VCV_CRVE_SE_CR0,
            UWLS_estimate, UWLS_CRVE_SE_CR1, UWLS_CRVE_SE_CR0,
            crve_cluster_n = n_cluster, crve_df,
            n_negative_weight, prop_negative_weight,
            observed_effect_min, observed_effect_max,
            FE_VCV_within_observed_range, UWLS_within_observed_range,
            observed_effects_straddle_zero,
            FE_VCV_ci_includes_zero, UWLS_ci_includes_zero),
    by = "MA_model") |>
  dplyr::left_join(
    d |> dplyr::select(MA_model, reversal_se_model, reversal_var_model, reversal_c3,
            shrinkage_magnitude, change_in_assumed_effect,
            t_c3_fixed_se, t_c3_own_se,
            corrected_ci_includes_zero, uncorrected_ci_includes_zero,
            ma_power_c3, ma_type_M_c3, ma_type_S_c3),
    by = "MA_model") |>
  dplyr::left_join(
    BR |> dplyr::transmute(MA_model, reversal_FE_VCV, reversal_UWLS),
    by = "MA_model") |>
  dplyr::mutate(
    crve_variant = "CR1_sandwich_cluster_study_ID",
    crve_verification_status = "single_derivation",
    point_estimate_verification_status = "two_derivations"
  )

stopifnot(nrow(per_ma) == 48L, sum(per_ma$k) == 5740L)
write_revision(per_ma, "per_meta_analysis_estimates.csv")

# --- 2. reversal counts ------------------------------------------------------
# Restricted to this corpus (48 models), which the scripts fully regenerate.
# External-corpus comparisons are documented in results/revision/README.md with
# their verification status rather than tabulated here, because reproducing them
# requires downloading another project's data.
rev_counts <- tibble::tribble(
  ~effect_estimator,           ~yang_equation, ~estimator_family,  ~n_reversal, ~n_effect_size_in_reversals, ~verification_status,
  "beta0_se_model",            "Eq. 2",        "Egger / PET-like", sum(d$reversal_se_model),  sum(o$k[d$reversal_se_model]),  "two_derivations",
  "beta0_var_model",           "Eq. 3",        "PEESE-like",       sum(d$reversal_var_model), sum(o$k[d$reversal_var_model]), "two_derivations",
  "yang2023_gated_beta0_c3",   "Eq. 3 + gate", "reported",         sum(d$reversal_c3),        sum(o$k[d$reversal_c3]),        "two_derivations",
  "yang2024_FE_VCV",           "2024 step 1",  "bias-robust GLS",  sum(BR$reversal_FE_VCV),   sum(o$k[BR$reversal_FE_VCV]),   "two_derivations",
  "yang2024_UWLS",             "2024 step 1",  "bias-robust WLS",  sum(BR$reversal_UWLS),     sum(o$k[BR$reversal_UWLS]),     "two_derivations"
) |>
  dplyr::mutate(n_meta_analysis = 48L, pct_reversal = round(100 * n_reversal / 48, 1),
                rho = RHO_DEFAULT, corpus = "animal_cognition_48", .before = 1)
write_revision(rev_counts, "reversal_counts.csv")

# --- 3 and 4. sensitivity summaries -----------------------------------------
col_order <- c("aggregation", "effect_estimator", "se_source", "se_method", "weighting", "metric",
               "geometric_mean", "ci_lower", "ci_upper", "raw_median", "raw_min", "raw_max",
               "arithmetic_mean_unweighted", "arithmetic_mean_kweighted",
               "n_unit", "verification_status")
write_revision(dplyr::select(S$primary, dplyr::all_of(col_order)), "primary_level_sensitivity.csv")
write_revision(dplyr::select(S$ma,      dplyr::all_of(col_order)), "meta_analysis_level_sensitivity.csv")

message("\ncanonical tables in results/revision/:")
print(sort(list.files(REV_OUT, pattern = "[.]csv$")))
