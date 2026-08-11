# R/revision/02_overshoot_diagnostics.R ----------------------------------------
# Step 2: diagnostics for over-correction in the Yang-2023-style step.
#
# The question is Reviewer 1's and the PI's: does the bias correction itself drive
# the reported power / Type M / Type S results? Two mechanisms are separated here.
#
#   (a) SIGN REVERSAL. The corrected mean can come out with the opposite sign to
#       the uncorrected mean. The gate compares magnitudes only and is blind to a
#       zero crossing, so a small reversal passes it.
#   (b) DIVERGENCE AS mu -> 0. Type M tends to 2.3378/|t| as the assumed effect
#       approaches zero, so a corrected mean near zero produces an arbitrarily
#       large Type M irrespective of how it was obtained.
#
# (b) is the reason (a) matters: the SMALLEST reversals produce the LARGEST Type M.
# Both are recorded per meta-analysis so the relationship is inspectable.

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 02: overshoot diagnostics ==")
o <- readRDS(file.path(REV_TMP, "original_estimates.rds"))$original

diag <- o |>
  dplyr::mutate(
    reversal_se_model  = sign(beta0_se_model)  != sign(beta0),
    reversal_var_model = sign(beta0_var_model) != sign(beta0),
    reversal_c3        = sign(beta0_c3)        != sign(beta0),
    # shrinkage as the submitted analysis defines it: a magnitude difference, which
    # is structurally blind to a sign reversal. MA13_02 is the clearest case:
    # -0.379 -> +0.372 scores as 0.008.
    shrinkage_magnitude = abs(beta0) - abs(beta0_c3),
    change_in_assumed_effect = abs(beta0 - beta0_c3),
    # t on the pipeline's own pairing: corrected mean with the UNcorrected SE. This
    # is Yang 2023's deliberate choice and is what drives the reported Type M. It is
    # NOT a test of whether the reversed effect is real.
    t_c3_fixed_se = abs(beta0_c3) / se_beta0,
    # t on the corrected estimate's own SE, which is the inferential question.
    t_c3_own_se   = abs(beta0_c3) / se_beta0_var_model,
    corrected_ci_includes_zero = ci_lb_beta0_var_model < 0 & ci_ub_beta0_var_model > 0,
    uncorrected_ci_includes_zero = ci_lb_beta0 < 0 & ci_ub_beta0 > 0,
    ma_power_c3  = power_two_tailed_cf(beta0_c3, se_beta0),
    ma_type_M_c3 = type_M_cf(beta0_c3, se_beta0),
    ma_type_S_c3 = type_S_cf(beta0_c3, se_beta0)
  ) |>
  dplyr::select(MA_model, effect_size_type, k, n_study_cluster, scenario,
                beta0, beta0_se_model, beta0_var_model, beta0_c3,
                dplyr::starts_with("reversal_"), shrinkage_magnitude,
                change_in_assumed_effect, t_c3_fixed_se, t_c3_own_se,
                corrected_ci_includes_zero, uncorrected_ci_includes_zero,
                ma_power_c3, ma_type_M_c3, ma_type_S_c3)


saveRDS(diag, file.path(REV_TMP, "overshoot_diagnostics.rds"))

message(sprintf("reversals: SE model %d | variance model %d | reported (gated) %d  of 48",
        sum(diag$reversal_se_model), sum(diag$reversal_var_model), sum(diag$reversal_c3)))
message(sprintf("corrected CI includes zero: %d of 48 (%d effect sizes)",
        sum(diag$corrected_ci_includes_zero), sum(diag$k[diag$corrected_ci_includes_zero])))
message(sprintf("reported reversals with |t| on their OWN se >= 1.96: %d of %d",
        sum(diag$reversal_c3 & diag$t_c3_own_se >= 1.96), sum(diag$reversal_c3)))
message(sprintf("meta-analysis-level Type M above 20: %d (max %.1f)",
        sum(diag$ma_type_M_c3 > 20), max(diag$ma_type_M_c3)))
