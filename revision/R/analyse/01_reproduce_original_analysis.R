# revision/R/analyse/01_reproduce_original_analysis.R -----------------------------------
# Step 1 of the revision workflow: reproduce the Yang-2023-style estimates that the
# submitted analysis is built on, so every later comparison has a fixed baseline.
#
# Produces, per meta-analysis model:
#   beta0             uncorrected pooled mean, intercept of the multilevel model
#   beta0_se_model    intercept of the sampling-ERROR moderator model   (Yang 2023 Eq. 2)
#   beta0_var_model   intercept of the sampling-VARIANCE moderator model (Yang 2023 Eq. 3)
#   beta0_c3          the REPORTED corrected mean: beta0_var_model when
#                     |beta0| > |beta0_var_model|, else beta0 (a magnitude-only gate)
#
# Scenario assignment follows Yang 2023 and `S2_v2.R`: it uses the SIGN of
# beta0*beta1 and beta0*beta2 only, with no significance test. That is the paper's
# stated method, not an error in the port; see revision/results/README.md.
#
# `S2_v2.R` is not sourced. The estimates are built by the already-committed
# revision/R/reproduce/01_estimates.R, sourced into a PRIVATE environment because its bottom-of-file
# guard fires when sourced into the global environment and would overwrite
# revision/results/reproduce/.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))

message("== 01: reproducing the Yang-2023-style estimates ==")

dat <- load_datasets()
idx <- dataset_index(dat)
check_hierarchy(dat, idx)

env <- new.env()
sys.source(here::here("revision", "R", "reproduce", "01_estimates.R"), envir = env)

est <- env$build_estimates(dat, idx, legacy = FALSE)

original <- est |>
  dplyr::transmute(
    MA_model          = case,
    effect_size_type  = es_type,
    source_paper      = paper,
    uses_effective_n  = uses_ess,
    k                 = k,
    n_study_cluster   = n_study_id,
    scenario          = scenario,
    beta0             = beta0,
    se_beta0          = se_beta0,
    ci_lb_beta0       = ci_lb,
    ci_ub_beta0       = ci_ub,
    beta0_se_model    = beta0_c,
    beta0_var_model   = beta0_c2,
    se_beta0_var_model = se_beta0_c2,
    ci_lb_beta0_var_model = ci_lb_c2,
    ci_ub_beta0_var_model = ci_ub_c2,
    beta0_c3          = beta0_c3,
    gate_selects_var_model = gate_selects_c2,
    has_small_study_effect = has_small_study_effect,
    has_decline_effect     = has_decline_effect
  )

stopifnot(nrow(original) == 48L, sum(original$k) == 5740L)
saveRDS(list(dat = dat, original = original), file.path(REV_TMP, "original_estimates.rds"))

message(sprintf("scenarios: %s", paste(sprintf("s%d=%d", 1:4, tabulate(original$scenario, 4)), collapse = " ")))
message(sprintf("gate selects the variance-model estimate for %d of 48",
                sum(original$gate_selects_var_model)))
