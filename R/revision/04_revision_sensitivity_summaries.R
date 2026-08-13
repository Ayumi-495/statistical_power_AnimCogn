# R/revision/04_revision_sensitivity_summaries.R -------------------------------
# Step 4: aggregate power / Type M / Type S under each retained assumed effect.
#
# PROVENANCE IS THE POINT OF THIS SCRIPT. Every row carries four columns that fix
# what the number means, because `FE+VCV with se_beta0` and `FE+VCV with its own
# CRVE SE` answer different questions and must never be conflated:
#
#   effect_estimator    which assumed effect mu was used
#   se_source           which standard error was paired with it
#   se_method           how that standard error was computed (e.g. CR2_Satterthwaite),
#                       NA when the standard error is a plain model standard error
#   crit_value_method   which critical value defines the metric (z = 1.96 throughout
#                       the reported rows; see 00_revision_functions.R)
#   role                what the row is FOR. This encodes the decisions taken, so
#                       that no reader has to infer them:
#                         primary                   the submitted Yang-2023 analysis
#                         reported_sensitivity      FE + VCV with its own CR2 SE
#                         supplementary             UWLS
#                         reference_uncorrected     no bias correction
#                         diagnostic_*              NOT for reporting; retained so a
#                                                   choice is visible, not assumed
#   aggregation         meta_analysis_level or primary_study_level
#   weighting           k (effect sizes per meta-analysis) or unweighted
#
# At the PRIMARY-STUDY level the SE is always the individual effect size's own
# sampling error, so se_source does not vary; mu is constant within a meta-analysis.
# At the META-ANALYSIS level the SE is a property of the meta-analysis, so both
# conventions are reported.

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 04: sensitivity summaries ==")
S  <- readRDS(file.path(REV_TMP, "original_estimates.rds"))
BR <- readRDS(file.path(REV_TMP, "bias_robust.rds"))
dat <- S$dat; o <- S$original
L <- all_datasets(dat)
stopifnot(identical(names(L), o$MA_model), identical(BR$MA_model, o$MA_model))

METRICS <- c("power", "type_M", "type_S")
metric_fun <- list(power = power_two_tailed_cf, type_M = type_M_cf, type_S = type_S_cf)

# --- specifications retained -------------------------------------------------
# Only assumed effects the team has decided to retain. No exploratory variants.
#
# The meta-analysis-level pairing follows the decision taken on the Yang-2024
# estimator: it is paired with ITS OWN cluster-robust standard error, not with the
# original `se_beta0`, because the two are outputs of the same fitted analysis and
# combining the new estimate with the old standard error creates a hybrid quantity
# that neither model estimates. The hybrid rows are retained below, explicitly
# labelled as diagnostics, so the size of that distinction stays visible.
#
# `se_method` records HOW the standard error was obtained. CR2_Satterthwaite is the
# specification Yang et al. (2024) themselves use, verified from their tutorial, their
# 448-model analysis script and metafor's source, and externally validated against
# their published worked example (06_validate_yang2024_reference.R). The hand-written
# CR1 sandwich is kept as a diagnostic only.
crit_satt <- stats::qt(0.975, BR$FE_VCV_CRVE_df_Satterthwaite)
z_crit    <- rep(CRIT, nrow(o))

ma_specs <- tibble::tribble(
  ~effect_estimator,          ~mu,                ~se,                          ~se_source,  ~se_method,           ~crit,      ~crit_value_method,        ~role,                            ~verification_status,
  "uncorrected_beta0",        o$beta0,            o$se_beta0,                   "se_beta0",  NA_character_,        z_crit,     "z_1.96",                  "reference_uncorrected",          "two_derivations",
  "yang2023_gated_beta0_c3",  o$beta0_c3,         o$se_beta0,                   "se_beta0",  NA_character_,        z_crit,     "z_1.96",                  "primary",                        "two_derivations",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, BR$FE_VCV_CRVE_SE_CR2,        "own_CRVE",  "CR2_Satterthwaite",  z_crit,     "z_1.96",                  "reported_sensitivity",           "two_derivations",
  "yang2024_UWLS",            BR$UWLS_estimate,   BR$UWLS_CRVE_SE_CR2_naive_t,  "own_CRVE",  "CR2_naive_t",        z_crit,     "z_1.96",                  "supplementary",                  "two_derivations",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, BR$FE_VCV_CRVE_SE_CR2,        "own_CRVE",  "CR2_Satterthwaite",  crit_satt,  "t_Satterthwaite_per_model", "diagnostic_critical_value",     "two_derivations",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, BR$FE_VCV_CRVE_SE_CR1,        "own_CRVE",  "CR1",                z_crit,     "z_1.96",                  "diagnostic_crve_variant",        "single_derivation",
  "yang2024_UWLS",            BR$UWLS_estimate,   BR$UWLS_CRVE_SE_CR1,          "own_CRVE",  "CR1",                z_crit,     "z_1.96",                  "diagnostic_crve_variant",        "single_derivation",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, o$se_beta0,                   "se_beta0",  NA_character_,        z_crit,     "z_1.96",                  "diagnostic_hybrid_not_reported", "two_derivations",
  "yang2024_UWLS",            BR$UWLS_estimate,   o$se_beta0,                   "se_beta0",  NA_character_,        z_crit,     "z_1.96",                  "diagnostic_hybrid_not_reported", "two_derivations"
)

ma_rows <- purrr::list_rbind(lapply(seq_len(nrow(ma_specs)), function(i) {
  sp <- ma_specs[i, ]
  purrr::list_rbind(lapply(METRICS, function(mt) {
    v <- metric_fun[[mt]](sp$mu[[1]], sp$se[[1]], crit = sp$crit[[1]])
    aggregate_ma(v, o$k, mt) |>
      dplyr::mutate(metric = mt, effect_estimator = sp$effect_estimator,
                    se_source = sp$se_source, se_method = sp$se_method,
                    crit_value_method = sp$crit_value_method, role = sp$role,
                    aggregation = "meta_analysis_level",
                    weighting = "k_effect_sizes", verification_status = sp$verification_status,
                    .before = 1)
  }))
}))

# primary-study level: mu from the meta-analysis, se = each effect size's own sei
#
# THE CHOICE OF CRVE CANNOT REACH THIS TABLE. Here `se` is each primary effect size's
# own sampling error - that pairing IS the design analysis, asking what precision an
# individual study had against an assumed true effect - and only `mu` comes from the
# meta-analysis. A cluster-robust standard error is a property of the meta-analytic
# mean, so substituting it here would answer a different question entirely. The
# decision to pair the Yang-2024 estimate with its own robust SE therefore applies to
# the meta-analysis-level summary ONLY, and these rows are unchanged, to the last
# digit, by the move from CR1 to CR2.
primary_specs <- tibble::tribble(
  ~effect_estimator,          ~mu,                ~role,                   ~verification_status,
  "uncorrected_beta0",        o$beta0,            "reference_uncorrected", "two_derivations",
  "yang2023_gated_beta0_c3",  o$beta0_c3,         "primary",               "two_derivations",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, "reported_sensitivity",  "two_derivations",
  "yang2024_UWLS",            BR$UWLS_estimate,   "supplementary",         "two_derivations"
)

primary_rows <- purrr::list_rbind(lapply(seq_len(nrow(primary_specs)), function(i) {
  sp <- primary_specs[i, ]; mu <- sp$mu[[1]]
  d <- purrr::list_rbind(lapply(seq_along(L), function(j) {
    x <- L[[j]]
    tibble::tibble(study_ID = as.character(x$study_ID),
                   power  = power_two_tailed_cf(mu[j], x$sei),
                   type_M = type_M_cf(mu[j], x$sei),
                   type_S = type_S_cf(mu[j], x$sei))
  }))
  stopifnot(nrow(d) == 5740L)
  purrr::list_rbind(lapply(METRICS, function(mt) {
    aggregate_primary(d[[mt]], d$study_ID, mt) |>
      dplyr::mutate(metric = mt, effect_estimator = sp$effect_estimator,
                    se_source = "own_sampling_error_per_effect_size",
                    se_method = NA_character_, crit_value_method = "z_1.96",
                    role = sp$role,
                    aggregation = "primary_study_level", weighting = "unweighted",
                    verification_status = sp$verification_status, .before = 1)
  }))
}))

# --- equally weighted meta-analysis-level summaries ---------------------------
# Requested as a secondary descriptive summary alongside the k-weighted one, because
# the k-weighted aggregate is highly sensitive to a single large meta-analysis (see
# 07_influence_loo.R and results/revision/README.md section 3).
#
# `weights = 1` gives each of the 48 meta-analyses equal say, rather than each effect
# size equal say. It answers a different question - "what is the typical meta-analysis
# like?" rather than "what is the typical effect-size estimate like?" - so it is a
# second descriptive summary, NOT a robustness check on the first. Note in particular
# that it does not simply confirm the k-weighted picture: it moves the Yang-2023 and
# FE + VCV summaries in OPPOSITE directions.
eq_specs <- dplyr::filter(ma_specs, role %in%
  c("reference_uncorrected", "primary", "reported_sensitivity", "supplementary"))

eq_rows <- purrr::list_rbind(lapply(seq_len(nrow(eq_specs)), function(i) {
  sp <- eq_specs[i, ]
  purrr::list_rbind(lapply(METRICS, function(mt) {
    v <- metric_fun[[mt]](sp$mu[[1]], sp$se[[1]], crit = sp$crit[[1]])
    aggregate_ma(v, rep(1, length(v)), mt) |>
      dplyr::mutate(metric = mt, effect_estimator = sp$effect_estimator,
                    se_source = sp$se_source, se_method = sp$se_method,
                    crit_value_method = sp$crit_value_method,
                    role = "secondary_descriptive",
                    aggregation = "meta_analysis_level",
                    weighting = "equal_per_meta_analysis",
                    verification_status = sp$verification_status, .before = 1)
  }))
}))
ma_rows <- dplyr::bind_rows(ma_rows, eq_rows)

saveRDS(list(ma = ma_rows, primary = primary_rows), file.path(REV_TMP, "summaries.rds"))
message(sprintf("meta-analysis-level rows: %d | primary-study-level rows: %d",
                nrow(ma_rows), nrow(primary_rows)))

# --- what the CRVE choice does, printed rather than left to be looked up -------
pick <- function(df, est, sem, mt, rl = NULL) {
  x <- dplyr::filter(df, effect_estimator == est, metric == mt,
                     if (is.na(sem)) is.na(se_method) else se_method == sem)
  if (!is.null(rl)) x <- dplyr::filter(x, role == rl)
  x$geometric_mean[1]
}
message(sprintf("\nMA-level power | submitted Yang-2023 %.5f | FE+VCV own CR2 %.5f (REPORTED) | FE+VCV own CR1 %.5f (diagnostic) | FE+VCV hybrid se_beta0 %.5f (diagnostic) | uncorrected %.5f",
  pick(ma_rows, "yang2023_gated_beta0_c3", NA, "power"),
  pick(ma_rows, "yang2024_FE_VCV", "CR2_Satterthwaite", "power", "reported_sensitivity"),
  pick(ma_rows, "yang2024_FE_VCV", "CR1", "power"),
  pick(ma_rows, "yang2024_FE_VCV", NA, "power", "diagnostic_hybrid_not_reported"),
  pick(ma_rows, "uncorrected_beta0", NA, "power")))
message(sprintf("MA-level power, same CR2 SE but Satterthwaite critical value instead of z: %.5f (diagnostic only)",
  pick(ma_rows, "yang2024_FE_VCV", "CR2_Satterthwaite", "power", "diagnostic_critical_value")))
message(sprintf("primary-level power | submitted %.5f | FE+VCV %.5f | UWLS %.5f | uncorrected %.5f  (invariant to the CRVE choice)",
  pick(primary_rows, "yang2023_gated_beta0_c3", NA, "power"),
  pick(primary_rows, "yang2024_FE_VCV", NA, "power"),
  pick(primary_rows, "yang2024_UWLS", NA, "power"),
  pick(primary_rows, "uncorrected_beta0", NA, "power")))

eqp <- function(est) dplyr::filter(eq_rows, effect_estimator == est, metric == "power")$geometric_mean[1]
message(sprintf("MA-level power, EQUAL weight per meta-analysis | submitted %.5f | FE+VCV own CR2 %.5f | UWLS %.5f | uncorrected %.5f",
  eqp("yang2023_gated_beta0_c3"), eqp("yang2024_FE_VCV"), eqp("yang2024_UWLS"), eqp("uncorrected_beta0")))
message(sprintf("  -> the primary/sensitivity gap is %+.5f k-weighted and %+.5f equally weighted: the two weightings move them in OPPOSITE directions",
  pick(ma_rows, "yang2024_FE_VCV", "CR2_Satterthwaite", "power", "reported_sensitivity") -
    pick(ma_rows, "yang2023_gated_beta0_c3", NA, "power"),
  eqp("yang2024_FE_VCV") - eqp("yang2023_gated_beta0_c3")))
