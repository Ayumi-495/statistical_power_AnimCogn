# R/revision/04_revision_sensitivity_summaries.R -------------------------------
# Step 4: aggregate power / Type M / Type S under each retained assumed effect.
#
# PROVENANCE IS THE POINT OF THIS SCRIPT. Every row carries four columns that fix
# what the number means, because `FE+VCV with se_beta0` and `FE+VCV with its own
# CRVE SE` answer different questions and must never be conflated:
#
#   effect_estimator  which assumed effect mu was used
#   se_source         which standard error was paired with it
#   se_method         how that standard error was computed (e.g. CR1), NA when the
#                     standard error is a plain model standard error
#   aggregation       meta_analysis_level or primary_study_level
#   weighting         k (effect sizes per meta-analysis) or unweighted
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
# `se_method` records HOW the standard error was obtained. Only the CR1 sandwich is
# implemented in this repository (see 00_revision_functions.R); a CR2 / Satterthwaite
# variant was used in an earlier pass but has not been independently reproduced here,
# so it is not written to the canonical tables. It is noted in
# results/revision/README.md as a provisional alternative.
ma_specs <- tibble::tribble(
  ~effect_estimator,            ~mu,                  ~se,                    ~se_source,  ~se_method,  ~verification_status,
  "uncorrected_beta0",          o$beta0,              o$se_beta0,             "se_beta0",  NA_character_, "two_derivations",
  "yang2023_gated_beta0_c3",    o$beta0_c3,           o$se_beta0,             "se_beta0",  NA_character_, "two_derivations",
  "yang2024_FE_VCV",            BR$FE_VCV_estimate,   o$se_beta0,             "se_beta0",  NA_character_, "two_derivations",
  "yang2024_UWLS",              BR$UWLS_estimate,     o$se_beta0,             "se_beta0",  NA_character_, "two_derivations",
  "yang2024_FE_VCV",            BR$FE_VCV_estimate,   BR$FE_VCV_CRVE_SE_CR1,  "own_CRVE",  "CR1",         "single_derivation",
  "yang2024_UWLS",              BR$UWLS_estimate,     BR$UWLS_CRVE_SE_CR1,    "own_CRVE",  "CR1",         "single_derivation"
)

ma_rows <- purrr::list_rbind(lapply(seq_len(nrow(ma_specs)), function(i) {
  sp <- ma_specs[i, ]
  purrr::list_rbind(lapply(METRICS, function(mt) {
    v <- metric_fun[[mt]](sp$mu[[1]], sp$se[[1]])
    aggregate_ma(v, o$k, mt) |>
      dplyr::mutate(metric = mt, effect_estimator = sp$effect_estimator,
                    se_source = sp$se_source, se_method = sp$se_method,
                    aggregation = "meta_analysis_level",
                    weighting = "k_effect_sizes", verification_status = sp$verification_status,
                    .before = 1)
  }))
}))

# primary-study level: mu from the meta-analysis, se = each effect size's own sei
primary_specs <- tibble::tribble(
  ~effect_estimator,          ~mu,                ~verification_status,
  "uncorrected_beta0",        o$beta0,            "two_derivations",
  "yang2023_gated_beta0_c3",  o$beta0_c3,         "two_derivations",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, "two_derivations",
  "yang2024_UWLS",            BR$UWLS_estimate,   "two_derivations"
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
                    se_method = NA_character_,
                    aggregation = "primary_study_level", weighting = "unweighted",
                    verification_status = sp$verification_status, .before = 1)
  }))
}))

saveRDS(list(ma = ma_rows, primary = primary_rows), file.path(REV_TMP, "summaries.rds"))
message(sprintf("meta-analysis-level rows: %d | primary-study-level rows: %d",
                nrow(ma_rows), nrow(primary_rows)))
