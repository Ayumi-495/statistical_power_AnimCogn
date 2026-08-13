# R/revision/07_influence_loo.R -------------------------------------------------
# Step 7: leave-one-out influence on the meta-analysis-level summaries.
#
# WHY ALL 48 AND NOT JUST MA09. The finding that prompted this is that dropping
# `MA09` moves the k-weighted FE + VCV summary from 0.479 to 0.693. Reporting that
# one number alone would invite the obvious question - why single out that model? -
# and the honest answer is that we looked at it because it is the largest. Computing
# the full leave-one-out puts MA09 where it belongs: as the largest of 48 influence
# values, on a scale the reader can see, rather than as a special case.
#
# This is an INFLUENCE DIAGNOSTIC, not a robustness check and not an argument for
# exclusion. Nothing is dropped from any reported analysis. Rows carry
# role = "influence_check".
#
# Scope: meta-analysis level only. The primary-study-level summaries aggregate 5,740
# units with a study random effect and have no comparable single-unit leverage; the
# k-weighted 48-unit aggregate is the quantity with the problem.

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 07: leave-one-out influence, meta-analysis level ==")
o  <- readRDS(file.path(REV_TMP, "original_estimates.rds"))$original
BR <- readRDS(file.path(REV_TMP, "bias_robust.rds"))
stopifnot(identical(o$MA_model, BR$MA_model), nrow(o) == 48L)

METRICS <- c("power", "type_M", "type_S")
metric_fun <- list(power = power_two_tailed_cf, type_M = type_M_cf, type_S = type_S_cf)

# The four reported specifications, each with the standard error it is reported with.
# The critical value is z = 1.96 throughout, as decided; see 00_revision_functions.R.
specs <- tibble::tribble(
  ~effect_estimator,          ~mu,                ~se,                          ~se_source,  ~se_method,          ~role_of_spec,
  "uncorrected_beta0",        o$beta0,            o$se_beta0,                   "se_beta0",  NA_character_,       "reference_uncorrected",
  "yang2023_gated_beta0_c3",  o$beta0_c3,         o$se_beta0,                   "se_beta0",  NA_character_,       "primary",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, BR$FE_VCV_CRVE_SE_CR2,        "own_CRVE",  "CR2_Satterthwaite", "reported_sensitivity",
  "yang2024_UWLS",            BR$UWLS_estimate,   BR$UWLS_CRVE_SE_CR2_naive_t,  "own_CRVE",  "CR2_naive_t",       "supplementary"
)

# k-weighted geometric mean, matching aggregate_ma()'s point estimate. Computed here
# as a weighted mean of logs rather than via lm(), which is a second route to the same
# quantity; agreement with the canonical tables is asserted below.
wgeo <- function(v, w, off) exp(sum(w * log(v + off)) / sum(w)) - off

loo <- purrr::list_rbind(lapply(seq_len(nrow(specs)), function(i) {
  sp <- specs[i, ]
  purrr::list_rbind(lapply(METRICS, function(mt) {
    off  <- offset_for(mt)
    v    <- metric_fun[[mt]](sp$mu[[1]], sp$se[[1]])
    base <- wgeo(v, o$k, off)
    vals <- vapply(seq_len(48), function(j) wgeo(v[-j], o$k[-j], off), numeric(1))
    tibble::tibble(
      effect_estimator = sp$effect_estimator, se_source = sp$se_source,
      se_method = sp$se_method, crit_value_method = "z_1.96",
      spec_role = sp$role_of_spec, metric = mt,
      dropped_MA_model = o$MA_model, dropped_k = o$k,
      dropped_pct_of_k = 100 * o$k / sum(o$k),
      summary_all_48 = base, summary_without = vals,
      abs_change = vals - base, pct_change = 100 * (vals / base - 1)
    ) |>
      dplyr::mutate(influence_rank = rank(-abs(pct_change), ties.method = "min"))
  }))
})) |>
  dplyr::mutate(aggregation = "meta_analysis_level", weighting = "k_effect_sizes",
                role = "influence_check", verification_status = "two_derivations",
                .before = 1)

# --- verification gate: the baselines must match the canonical table ----------
S  <- readRDS(file.path(REV_TMP, "summaries.rds"))
canon <- dplyr::filter(S$ma, aggregation == "meta_analysis_level",
                       weighting == "k_effect_sizes",
                       role %in% c("reference_uncorrected", "primary",
                                   "reported_sensitivity", "supplementary"))
chk <- loo |>
  dplyr::distinct(effect_estimator, se_method, metric, summary_all_48) |>
  dplyr::inner_join(dplyr::select(canon, effect_estimator, se_method, metric, geometric_mean),
                    by = c("effect_estimator", "se_method", "metric"))
stopifnot(nrow(chk) == 12L)
d <- max(abs(chk$summary_all_48 - chk$geometric_mean))
message(sprintf("baseline agreement with meta_analysis_level_sensitivity.csv (lm vs weighted mean of logs): max|diff| = %.3g", d))
if (d > 1e-10) stop("Leave-one-out baselines disagree with the canonical summaries.")

write_revision(loo, "loo_influence.csv")

# --- what it shows -----------------------------------------------------------
pw <- dplyr::filter(loo, metric == "power")
for (r in c("primary", "reported_sensitivity")) {
  x <- dplyr::filter(pw, spec_role == r) |> dplyr::arrange(dplyr::desc(abs(pct_change)))
  message(sprintf("\n%s (%s): summary %.5f", r, x$effect_estimator[1], x$summary_all_48[1]))
  for (j in 1:3)
    message(sprintf("   drop %-11s (k=%4d, %4.1f%% of k) -> %.5f  (%+.1f%%)",
            sub("[.]csv", "", x$dropped_MA_model[j]), x$dropped_k[j],
            x$dropped_pct_of_k[j], x$summary_without[j], x$pct_change[j]))
  message(sprintf("   median |change| across all 48: %.1f%% | n above 10%%: %d | above 20%%: %d",
          stats::median(abs(x$pct_change)), sum(abs(x$pct_change) > 10), sum(abs(x$pct_change) > 20)))
}
