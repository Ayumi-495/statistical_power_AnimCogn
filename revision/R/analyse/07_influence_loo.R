# revision/R/analyse/07_influence_loo.R -------------------------------------------------
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
# 48-unit aggregate is the quantity with the problem.
#
# BOTH WEIGHTINGS, since 2026-08-17. This script previously computed the k-weighted
# summary alone, which made it the odd one out: `15_leave_one_paper_out.R` has always
# covered both, and the equally weighted summary became a reported sensitivity analysis
# on 2026-08-15. Leaving one of the two influence analyses k-weighted-only meant the
# same question - how much does one unit move the summary? - was answered with different
# coverage depending on whether the unit was a model or a paper.
#
# The distinction matters here more than anywhere else in the workflow, because the
# leverage this script measures is largely a property of the weighting rather than of
# any model: k-weighting concentrates 22.6% of the weight on MA09, and equal weighting
# by construction gives every model 1/48.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))

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

# Weighted geometric mean, matching aggregate_ma()'s point estimate. Computed here as a
# weighted mean of logs rather than via lm(), which is a second route to the same
# quantity; agreement with the canonical tables is asserted below for both weightings.
wgeo <- function(v, w, off) exp(sum(w * log(v + off)) / sum(w)) - off

# The two summaries, as weight vectors over the retained models. `k_effect_sizes`
# describes the typical effect-size estimate; `equal_per_meta_analysis` describes the
# typical meta-analytic model. They answer different questions and are not alternative
# estimates of one quantity, so both are carried through rather than reconciled.
WEIGHTINGS <- list(
  k_effect_sizes          = function(keep) o$k[keep],
  equal_per_meta_analysis = function(keep) rep(1, sum(keep))
)

loo <- purrr::list_rbind(lapply(seq_len(nrow(specs)), function(i) {
  sp <- specs[i, ]
  purrr::list_rbind(lapply(METRICS, function(mt) {
    off <- offset_for(mt)
    v   <- metric_fun[[mt]](sp$mu[[1]], sp$se[[1]])
    purrr::list_rbind(lapply(names(WEIGHTINGS), function(wn) {
      wf   <- WEIGHTINGS[[wn]]
      base <- wgeo(v, wf(rep(TRUE, 48L)), off)
      vals <- vapply(seq_len(48), function(j) {
        keep <- seq_len(48) != j
        wgeo(v[keep], wf(keep), off)
      }, numeric(1))
      tibble::tibble(
        weighting = wn,
        effect_estimator = sp$effect_estimator, se_source = sp$se_source,
        se_method = sp$se_method, crit_value_method = "z_1.96",
        spec_role = sp$role_of_spec, metric = mt,
        dropped_MA_model = o$MA_model, dropped_k = o$k,
        dropped_pct_of_k = 100 * o$k / sum(o$k),
        summary_all_48 = base, summary_without = vals,
        abs_change = vals - base, pct_change = 100 * (vals / base - 1)
      ) |>
        # ranked within a weighting: the most influential model under k-weighting is not
        # the most influential model under equal weighting, and that is the finding.
        dplyr::mutate(influence_rank = rank(-abs(pct_change), ties.method = "min"))
    }))
  }))
})) |>
  dplyr::mutate(aggregation = "meta_analysis_level",
                role = "influence_check", verification_status = "two_derivations",
                .before = 1)

# --- verification gate: the baselines must match the canonical table ----------
# The equally weighted rows are filed as `secondary_descriptive` in the canonical table
# whatever specification produced them, so the join cannot key on `role`. Dropping the
# `diagnostic_*` rows leaves exactly the 24 reportable ones, on which
# (effect_estimator, se_method, metric, weighting) is unique.
S  <- readRDS(file.path(REV_TMP, "summaries.rds"))
canon <- dplyr::filter(S$ma, aggregation == "meta_analysis_level",
                       !startsWith(role, "diagnostic"),
                       weighting %in% names(WEIGHTINGS))
stopifnot(nrow(canon) == 24L,
          !any(duplicated(canon[, c("effect_estimator", "se_method", "metric", "weighting")])))
chk <- loo |>
  dplyr::distinct(effect_estimator, se_method, metric, weighting, summary_all_48) |>
  dplyr::inner_join(dplyr::select(canon, effect_estimator, se_method, metric, weighting,
                                  geometric_mean),
                    by = c("effect_estimator", "se_method", "metric", "weighting"))
stopifnot(nrow(chk) == 24L)
d <- max(abs(chk$summary_all_48 - chk$geometric_mean))
message(sprintf("baseline agreement with meta_analysis_level_sensitivity.csv over %d cells (lm vs weighted mean of logs): max|diff| = %.3g",
        nrow(chk), d))
if (d > 1e-10) stop("Leave-one-out baselines disagree with the canonical summaries.")

# every (model x specification x metric x weighting) exactly once
expected <- 48L * nrow(specs) * length(METRICS) * length(WEIGHTINGS)
stopifnot(nrow(loo) == expected)
combo <- table(loo$dropped_MA_model, loo$effect_estimator, loo$metric, loo$weighting)
if (!all(combo == 1L)) stop("incomplete or duplicated combinations in the leave-one-model-out grid")
message(sprintf("rows: %d = 48 models x %d specifications x %d metrics x %d weightings",
        nrow(loo), nrow(specs), length(METRICS), length(WEIGHTINGS)))

write_revision(loo, "loo_influence.csv")

# --- what it shows -----------------------------------------------------------
pw <- dplyr::filter(loo, metric == "power")
for (r in c("primary", "reported_sensitivity")) {
  for (wn in names(WEIGHTINGS)) {
    x <- dplyr::filter(pw, spec_role == r, weighting == wn) |>
      dplyr::arrange(dplyr::desc(abs(pct_change)))
    message(sprintf("\n%s (%s), %s: summary %.5f",
            r, x$effect_estimator[1], wn, x$summary_all_48[1]))
    for (j in 1:3)
      message(sprintf("   drop %-11s (k=%4d, %4.1f%% of k) -> %.5f  (%+.1f%%)",
              sub("[.]csv", "", x$dropped_MA_model[j]), x$dropped_k[j],
              x$dropped_pct_of_k[j], x$summary_without[j], x$pct_change[j]))
    message(sprintf("   median |change| across all 48: %.1f%% | n above 10%%: %d | above 20%%: %d",
            stats::median(abs(x$pct_change)), sum(abs(x$pct_change) > 10),
            sum(abs(x$pct_change) > 20)))
  }
}
