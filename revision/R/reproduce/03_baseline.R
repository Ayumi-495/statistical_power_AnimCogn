# 03_baseline.R --------------------------------------------------------------
# Baseline run, reproduction gate, and the effect of the two baseline defects.
#
# Three sets of results are produced:
#
#   legacy     the submitted manuscript's specification, including the
#              composite-moderator typo at S2_v2.R:544,554. Used only for the
#              reproduction gate -- if this does not match the manuscript, the
#              pipeline has drifted and nothing downstream can be trusted.
#   revision   the intended specification. The basis for every revision
#              analysis and for what should be reported.
#   summaries  all three candidate summary statistics side by side, at both
#              levels, so the reporting decision is inspectable rather than
#              asserted.
#
# Expected residual deviations in the gate, both explained in the audit:
#   - Type M is Monte Carlo; S2_v2.R set no seed, so exact agreement is not
#     expected and the tolerances reflect that.
#   - primary-study power reproduces as 17.15% against 17.20% reported. The
#     reported mean (23.1%) reproduces exactly, and 17.15/23.10 are mutually
#     consistent under a single intercept while 17.20/23.1 are not, so this
#     most likely reflects an lme4/R version difference. Regenerating Table S1
#     in the declared environment (R 4.4.2, lme4 1.1.37, metafor 4.6.0) settles
#     it.

source(here::here("revision", "R", "reproduce", "02_metrics.R"))

e <- readRDS(out_dir("estimates.rds"))
dat <- e$dat

run_level_set <- function(est, tag) {
  message("computing metrics: ", tag)
  ma_unc <- ma_metrics(est, "beta0",    "uncorrected")
  ma_cor <- ma_metrics(est, "beta0_c3", "corrected")
  pr_unc <- primary_metrics(dat, est, "beta0",    "uncorrected")
  pr_cor <- primary_metrics(dat, est, "beta0_c3", "corrected")
  stopifnot("primary-level rows lost" = nrow(pr_unc) == EXPECTED_N_ROWS,
            "primary-level rows lost" = nrow(pr_cor) == EXPECTED_N_ROWS)
  list(ma_unc = ma_unc, ma_cor = ma_cor, pr_unc = pr_unc, pr_cor = pr_cor,
       summary = purrr::list_rbind(list(
         summarise_all(ma = ma_unc, primary = pr_unc),
         summarise_all(ma = ma_cor, primary = pr_cor))) |>
         dplyr::mutate(spec = tag, .before = 1))
}

legacy   <- run_level_set(e$est_legacy, "legacy (S2_v2.R as written)")
revision <- run_level_set(e$est,        "corrected specification")

baseline <- dplyr::bind_rows(legacy$summary, revision$summary)
readr::write_csv(baseline, out_dir("baseline_summaries.csv"))
saveRDS(list(legacy = legacy, revision = revision), out_dir("baseline_metrics.rds"))

# per-meta-analysis distributions of the primary-study metrics (Figure 2 bars),
# on the corrected specification
per_ma <- purrr::list_rbind(list(revision$pr_unc, revision$pr_cor)) |>
  dplyr::group_by(scenario_label, case, es_type, paper) |>
  dplyr::summarise(dplyr::across(c(power, type_M, type_S),
                    list(min = min, q1 = ~stats::quantile(.x, .25), median = stats::median,
                         q3 = ~stats::quantile(.x, .75), max = max, mean = mean)),
                   n_rows = dplyr::n(), .groups = "drop")
readr::write_csv(per_ma, out_dir("primary_metrics_by_meta_analysis.csv"))

# reproduction gate ----------------------------------------------------------
get <- function(tbl, lvl, met, scen, col = "model_median") {
  r <- tbl[tbl$level == lvl & tbl$metric == met & tbl$scenario_label == scen, ]
  if (lvl == "primary-study") r <- r[r$grouping == "study", ]
  r[[col]][1]
}
g <- function(...) get(legacy$summary, ...)

checks <- tibble::tribble(
  ~what,                                  ~observed,                                     ~reported, ~tol,
  "MA power, uncorrected (median)",        g("meta-analysis", "power",  "uncorrected"),   0.822,     0.002,
  "MA power, uncorrected (legacy mean)",   g("meta-analysis", "power",  "uncorrected", "legacy_lognormal_mean"), 1.14, 0.01,
  "MA power, corrected (median)",          g("meta-analysis", "power",  "corrected"),     0.449,     0.010,
  "MA type M, uncorrected (median)",       g("meta-analysis", "type_M", "uncorrected"),   1.11,      0.03,
  "MA type M, corrected (median)",         g("meta-analysis", "type_M", "corrected"),     2.03,      0.20,
  "MA type S, uncorrected (median)",       g("meta-analysis", "type_S", "uncorrected"),   0.0006,    0.0005,
  "MA type S, corrected (median)",         g("meta-analysis", "type_S", "corrected"),     0.0121,    0.0030,
  "primary power, uncorrected (median)",   g("primary-study", "power",  "uncorrected"),   0.172,     0.001,
  "primary power, corrected (median)",     g("primary-study", "power",  "corrected"),     0.0906,    0.0030,
  "primary type M, uncorrected (median)",  g("primary-study", "type_M", "uncorrected"),   2.86,      0.05,
  "primary type M, corrected (median)",    g("primary-study", "type_M", "corrected"),     7.79,      0.40,
  "primary type S, uncorrected (median)",  g("primary-study", "type_S", "uncorrected"),   0.0269,    0.0010,
  "primary type S, corrected (median)",    g("primary-study", "type_S", "corrected"),     0.0985,    0.0040
) |>
  dplyr::mutate(diff = observed - reported, ok = abs(diff) <= tol)

readr::write_csv(checks, out_dir("reproduction_check.csv"))
cat("\n--- reproduction gate: legacy specification vs submitted manuscript ---\n")
print(as.data.frame(checks), digits = 4)

if (all(checks$ok)) {
  message("\nreproduction gate PASSED for all ", nrow(checks), " checks")
} else {
  warning("reproduction gate FAILED for: ",
          paste(checks$what[!checks$ok], collapse = "; "),
          "\nDo not use downstream results until this is understood.")
}

# what the corrected specification changes -----------------------------------
delta <- dplyr::inner_join(
  legacy$summary |> dplyr::filter(is.na(grouping) | grouping == "study") |>
    dplyr::select(level, metric, scenario_label, legacy = model_median),
  revision$summary |> dplyr::filter(is.na(grouping) | grouping == "study") |>
    dplyr::select(level, metric, scenario_label, revised = model_median),
  by = c("level", "metric", "scenario_label")) |>
  dplyr::mutate(abs_change = revised - legacy,
                pct_change = 100 * (revised / legacy - 1))
readr::write_csv(delta, out_dir("legacy_vs_corrected_specification.csv"))
cat("\n--- effect of fixing the composite-moderator typo (MA09 only) ---\n")
print(as.data.frame(delta), digits = 4)

# the >100% legacy mean, isolated for the reviewer response -------------------
legacy_mean <- baseline |>
  dplyr::filter(level == "meta-analysis", metric == "power") |>
  dplyr::select(spec, scenario_label, model_median, ci_lower_raw, ci_lower_floored,
                ci_upper, legacy_lognormal_mean, arith_mean_unweighted,
                arith_mean_kweighted, n_at_bound, n_units)
readr::write_csv(legacy_mean, out_dir("summary_statistic_comparison_power.csv"))
cat("\n--- candidate summaries, meta-analysis-level power ---\n")
print(as.data.frame(legacy_mean), digits = 4)
