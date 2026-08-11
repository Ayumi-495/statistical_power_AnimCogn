# 04_sensitivity.R -----------------------------------------------------------
# The required reviewer/PI sensitivity analyses that do not need refitting:
#
#   B1  publication-bias correction sensitivity           (Reviewer 1)
#   B2  leave-one-meta-analytical-paper-out               (PI)
#   B3  optimistic assumed effect: farther-from-zero CI   (PI, Reviewer 2 3.6)
#   B4  external/conventional assumed-effect scenarios    (Reviewer 2 4.3, PI)
#
# B1, B3, and B4 differ only in the assumed effect mu, so they share one
# implementation. B2 varies which units are aggregated, not mu.
#
# All results use the corrected specification (see 01_estimates.R) and the
# model-based median with its 95% CI. The legacy lognormal mean is not reported
# anywhere here: it has no upper bound and returns 1.137 for meta-analysis-level
# power, which is impossible for a probability.

source(here::here("R", "02_metrics.R"))

e <- readRDS(out_dir("estimates.rds"))
dat <- e$dat; est <- e$est

# B1: correction sensitivity --------------------------------------------------
# Four assumed effects already present in the pipeline, so no new correction
# specification is invented:
#   beta0     uncorrected pooled mean
#   beta0_c   intercept of the sampling-ERROR moderator model
#   beta0_c2  intercept of the sampling-VARIANCE moderator model
#   beta0_c3  the reported corrected mean: beta0_c2 when |beta0| > |beta0_c2|,
#             otherwise beta0 (the one-directional gate)
message("B1: bias-correction sensitivity (4 proxies)")

B1_PROXIES <- c(uncorrected = "beta0", corrected_error_model = "beta0_c",
                corrected_variance_model = "beta0_c2", corrected_reported = "beta0_c3")

b1 <- purrr::list_rbind(lapply(names(B1_PROXIES), function(lab) {
  col <- B1_PROXIES[[lab]]
  summarise_all(ma = ma_metrics(est, col, lab),
                primary = primary_metrics(dat, est, col, lab))
})) |> dplyr::mutate(analysis = "B1 correction sensitivity", .before = 1)

readr::write_csv(b1, out_dir("B1_correction_sensitivity.csv"))

# Is the corrected result driven by a few strongly corrected meta-analyses?
# Reported as the distribution of shrinkage, not just the pooled change.
b1_units <- est |>
  dplyr::mutate(
    shrinkage_pct = 100 * (1 - abs(beta0_c3) / abs(beta0)),
    corrected_toward_zero = shrinkage > 0
  ) |>
  dplyr::select(case, es_type, paper, k, n_study_id, scenario,
                has_small_study_effect, has_decline_effect, gate_selects_c2,
                beta0, beta0_c, beta0_c2, beta0_c3, shrinkage, shrinkage_pct,
                corrected_toward_zero, se_beta0, se_beta0_c2) |>
  dplyr::arrange(dplyr::desc(shrinkage_pct))
readr::write_csv(b1_units, out_dir("B1_shrinkage_by_meta_analysis.csv"))

b1_concentration <- tibble::tibble(
  n_models = nrow(est),
  n_corrected_toward_zero = sum(b1_units$corrected_toward_zero),
  n_unchanged = sum(!b1_units$corrected_toward_zero),
  n_shrinkage_over_50pct = sum(b1_units$shrinkage_pct > 50, na.rm = TRUE),
  n_shrinkage_over_90pct = sum(b1_units$shrinkage_pct > 90, na.rm = TRUE),
  median_shrinkage_pct = stats::median(b1_units$shrinkage_pct, na.rm = TRUE),
  # share of the total k-weighted shrinkage contributed by the top 5 models
  top5_share_of_kweighted_shrinkage =
    sum((b1_units$shrinkage * b1_units$k)[seq_len(5)], na.rm = TRUE) /
    sum(b1_units$shrinkage * b1_units$k, na.rm = TRUE),
  # SE inflation of the variance-moderator intercept, relevant to D7
  median_se_ratio_c2_over_beta0 = stats::median(b1_units$se_beta0_c2 / b1_units$se_beta0, na.rm = TRUE),
  max_se_ratio_c2_over_beta0 = max(b1_units$se_beta0_c2 / b1_units$se_beta0, na.rm = TRUE)
)
readr::write_csv(b1_concentration, out_dir("B1_shrinkage_concentration.csv"))

# B3 and B4: assumed-effect scenarios (primary-study level) -------------------
# B3 uses the t-based confidence limit farther from zero, read from the fitted
# object. It is an intentionally optimistic scenario, NOT a better estimate of
# the underlying effect.
#
# B4 uses externally specified values, per metric, and they are never pooled
# across metrics -- SMD, Fisher's z, and lnRR share no common scale.
#   SMD  Cohen's d = 0.2 / 0.5 / 0.8, used directly.
#   Zr   r = 0.1 / 0.3 / 0.5 converted with atanh(): 0.100 / 0.310 / 0.549.
#        Cohen's d thresholds must never be applied to a z scale.
#   lnRR no conventional benchmark exists and none is invented. Interpretable
#        response ratios are used instead: 10% / 25% / 50% proportional change,
#        i.e. log(1.10) / log(1.25) / log(1.50) = 0.095 / 0.223 / 0.405. These
#        are reference scenarios, not conventions, and must be labelled as such.
message("B3: optimistic assumed effect (farther-from-zero CI limit)")

b3 <- summarise_all(primary = primary_metrics(dat, est, "mu_optimistic", "optimistic_CI_limit"),
                    ma = ma_metrics(est, "mu_optimistic", "optimistic_CI_limit")) |>
  dplyr::mutate(analysis = "B3 optimistic assumed effect", .before = 1)

# scenario aggressiveness, to be reported alongside the result
b3_scale <- tibble::tibble(
  n_models = nrow(est),
  median_ratio_mu_opt_to_beta0 = stats::median(abs(est$mu_optimistic / est$beta0)),
  min_ratio = min(abs(est$mu_optimistic / est$beta0)),
  max_ratio = max(abs(est$mu_optimistic / est$beta0)),
  n_ci_includes_zero = sum(est$ci_includes_zero),
  sign_preserved_all = all(sign(est$mu_optimistic) == sign(est$beta0)),
  ddf_min = min(est$ddf), ddf_max = max(est$ddf),
  max_abs_diff_t_vs_normal_ci =
    max(abs(est$ci_lb - (est$beta0 - stats::qnorm(0.975) * est$se_beta0)))
)
readr::write_csv(b3_scale, out_dir("B3_scenario_scale.csv"))

# per-row change distribution and threshold counts
pr_base <- primary_metrics(dat, est, "beta0", "baseline")
pr_opt  <- primary_metrics(dat, est, "mu_optimistic", "optimistic_CI_limit")
stopifnot(identical(pr_base$case, pr_opt$case), identical(pr_base$sei, pr_opt$sei))

b3_rows <- tibble::tibble(
  case = pr_base$case, es_type = pr_base$es_type, study_ID = pr_base$study_ID,
  sei = pr_base$sei,
  power_base = pr_base$power, power_opt = pr_opt$power,
  type_M_base = pr_base$type_M, type_M_opt = pr_opt$type_M,
  type_S_base = pr_base$type_S, type_S_opt = pr_opt$type_S
) |>
  dplyr::mutate(power_change_pp = 100 * (power_opt - power_base))
readr::write_csv(b3_rows, out_dir("B3_per_effect_size.csv"))

b3_dist <- tibble::tibble(
  quantile = c("min", "q1", "median", "q3", "max"),
  power_change_pp = as.numeric(stats::quantile(b3_rows$power_change_pp, c(0, .25, .5, .75, 1)))
)
b3_thresholds <- tibble::tibble(
  threshold = c(0.5, 0.8),
  n_baseline = c(sum(b3_rows$power_base >= 0.5), sum(b3_rows$power_base >= 0.8)),
  n_optimistic = c(sum(b3_rows$power_opt >= 0.5), sum(b3_rows$power_opt >= 0.8)),
  n_total = nrow(b3_rows)
)
readr::write_csv(b3_dist, out_dir("B3_change_distribution.csv"))
readr::write_csv(b3_thresholds, out_dir("B3_threshold_counts.csv"))

# by metric, since the three are on different scales
b3_by_metric <- purrr::list_rbind(lapply(c("lnRR", "SMD", "Zr"), function(mt) {
  base_m <- pr_base[pr_base$es_type == mt, ]; opt_m <- pr_opt[pr_opt$es_type == mt, ]
  dplyr::bind_rows(
    summarise_all(primary = dplyr::mutate(base_m, scenario_label = "baseline")),
    summarise_all(primary = dplyr::mutate(opt_m, scenario_label = "optimistic_CI_limit"))
  ) |> dplyr::mutate(es_type = mt, .before = 1)
}))
readr::write_csv(b3_by_metric, out_dir("B3_by_metric.csv"))

message("B4: external assumed-effect scenarios (per metric, never pooled)")

B4_SCENARIOS <- tibble::tribble(
  ~es_type, ~label,          ~mu,          ~kind,
  "SMD",    "small (d=0.2)",  0.2,          "Cohen convention",
  "SMD",    "medium (d=0.5)", 0.5,          "Cohen convention",
  "SMD",    "large (d=0.8)",  0.8,          "Cohen convention",
  "Zr",     "small (r=0.1)",  atanh(0.1),   "Cohen convention, atanh-converted",
  "Zr",     "medium (r=0.3)", atanh(0.3),   "Cohen convention, atanh-converted",
  "Zr",     "large (r=0.5)",  atanh(0.5),   "Cohen convention, atanh-converted",
  "lnRR",   "10% change",     log(1.10),    "reference scenario (not a convention)",
  "lnRR",   "25% change",     log(1.25),    "reference scenario (not a convention)",
  "lnRR",   "50% change",     log(1.50),    "reference scenario (not a convention)"
)
readr::write_csv(B4_SCENARIOS, out_dir("B4_scenario_definitions.csv"))

b4 <- purrr::list_rbind(lapply(seq_len(nrow(B4_SCENARIOS)), function(j) {
  s <- B4_SCENARIOS[j, ]
  est_m <- est |> dplyr::filter(es_type == s$es_type) |> dplyr::mutate(mu_fixed = s$mu)
  dat_m <- list(lnRR = list(), SMD = list(), Zr = list(),
                smd_uses_ess = logical(0))
  dat_m[[s$es_type]] <- all_datasets(dat)[est$es_type == s$es_type]
  pr <- primary_metrics(dat_m, est_m, "mu_fixed", s$label)
  summarise_all(primary = pr) |>
    dplyr::mutate(es_type = s$es_type, assumed_mu = s$mu, kind = s$kind, .before = 1)
}))
readr::write_csv(b4, out_dir("B4_external_scenarios.csv"))

# B2: leave-one-meta-analytical-paper-out ------------------------------------
# Re-aggregation only: every model fit is held fixed, and one source paper's
# models (1-4 each) are dropped from the aggregation at a time. 28 iterations.
message("B2: leave-one-meta-analytical-paper-out (28 iterations)")

bm <- readRDS(out_dir("baseline_metrics.rds"))$revision

loo_paper <- purrr::list_rbind(lapply(sort(unique(est$paper)), function(p) {
  ma_u <- bm$ma_unc[bm$ma_unc$paper != p, ]
  ma_c <- bm$ma_cor[bm$ma_cor$paper != p, ]
  pr_u <- bm$pr_unc[bm$pr_unc$paper != p, ]
  pr_c <- bm$pr_cor[bm$pr_cor$paper != p, ]
  dplyr::bind_rows(
    summarise_all(ma = ma_u, primary = pr_u),
    summarise_all(ma = ma_c, primary = pr_c)
  ) |>
    dplyr::mutate(paper_left_out = p,
                  n_models_dropped = sum(est$paper == p),
                  n_rows_dropped = sum(bm$pr_unc$paper == p), .before = 1)
}))
readr::write_csv(loo_paper, out_dir("B2_leave_one_paper_out.csv"))

# influence: deviation from the all-papers result
full <- bm$summary |> dplyr::filter(is.na(grouping) | grouping == "study") |>
  dplyr::select(level, metric, scenario_label, full = model_median)
b2_influence <- loo_paper |>
  dplyr::filter(is.na(grouping) | grouping == "study") |>
  dplyr::select(paper_left_out, n_models_dropped, n_rows_dropped,
                level, metric, scenario_label, model_median) |>
  dplyr::inner_join(full, by = c("level", "metric", "scenario_label")) |>
  dplyr::mutate(abs_change = model_median - full,
                pct_change = 100 * (model_median / full - 1)) |>
  dplyr::arrange(level, metric, scenario_label, dplyr::desc(abs(pct_change)))
readr::write_csv(b2_influence, out_dir("B2_influence.csv"))

b2_max <- b2_influence |>
  dplyr::group_by(level, metric, scenario_label) |>
  dplyr::slice_max(abs(pct_change), n = 1) |>
  dplyr::ungroup()
readr::write_csv(b2_max, out_dir("B2_most_influential_paper.csv"))

message("done. outputs written to outputs/")
cat("\n--- B3: optimistic assumed effect, primary-study level ---\n")
print(as.data.frame(dplyr::filter(b3, level == "primary-study", grouping == "study") |>
        dplyr::select(metric, scenario_label, model_median, ci_lower_floored, ci_upper)), digits = 4)
cat("\n--- B3 scenario scale ---\n"); print(as.data.frame(b3_scale), digits = 4)
cat("\n--- B2: largest single-paper influence ---\n")
print(as.data.frame(dplyr::select(b2_max, level, metric, scenario_label,
        paper_left_out, full, model_median, pct_change)), digits = 4)
