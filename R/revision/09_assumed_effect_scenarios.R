# R/revision/09_assumed_effect_scenarios.R -------------------------------------
# Step 9: how much do the results depend on the assumed underlying effect?
#
# Two scenario families, both migrated here from the earlier audit pipeline so that
# they use the same closed-form metrics as everything else in results/revision/.
# The audit versions used a seeded Monte Carlo Type M (10,000 draws, 3 decimals);
# those values must not be tabulated beside closed-form ones at three significant
# figures, which is why they are recomputed rather than copied.
#
#   OPTIMISTIC (internal upper bound). The assumed effect is replaced by the limit of
#   the 95% confidence interval of the uncorrected pooled mean that lies FARTHER FROM
#   ZERO. It is deliberately the most favourable effect the data support, not a better
#   estimate. Intervals are read from the fitted models (t-based, `qt(0.975, k-1)`),
#   never reconstructed as 1.96 * se: `ddf` runs from 3 to 1,296 here.
#
#   EXTERNAL (conventional values). The assumed effect is set from outside the corpus
#   entirely. This is what Gelman & Carlin (2014) ask for and what Reviewer 2 cites:
#   a design analysis should use an effect size external to the data being evaluated.
#   SMD: Cohen's d = 0.2 / 0.5 / 0.8. Zr: Cohen's r = 0.1 / 0.3 / 0.5 expressed on the
#   analysis scale via atanh. lnRR: 10 / 25 / 50 per cent change, which are NOT
#   conventions and are labelled as reference values throughout.
#
# THE THREE LADDERS ARE NOT COMMENSURABLE. Cohen's d = 0.5 corresponds to r = 0.243,
# i.e. Zr = 0.248, whereas the Zr ladder's "medium" rung is 0.310 (equivalent to
# d = 0.62); "large" is d = 0.8 against d = 1.27. Results are therefore reported per
# metric and never pooled across metrics.
#
# WHY THE EXTERNAL SCENARIO IS ALSO RUN AT THE META-ANALYSIS LEVEL. Reviewer 2's
# objection is that meta-analysis-level power is self-referential and carries no
# information beyond the observed effect-to-uncertainty ratio. That is correct as
# stated: power is a monotone function of |mu|/se, and with mu and se from the same
# fit the ranking is fixed. The external scenario is the analysis that breaks it, and
# it changes the interpretation of the headline number - see the message printed at
# the end of this script and results/revision/README.md.

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 09: assumed-effect scenarios ==")
S <- readRDS(file.path(REV_TMP, "original_estimates.rds"))
o <- S$original; L <- all_datasets(S$dat)
stopifnot(identical(names(L), o$MA_model), nrow(o) == 48L)

METRICS <- c("power", "type_M", "type_S")
metric_fun <- list(power = power_two_tailed_cf, type_M = type_M_cf, type_S = type_S_cf)

# --- scenario definitions -----------------------------------------------------
mu_optimistic <- ifelse(o$beta0 >= 0, o$ci_ub_beta0, o$ci_lb_beta0)
stopifnot(all(sign(mu_optimistic) == sign(o$beta0)))   # the scenario must not flip direction
ratio <- abs(mu_optimistic) / abs(o$beta0)
message(sprintf("optimistic assumed effect: median %.2fx beta0 (range %.2f-%.2f); %d of 48 intervals include zero",
        stats::median(ratio), min(ratio), max(ratio),
        sum(o$ci_lb_beta0 < 0 & o$ci_ub_beta0 > 0)))

external <- tibble::tribble(
  ~effect_size_type, ~scenario_label,  ~mu,             ~convention_status,
  "SMD",  "small (d = 0.2)",   0.2,                     "Cohen convention",
  "SMD",  "medium (d = 0.5)",  0.5,                     "Cohen convention",
  "SMD",  "large (d = 0.8)",   0.8,                     "Cohen convention",
  "Zr",   "small (r = 0.1)",   atanh(0.1),              "Cohen convention, atanh-transformed to the analysis scale",
  "Zr",   "medium (r = 0.3)",  atanh(0.3),              "Cohen convention, atanh-transformed to the analysis scale",
  "Zr",   "large (r = 0.5)",   atanh(0.5),              "Cohen convention, atanh-transformed to the analysis scale",
  "lnRR", "10% change",        log(1.10),               "reference value, not a convention",
  "lnRR", "25% change",        log(1.25),               "reference value, not a convention",
  "lnRR", "50% change",        log(1.50),               "reference value, not a convention"
)

# --- primary-study level ------------------------------------------------------
# mu is constant within a meta-analysis; se is each effect size's own sampling error.
long <- purrr::list_rbind(lapply(seq_along(L), function(i) {
  x <- L[[i]]
  tibble::tibble(MA_model = o$MA_model[i], effect_size_type = o$effect_size_type[i],
                 study_ID = as.character(x$study_ID), sei = x$sei,
                 mu_uncorrected = o$beta0[i], mu_optimistic = mu_optimistic[i])
}))
stopifnot(nrow(long) == 5740L)
# Clustering matches 04_revision_sensitivity_summaries.R: prefixed by meta-analysis,
# adopted as the primary definition on 2026-08-15.
long$cluster <- namespaced_study_id(long$MA_model, long$study_ID)

summarise_primary <- function(mu, se, cl, label, grouping, extra = list()) {
  purrr::list_rbind(lapply(METRICS, function(mt) {
    v <- metric_fun[[mt]](mu, se)
    aggregate_primary(v, cl, mt) |>
      dplyr::mutate(metric = mt, scenario = label, grouping = grouping, .before = 1)
  })) |> dplyr::mutate(!!!extra)
}

prim <- dplyr::bind_rows(
  summarise_primary(long$mu_uncorrected, long$sei, long$cluster, "uncorrected pooled mean", "all metrics",
                    list(scenario_family = "baseline", convention_status = NA_character_)),
  summarise_primary(long$mu_optimistic,  long$sei, long$cluster, "optimistic (95% CI limit farther from zero)", "all metrics",
                    list(scenario_family = "optimistic", convention_status = NA_character_)),
  purrr::list_rbind(lapply(seq_len(nrow(external)), function(j) {
    e <- external[j, ]; ix <- long$effect_size_type == e$effect_size_type
    summarise_primary(rep(e$mu, sum(ix)), long$sei[ix], long$cluster[ix],
                      e$scenario_label, e$effect_size_type,
                      list(scenario_family = "external", convention_status = e$convention_status))
  }))
) |>
  dplyr::mutate(aggregation = "primary_study_level", weighting = "unweighted",
                se_source = "own_sampling_error_per_effect_size", crit_value_method = "z_1.96",
                clustering_unit = "meta_analysis_x_study_ID", .before = 1)

# per-metric baseline and optimistic rows, so the ladder is comparable within a metric
prim_by_metric <- purrr::list_rbind(lapply(c("lnRR", "SMD", "Zr"), function(tt) {
  ix <- long$effect_size_type == tt
  dplyr::bind_rows(
    summarise_primary(long$mu_uncorrected[ix], long$sei[ix], long$cluster[ix],
                      "uncorrected pooled mean", tt, list(scenario_family = "baseline", convention_status = NA_character_)),
    summarise_primary(long$mu_optimistic[ix],  long$sei[ix], long$cluster[ix],
                      "optimistic (95% CI limit farther from zero)", tt,
                      list(scenario_family = "optimistic", convention_status = NA_character_)))
})) |>
  dplyr::mutate(aggregation = "primary_study_level", weighting = "unweighted",
                se_source = "own_sampling_error_per_effect_size", crit_value_method = "z_1.96",
                clustering_unit = "meta_analysis_x_study_ID", .before = 1)

# --- meta-analysis level ------------------------------------------------------
summarise_ma <- function(mu, se, k, label, grouping, extra = list()) {
  purrr::list_rbind(lapply(METRICS, function(mt) {
    v <- metric_fun[[mt]](mu, se)
    dplyr::bind_rows(
      aggregate_ma(v, k, mt)                  |> dplyr::mutate(weighting = "k_effect_sizes"),
      aggregate_ma(v, rep(1, length(v)), mt)  |> dplyr::mutate(weighting = "equal_per_meta_analysis")
    ) |> dplyr::mutate(metric = mt, scenario = label, grouping = grouping, .before = 1)
  })) |> dplyr::mutate(!!!extra)
}

ma <- dplyr::bind_rows(
  purrr::list_rbind(lapply(c("lnRR", "SMD", "Zr"), function(tt) {
    ix <- o$effect_size_type == tt
    summarise_ma(o$beta0[ix], o$se_beta0[ix], o$k[ix], "uncorrected pooled mean", tt,
                 list(scenario_family = "baseline", convention_status = NA_character_))
  })),
  purrr::list_rbind(lapply(seq_len(nrow(external)), function(j) {
    e <- external[j, ]; ix <- o$effect_size_type == e$effect_size_type
    summarise_ma(rep(e$mu, sum(ix)), o$se_beta0[ix], o$k[ix], e$scenario_label, e$effect_size_type,
                 list(scenario_family = "external", convention_status = e$convention_status))
  }))
) |>
  dplyr::mutate(aggregation = "meta_analysis_level", se_source = "se_beta0",
                crit_value_method = "z_1.96", clustering_unit = NA_character_, .before = 1)

out <- dplyr::bind_rows(prim, prim_by_metric, ma) |>
  dplyr::mutate(verification_status = "single_derivation_migrated_from_audit_pipeline")
write_revision(out, "assumed_effect_scenarios.csv")

# --- the optimistic scenario at the meta-analysis level is deterministic -------
# Recorded because an earlier note claimed it was "near-vacuous because power -> 1".
# That is wrong. At this level se = se_beta0, so
#     |mu_opt| / se_beta0 = |beta0| / se_beta0 + qt(0.975, ddf),
# an exact additive shift of the Wald statistic by a known constant. The scenario is
# uninformative because it is DETERMINISTIC, not because it saturates: the arithmetic
# floor is pnorm(qt(0.975, ddf) - 1.96), and the observed minimum is well below 1.
lam_base <- abs(o$beta0) / o$se_beta0
lam_opt  <- abs(mu_optimistic) / o$se_beta0
tcrit    <- stats::qt(0.975, o$k - 1)
message(sprintf("optimistic at MA level: max|lambda_opt - (lambda_base + t_crit)| = %.3g (exact additive shift)",
        max(abs(lam_opt - (lam_base + tcrit)))))
po <- power_two_tailed_cf(mu_optimistic, o$se_beta0)
message(sprintf("  power: min %.3f, median %.4f, max %.4f | arithmetic floor %.4f-%.4f | NOT 'power -> 1'",
        min(po), stats::median(po), max(po), min(stats::pnorm(tcrit - CRIT)), max(stats::pnorm(tcrit - CRIT))))

# --- what the external scenario does to the meta-analysis-level interpretation --
pw <- function(d) d$geometric_mean[d$metric == "power"]
mtxt <- sapply(c("SMD", "Zr", "lnRR"), function(tt) {
  m <- dplyr::filter(ma, grouping == tt, weighting == "k_effect_sizes", metric == "power",
                     scenario %in% c("medium (d = 0.5)", "medium (r = 0.3)", "25% change"))
  sprintf("%s %.3f", tt, m$geometric_mean[1])
})
message("\nMETA-ANALYSIS-LEVEL POWER AGAINST AN EXTERNAL *MEDIUM* EFFECT: ", paste(mtxt, collapse = " | "))
message("  against a reported corrected summary of 0.390. The models are precise; the low")
message("  corrected value reflects corrected means near zero, not imprecision. This is the")
message("  substantive answer to the self-reference objection.")

message(sprintf("\nprimary-study level, pooled: uncorrected %.5f -> optimistic %.5f (power)",
        pw(dplyr::filter(prim, scenario_family == "baseline", metric == "power")),
        pw(dplyr::filter(prim, scenario_family == "optimistic", metric == "power"))))
nf <- sum(out$summary_dominated_by_offset, na.rm = TRUE)
message(sprintf("rows where the Type S summary is dominated by the 0.025 offset and the raw quantiles must be reported instead: %d of %d",
        nf, nrow(out)))

# --- clustering-unit comparison, for the record ------------------------------
# `study_ID` is the raw identifier from each source dataset, so the same author-year
# string can occur in several meta-analyses and those rows are merged into one cluster:
# 130 identifiers covering 1,098 of the 5,740 rows. Prefixing the meta-analysis gives a
# consistent within-meta-analysis unit. Reported here as a comparison only; the decision
# was taken on 2026-08-15: the prefixed version is now the primary definition, and the
# raw figures below are kept only to document the size of the change.
n_collide <- sum(duplicated(unique(long[, c("MA_model","study_ID")])$study_ID))
message(sprintf("\nclustering comparison: %d raw identifiers occur in more than one meta-analysis",
        n_collide))
cl_raw <- long$study_ID; cl_pre <- namespaced_study_id(long$MA_model, long$study_ID)
for (lab in c("uncorrected", "optimistic")) {
  mu <- if (lab == "uncorrected") long$mu_uncorrected else long$mu_optimistic
  v <- power_two_tailed_cf(mu, long$sei)
  a <- aggregate_primary(v, cl_raw, "power")$geometric_mean
  b <- aggregate_primary(v, cl_pre, "power")$geometric_mean
  message(sprintf("  power, %-12s raw %.5f | prefixed %.5f  (raw retained only as a record of the change)", lab, a, b))
}
