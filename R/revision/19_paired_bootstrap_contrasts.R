# R/revision/19_paired_bootstrap_contrasts.R ------------------------------------
# Step 19: paired paper-level bootstrap contrasts between assumed-effect specifications.
#
# WHY THIS AND NOT A COMPARISON OF THE MARGINAL INTERVALS. The manuscript says power
# "declined from 82.2% to 39.0%". Both numbers are computed from the SAME 48 models with
# the SAME weights; only the assumed effect changes. That is a paired comparison, and
# asking whether two marginal 95% intervals overlap is the wrong tool for it - it ignores
# the pairing entirely and is conservative in a way that has no interpretation here.
#
# The right quantity is the difference recomputed inside each resample:
#
#     D_b = S_uncorrected,b  -  S_corrected,b
#
# with both terms built from the SAME resampled set of source papers. This script
# computes that for every pair of the four specifications and all three metrics.
#
# THE RESAMPLES ARE THE SAME ONES `18_ma_level_uncertainty.R` USES. Same seed, same
# construction, drawn before anything else, so the marginal intervals recomputed here
# reproduce 18's output exactly - which is asserted below. Without that, the pairing
# would be meaningless.
#
# ESTIMAND. Unchanged and stated for the record: the meta-analysis-level aggregate is a
# descriptive summary of the 48 models actually included, weighted by the number of
# effect-size estimates each contributes. The bootstrap here quantifies how much that
# descriptive contrast depends on which source papers happen to be in the corpus. It is
# not an inferential claim about a wider population of meta-analyses.
#
# WHAT IS REPORTED, AND ON WHICH SCALE.
#   diff_metric   S_a - S_b on the metric's own scale (probability, or ratio for Type M).
#                 Well defined for all three metrics; this is the headline.
#   diff_log      the difference of the two fitted log-scale intercepts, i.e. the scale
#                 the estimator is linear in. exp() of it is a ratio of the aggregates.
#                 For Type S the log scale carries the 0.025 offset, so that ratio is a
#                 ratio of (value + offset) and is flagged rather than reported as if it
#                 were a ratio of Type S errors.
#
# NOTHING DOWNSTREAM IS CHANGED. This writes results/revision/ma_level_paired_contrasts.csv
# and nothing else. Table S1 and the manuscript are untouched, as instructed.

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 19: paired paper-level bootstrap contrasts ==")
o  <- readRDS(file.path(REV_TMP, "original_estimates.rds"))$original
BR <- readRDS(file.path(REV_TMP, "bias_robust.rds"))
stopifnot(identical(o$MA_model, BR$MA_model), nrow(o) == 48L)

B    <- 20000L
SEED <- 20260815L          # identical to 18_ma_level_uncertainty.R, deliberately
METRICS <- c("power", "type_M", "type_S")
metric_fun <- list(power = power_two_tailed_cf, type_M = type_M_cf, type_S = type_S_cf)

specs <- tibble::tribble(
  ~effect_estimator,          ~label,                  ~mu,                ~se,
  "uncorrected_beta0",        "uncorrected",           o$beta0,            o$se_beta0,
  "yang2023_gated_beta0_c3",  "Yang 2023 corrected",   o$beta0_c3,         o$se_beta0,
  "yang2024_FE_VCV",          "Yang 2024 FE + VCV",    BR$FE_VCV_estimate, BR$FE_VCV_CRVE_SE_CR2,
  "yang2024_UWLS",            "Yang 2024 UWLS",        BR$UWLS_estimate,   BR$UWLS_CRVE_SE_CR2_naive_t
)
spec_role <- c(uncorrected_beta0 = "reference_uncorrected",
               yang2023_gated_beta0_c3 = "primary",
               yang2024_FE_VCV = "reported_sensitivity",
               yang2024_UWLS = "supplementary")

paper  <- o$source_paper
papers <- sort(unique(paper))
n_pap  <- length(papers)
k      <- o$k
stopifnot(n_pap == 28L)

# --- the resamples: same seed, same construction, drawn first -----------------
set.seed(SEED)
paper_index <- lapply(papers, function(p) which(paper == p))
M_paper <- matrix(0L, nrow = B, ncol = nrow(o))
for (b in seq_len(B)) {
  tab <- tabulate(sample.int(n_pap, n_pap, replace = TRUE), nbins = n_pap)
  for (p in which(tab > 0L)) M_paper[b, paper_index[[p]]] <- tab[p]
}
message(sprintf("B = %d resamples of %d source papers, seed = %d (shared with step 18)",
        B, n_pap, SEED))

# --- per specification and metric: the point estimate and the B replicates ----
# Both are on the log(value + offset) scale, because that is what the weighted lm
# intercept is; the metric scale is recovered by exp(.) - offset.
cell <- list()
for (i in seq_len(nrow(specs))) {
  sp <- specs[i, ]
  for (mt in METRICS) {
    off <- offset_for(mt)
    v   <- metric_fun[[mt]](sp$mu[[1]], sp$se[[1]])
    y   <- log(v + off)
    cell[[paste(sp$effect_estimator, mt)]] <- list(
      bhat = sum(k * y) / sum(k),
      rep  = as.vector((M_paper %*% (k * y)) / (M_paper %*% k)),
      off  = off)
  }
}

# gate: these point estimates must be the canonical ones
canon <- readr::read_csv(file.path(REV_OUT, "meta_analysis_level_sensitivity.csv"),
                         show_col_types = FALSE) |>
  dplyr::filter(weighting == "k_effect_sizes", role %in% spec_role) |>
  dplyr::select(effect_estimator, role, metric, canonical = geometric_mean)
for (i in seq_len(nrow(specs))) for (mt in METRICS) {
  ce <- cell[[paste(specs$effect_estimator[i], mt)]]
  want <- canon$canonical[canon$effect_estimator == specs$effect_estimator[i] &
                          canon$metric == mt &
                          canon$role == spec_role[[specs$effect_estimator[i]]]]
  stopifnot(length(want) == 1L, abs((exp(ce$bhat) - ce$off) - want) < 1e-12)
}
message("point estimates match the canonical table in all 12 cells")

# gate: the marginal percentile intervals must reproduce step 18 exactly, which is what
# guarantees these are the same resamples and therefore that the pairing is real
prev <- readr::read_csv(file.path(REV_OUT, "ma_level_uncertainty.csv"), show_col_types = FALSE) |>
  dplyr::filter(method == "bootstrap_paper_percentile")
for (i in seq_len(nrow(prev))) {
  ce <- cell[[paste(prev$effect_estimator[i], prev$metric[i])]]
  q  <- stats::quantile(ce$rep, c(0.025, 0.975), names = FALSE)
  stopifnot(abs(exp(q[1]) - ce$off - prev$ci_lower[i]) < 1e-12,
            abs(exp(q[2]) - ce$off - prev$ci_upper[i]) < 1e-12)
}
message("marginal percentile intervals reproduce step 18 exactly: the resamples are the same draws")

# --- jackknife over papers, for BCa -------------------------------------------
lopo <- readr::read_csv(file.path(REV_OUT, "leave_one_paper_out.csv"), show_col_types = FALSE) |>
  dplyr::filter(weighting == "k_effect_sizes")
jack_log <- function(est, mt) {
  j <- lopo |> dplyr::filter(effect_estimator == est, metric == mt) |>
    dplyr::arrange(dropped_source_paper)
  stopifnot(nrow(j) == n_pap)
  log(j$geometric_mean + offset_for(mt))
}

bca_limits <- function(theta_star, theta_hat, theta_jack, level = 0.95) {
  prop <- mean(theta_star < theta_hat)
  if (prop <= 0 || prop >= 1) return(c(NA_real_, NA_real_))
  z0 <- stats::qnorm(prop)
  d  <- mean(theta_jack) - theta_jack
  den <- 6 * (sum(d^2))^1.5
  if (den == 0) return(c(NA_real_, NA_real_))
  a <- sum(d^3) / den
  z <- stats::qnorm(c((1 - level) / 2, 1 - (1 - level) / 2))
  adj <- stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
  if (any(!is.finite(adj))) return(c(NA_real_, NA_real_))
  unname(stats::quantile(theta_star, adj, names = FALSE))
}

# --- every pair of specifications ---------------------------------------------
pairs <- utils::combn(specs$effect_estimator, 2, simplify = FALSE)
out <- purrr::list_rbind(lapply(pairs, function(pr) {
  a <- pr[1]; bq <- pr[2]
  purrr::list_rbind(lapply(METRICS, function(mt) {
    ca <- cell[[paste(a, mt)]]; cb <- cell[[paste(bq, mt)]]; off <- ca$off
    sa <- exp(ca$bhat) - off; sb <- exp(cb$bhat) - off

    # replicate-wise, so the pairing is preserved inside each resample
    rep_a <- exp(ca$rep) - off; rep_b <- exp(cb$rep) - off
    d_metric <- rep_a - rep_b                 # metric scale
    d_log    <- ca$rep - cb$rep               # log(value + offset) scale

    jd_log <- jack_log(a, mt) - jack_log(bq, mt)
    jd_met <- (exp(jack_log(a, mt)) - off) - (exp(jack_log(bq, mt)) - off)

    qm  <- stats::quantile(d_metric, c(0.025, 0.975), names = FALSE)
    ql  <- stats::quantile(d_log,    c(0.025, 0.975), names = FALSE)
    bm  <- bca_limits(d_metric, sa - sb, jd_met)
    bl  <- bca_limits(d_log, ca$bhat - cb$bhat, jd_log)
    obs <- sa - sb

    tibble::tibble(
      spec_a = a, spec_b = bq, metric = mt,
      point_a = sa, point_b = sb,
      diff_metric = obs,
      diff_metric_ci_lower = qm[1], diff_metric_ci_upper = qm[2],
      diff_metric_bca_lower = bm[1], diff_metric_bca_upper = bm[2],
      diff_log = ca$bhat - cb$bhat,
      diff_log_ci_lower = ql[1], diff_log_ci_upper = ql[2],
      diff_log_bca_lower = bl[1], diff_log_bca_upper = bl[2],
      ratio_offset_scale = exp(ca$bhat - cb$bhat),
      ratio_ci_lower = exp(ql[1]), ratio_ci_upper = exp(ql[2]),
      ratio_is_on_offset_scale = off > 0,
      # how often the resampled contrast takes the opposite sign to the observed one
      boot_prop_opposite_sign = if (obs >= 0) mean(d_metric < 0) else mean(d_metric > 0),
      excludes_zero_percentile = (qm[1] > 0 & qm[2] > 0) | (qm[1] < 0 & qm[2] < 0),
      excludes_zero_bca = !is.na(bm[1]) &&
        ((bm[1] > 0 & bm[2] > 0) | (bm[1] < 0 & bm[2] < 0)),
      boot_bias_metric = mean(d_metric) - obs)
  }))
})) |>
  dplyr::mutate(aggregation = "meta_analysis_level", weighting = "k_effect_sizes",
                cluster_unit = "source_paper", n_cluster = n_pap,
                crit_value_method = "z_1.96", B = B, seed = SEED,
                estimand = "finite_corpus_descriptive_summary_of_48_models",
                role = "supplementary_diagnostic", .before = 1)

stopifnot(nrow(out) == length(pairs) * length(METRICS))
write_revision(out, "ma_level_paired_contrasts.csv")

# --- what it shows ------------------------------------------------------------
lab <- setNames(specs$label, specs$effect_estimator)
message("\npaired differences, spec_a minus spec_b, on the metric's own scale")
message("percentile and BCa are both from the same 20,000 paper-level resamples\n")
for (mt in METRICS) {
  message(sprintf("--- %s ---", mt))
  d <- dplyr::filter(out, metric == mt)
  for (i in seq_len(nrow(d))) {
    fmt <- if (mt == "type_S") "%.5f" else "%.4f"
    message(sprintf(paste0("  %-22s - %-22s  ", fmt, "  percentile [", fmt, ", ", fmt,
                           "]  BCa [", fmt, ", ", fmt, "]  %s"),
            lab[[d$spec_a[i]]], lab[[d$spec_b[i]]], d$diff_metric[i],
            d$diff_metric_ci_lower[i], d$diff_metric_ci_upper[i],
            d$diff_metric_bca_lower[i], d$diff_metric_bca_upper[i],
            if (d$excludes_zero_percentile[i] && d$excludes_zero_bca[i]) "excludes 0"
            else if (d$excludes_zero_percentile[i]) "excludes 0 (percentile only)"
            else "INCLUDES 0"))
  }
}

hl <- dplyr::filter(out, spec_a == "uncorrected_beta0",
                    spec_b == "yang2023_gated_beta0_c3", metric == "power")
message(sprintf(
  "\nHEADLINE. Power, uncorrected minus Yang-2023 corrected: %.1f percentage points",
  100 * hl$diff_metric))
message(sprintf("  percentile 95%%: [%.1f, %.1f] pp | BCa: [%.1f, %.1f] pp",
        100 * hl$diff_metric_ci_lower, 100 * hl$diff_metric_ci_upper,
        100 * hl$diff_metric_bca_lower, 100 * hl$diff_metric_bca_upper))
message(sprintf("  ratio of the two aggregates %.2fx, 95%% [%.2f, %.2f]",
        hl$ratio_offset_scale, hl$ratio_ci_lower, hl$ratio_ci_upper))
message(sprintf("  resamples in which the difference takes the opposite sign: %.2f%%",
        100 * hl$boot_prop_opposite_sign))
message(sprintf("  -> %s",
        if (hl$excludes_zero_percentile && hl$excludes_zero_bca)
          "the paired contrast excludes zero under both interval types"
        else "the paired contrast INCLUDES zero: use descriptive wording, not inferential"))
