# revision/R/analyse/18_ma_level_uncertainty.R ------------------------------------------
# Step 18: how should the meta-analysis-level summary's uncertainty be estimated?
#
# THE POINT ESTIMATE IS NOT IN QUESTION and is not changed anywhere here. It stays
# the k-weighted geometric mean, exp(intercept) of `lm(log(value + offset) ~ 1,
# weights = k)`, with k a DESCRIPTIVE weight: it defines the estimand (the typical
# effect-size estimate) rather than claiming efficiency. Every interval below is
# computed from that same fitted model, so all of them share one point estimate,
# and the script asserts this is exact rather than approximate.
#
# WHAT IS IN QUESTION is the interval. `confint()` on that fit treats k as an
# inverse-variance weight, i.e. it assumes Var(log value_i) is proportional to 1/k_i.
# That assumption fails here: the weighted residual variance rises with k rather than
# falling. It also ignores a second dependence entirely - 12 of the 28 source papers
# contribute more than one of the 48 models, so the 48 residuals are not independent.
#
# Three intervals are therefore computed for every cell:
#
#   model_based              confint() on the weighted lm. df = 47. What Table S1
#                            currently reports. Retained as the reference.
#   CR2_paper_cluster        cluster-robust variance, CLUSTERED BY SOURCE PAPER (28),
#                            CR2 small-sample adjustment, Satterthwaite df.
#   bootstrap_paper          nonparametric cluster bootstrap over the 28 source
#                            papers: papers resampled with replacement, each carrying
#                            all of its meta-analytic models, the k-weighted geometric
#                            mean recomputed each time. Percentile and BCa.
#
# ON CLUSTERING BY PAPER RATHER THAN BY MODEL. The dependence that matters is at the
# paper level: models from one paper share authors, a laboratory, a corpus of primary
# studies and a set of analytical conventions. Clustering by model would treat the 48
# as independent and answer a question nobody asked.
#
# ON A CR2 OBJECTION THAT SHOULD NOT BE RE-LITIGATED. clubSandwich derives the CR2
# adjustment under the working model Var proportional to 1/k - the same assumption
# this script exists to distrust. That does not invalidate it. The sandwich estimator
# is consistent whatever the working model does; CR2 is a small-sample correction to
# it, not a claim that the working model is true. This is the standard use of CR2 with
# a weighted fit.
#
# ON THE 48-MODEL BOOTSTRAP. Kept as a labelled comparison row only, and marked
# `candidate = FALSE`, because it resamples models independently and so ignores exactly
# the source-paper dependence this script was written to capture. Empirically it turns
# out to be barely narrower than the paper bootstrap (median width ratio 1.74x against
# 1.75x). That is not a reason to prefer it: it is a consequence of k-weighting, which
# concentrates the weight on a handful of large models, 16 of the 28 papers contributing
# only one model each. The paper bootstrap is right on principle whether or not the two
# happen to agree here.
#
# WHAT THIS SCRIPT FOUND, so it is not lost if only the CSV is read. Under k-weighting
# the 28 papers are nowhere near 28 independent units. MA09 alone carries 22.6% of the
# weight and the top three carry 45.1%, giving a Kish effective number of papers of 9.3.
# The CR2 Satterthwaite degrees of freedom come out at 2.69 - identical in all 12 cells,
# because for a single coefficient that quantity is a function of the design (weights
# and cluster structure) and not of the outcome, which is verified in the gates below.
# So most of the difference between the CR2 interval and the model-based one is not the
# standard error at all: the CR2 SE is 1.2-3.2x the model-based SE, but the t multiplier
# is a further 3.40/2.01 = 1.69x. The CR2 SE and the bootstrap SE broadly agree.
#
# ON TYPE S. Every meta-analysis-level Type S summary is already flagged
# `summary_dominated_by_offset`: the observed values run from 7e-08 to 0.002 and the
# log-scale model needs an offset of 0.025 that exceeds nearly all of them. Bootstrap
# replicates therefore fall below the offset and `exp(mean) - offset` goes NEGATIVE,
# which a probability cannot be. Nothing is censored here. Limits are reported on the
# log(value + offset) scale as well as back-transformed, and any back-transformed
# limit below the metric's lower bound is flagged in `limit_below_bound` for a
# decision rather than quietly clamped.
#
# NOTHING DOWNSTREAM IS CHANGED. Table S1 still carries the model-based interval.
# This script only writes `revision/results/ma_level_uncertainty.csv`.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))

message("== 18: meta-analysis-level uncertainty, three methods ==")
o  <- readRDS(file.path(REV_TMP, "original_estimates.rds"))$original
BR <- readRDS(file.path(REV_TMP, "bias_robust.rds"))
stopifnot(identical(o$MA_model, BR$MA_model), nrow(o) == 48L)

B    <- 20000L
SEED <- 20260815L
METRICS <- c("power", "type_M", "type_S")
metric_fun <- list(power = power_two_tailed_cf, type_M = type_M_cf, type_S = type_S_cf)
# lower bound each metric can take, used only to FLAG limits, never to clamp them
metric_lower <- c(power = 0, type_M = 1, type_S = 0)

specs <- tibble::tribble(
  ~effect_estimator,          ~mu,                ~se,                          ~spec_role,
  "uncorrected_beta0",        o$beta0,            o$se_beta0,                   "reference_uncorrected",
  "yang2023_gated_beta0_c3",  o$beta0_c3,         o$se_beta0,                   "primary",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, BR$FE_VCV_CRVE_SE_CR2,        "reported_sensitivity",
  "yang2024_UWLS",            BR$UWLS_estimate,   BR$UWLS_CRVE_SE_CR2_naive_t,  "supplementary"
)

paper  <- o$source_paper
papers <- sort(unique(paper))
n_pap  <- length(papers)
stopifnot(n_pap == 28L)
message(sprintf("%d source papers over %d models; models per paper: %s",
        n_pap, nrow(o), paste(sprintf("%d paper(s) x %d model(s)",
        table(table(paper)), as.integer(names(table(table(paper))))), collapse = ", ")))

# --- resamples, drawn ONCE and shared by every cell ---------------------------
# Common random numbers: the same 20,000 resamples are used for all 12 cells, so
# differences between cells reflect the data and not the draw. Each resample is stored
# as a multiplicity vector over the 48 models, which turns every bootstrap replicate
# into one matrix product.
set.seed(SEED)
paper_index <- lapply(papers, function(p) which(paper == p))
M_paper <- matrix(0L, nrow = B, ncol = nrow(o))
for (b in seq_len(B)) {
  drawn <- sample.int(n_pap, n_pap, replace = TRUE)
  tab   <- tabulate(drawn, nbins = n_pap)
  for (p in which(tab > 0L)) M_paper[b, paper_index[[p]]] <- tab[p]
}
M_model <- matrix(0L, nrow = B, ncol = nrow(o))     # comparison only
for (b in seq_len(B)) M_model[b, ] <- tabulate(sample.int(nrow(o), nrow(o), replace = TRUE),
                                               nbins = nrow(o))
distinct_papers <- mean(rowSums(vapply(paper_index, function(ix) M_paper[, ix[1]] > 0L,
                                       logical(B))))
message(sprintf("bootstrap: B = %d, seed = %d; a resample contains %.1f distinct papers on average (of %d)",
        B, SEED, distinct_papers, n_pap))

# --- how many independent units does k-weighting actually leave? ---------------
# This is the finding that explains everything below, so it is computed and reported
# rather than left implicit in the degrees of freedom.
k_paper <- tapply(o$k, paper, sum)
kish_paper <- sum(k_paper)^2 / sum(k_paper^2)
kish_model <- sum(o$k)^2 / sum(o$k^2)
message(sprintf("effective units under k-weighting: Kish %.2f of %d papers (%.2f of %d models)",
        kish_paper, n_pap, kish_model, nrow(o)))
message(sprintf("  largest paper carries %.1f%% of the weight; the top three carry %.1f%%",
        100 * max(k_paper) / sum(k_paper),
        100 * sum(rev(sort(k_paper))[1:3]) / sum(k_paper)))

# jackknife for BCa: the delete-one-PAPER values of this same statistic, which
# `15_leave_one_paper_out.R` has already computed and gated against the canonical table
lopo <- readr::read_csv(file.path(REV_OUT, "leave_one_paper_out.csv"), show_col_types = FALSE) |>
  dplyr::filter(weighting == "k_effect_sizes")

bca_limits <- function(theta_star, theta_hat, theta_jack, level = 0.95) {
  prop <- mean(theta_star < theta_hat)
  if (prop <= 0 || prop >= 1) return(c(NA_real_, NA_real_))   # z0 not defined
  z0 <- stats::qnorm(prop)
  jbar <- mean(theta_jack); d <- jbar - theta_jack
  denom <- 6 * (sum(d^2))^1.5
  if (denom == 0) return(c(NA_real_, NA_real_))
  a <- sum(d^3) / denom
  z <- stats::qnorm(c((1 - level) / 2, 1 - (1 - level) / 2))
  adj <- stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
  if (any(!is.finite(adj))) return(c(NA_real_, NA_real_))
  unname(stats::quantile(theta_star, adj, names = FALSE))
}

# --- one cell -----------------------------------------------------------------
one_cell <- function(sp_i, mt) {
  sp  <- specs[sp_i, ]
  v   <- metric_fun[[mt]](sp$mu[[1]], sp$se[[1]])
  off <- offset_for(mt)
  y   <- log(v + off)
  k   <- o$k
  lb  <- metric_lower[[mt]]

  fit  <- stats::lm(y ~ 1, weights = k)
  bhat <- unname(stats::coef(fit)[1])
  point <- exp(bhat) - off

  rows <- list()
  add <- function(method, lo_log, hi_log, se, df, n_cluster, extra = list()) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      method = method,
      log_ci_lower = lo_log, log_ci_upper = hi_log,
      ci_lower = exp(lo_log) - off, ci_upper = exp(hi_log) - off,
      se_log = se, df = df, n_cluster = n_cluster, !!!extra)
  }

  # 1. model-based, exactly what Table S1 reports
  ci <- stats::confint(fit)
  add("model_based", ci[1], ci[2],
      summary(fit)$coefficients[1, 2], stats::df.residual(fit), nrow(o))

  # 2. CR2 clustered by source paper, Satterthwaite df
  V  <- clubSandwich::vcovCR(fit, cluster = paper, type = "CR2")
  ct <- clubSandwich::coef_test(fit, vcov = V, test = "Satterthwaite")
  se_cr2 <- ct$SE[1]
  df_cr2 <- ct[[grep("^df", names(ct))[1]]][1]
  tq <- stats::qt(0.975, df_cr2)
  add("CR2_paper_cluster", bhat - tq * se_cr2, bhat + tq * se_cr2, se_cr2, df_cr2, n_pap)

  # 3. cluster bootstrap over source papers
  boot_log <- function(M) as.vector((M %*% (k * y)) / (M %*% k))
  tp <- boot_log(M_paper)
  qp <- stats::quantile(tp, c(0.025, 0.975), names = FALSE)
  add("bootstrap_paper_percentile", qp[1], qp[2], stats::sd(tp), NA_real_, n_pap,
      list(B = B, seed = SEED, boot_bias_log = mean(tp) - bhat))

  jack <- lopo |>
    dplyr::filter(effect_estimator == sp$effect_estimator, metric == mt) |>
    dplyr::arrange(dropped_source_paper)
  stopifnot(nrow(jack) == n_pap,
            max(abs(jack$summary_all_28_papers - point)) < 1e-10)  # same statistic
  bca <- bca_limits(tp, bhat, log(jack$geometric_mean + off))
  add("bootstrap_paper_bca", bca[1], bca[2], stats::sd(tp), NA_real_, n_pap,
      list(B = B, seed = SEED, boot_bias_log = mean(tp) - bhat))

  # 4. the 48-model bootstrap, comparison only
  tm <- boot_log(M_model)
  qm <- stats::quantile(tm, c(0.025, 0.975), names = FALSE)
  add("bootstrap_model_percentile", qm[1], qm[2], stats::sd(tm), NA_real_, nrow(o),
      list(B = B, seed = SEED, boot_bias_log = mean(tm) - bhat))

  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      effect_estimator = sp$effect_estimator, spec_role = sp$spec_role, metric = mt,
      point_estimate = point,
      metric_lower_bound = lb,
      limit_below_bound = ci_lower < lb | ci_upper < lb,
      ci_width_log = log_ci_upper - log_ci_lower,
      offset_used = off,
      summary_dominated_by_offset = offset_flag(v, mt, point),
      candidate = method != "bootstrap_model_percentile",
      .before = 1)
}

out <- purrr::list_rbind(lapply(seq_len(nrow(specs)), function(i)
  purrr::list_rbind(lapply(METRICS, function(mt) one_cell(i, mt))))) |>
  dplyr::mutate(aggregation = "meta_analysis_level", weighting = "k_effect_sizes",
                cluster_unit = ifelse(n_cluster == n_pap, "source_paper", "meta_analysis_model"),
                crit_value_method = "z_1.96",
                kish_effective_papers = kish_paper,
                largest_paper_weight_share = max(k_paper) / sum(k_paper),
                mean_distinct_papers_per_resample = distinct_papers, .before = 1) |>
  dplyr::relocate(B, seed, boot_bias_log, .after = dplyr::last_col())

# --- gates --------------------------------------------------------------------
stopifnot(nrow(out) == nrow(specs) * length(METRICS) * 5L)

# the point estimate must be IDENTICAL across methods within a cell, not merely close
spread <- out |>
  dplyr::group_by(effect_estimator, metric) |>
  dplyr::summarise(d = diff(range(point_estimate)), .groups = "drop")
stopifnot(max(spread$d) == 0)
message("point estimate identical across all five methods in all 12 cells (exact)")

# and it must equal the canonical value in Table S1. The canonical table also holds
# `diagnostic_*` rows for the two Yang-2024 estimators (an alternative critical value,
# a CRVE variant, and a hybrid that is deliberately not reported), so the join has to
# match on `role` as well: those diagnostics are different quantities, not duplicates.
canon <- readr::read_csv(file.path(REV_OUT, "meta_analysis_level_sensitivity.csv"),
                         show_col_types = FALSE) |>
  dplyr::filter(weighting == "k_effect_sizes", role %in% specs$spec_role) |>
  dplyr::select(effect_estimator, spec_role = role, metric, canonical = geometric_mean)
cmp <- out |>
  dplyr::distinct(effect_estimator, spec_role, metric, point_estimate) |>
  dplyr::inner_join(canon, by = c("effect_estimator", "spec_role", "metric"))
stopifnot(nrow(cmp) == 12L)
d <- max(abs(cmp$point_estimate - cmp$canonical))
message(sprintf("point estimates match meta_analysis_level_sensitivity.csv over %d cells: max|diff| = %.3g",
        nrow(cmp), d))
if (d > 1e-12) stop("point estimates do not match the canonical summaries")

# the model-based interval must reproduce Table S1's interval
mb <- out |> dplyr::filter(method == "model_based") |>
  dplyr::inner_join(readr::read_csv(file.path(REV_OUT, "meta_analysis_level_sensitivity.csv"),
                                    show_col_types = FALSE) |>
                      dplyr::filter(weighting == "k_effect_sizes", role %in% specs$spec_role) |>
                      dplyr::select(effect_estimator, spec_role = role, metric,
                                    cl = ci_lower, cu = ci_upper),
                    by = c("effect_estimator", "spec_role", "metric"))
stopifnot(nrow(mb) == 12L,
          max(abs(mb$ci_lower - mb$cl)) < 1e-12, max(abs(mb$ci_upper - mb$cu)) < 1e-12)
message("model-based intervals reproduce the canonical table exactly")

# The CR2 Satterthwaite df is identical in all 12 cells. That is expected, not a bug:
# for a single coefficient it is a function of the weights and the cluster structure
# alone. Verified directly by refitting with outcomes that carry no signal - if this
# ever stops holding, the df is picking up something it should not.
df_cr2 <- dplyr::filter(out, method == "CR2_paper_cluster")$df
stopifnot(diff(range(df_cr2)) < 1e-9)
set.seed(SEED)
df_noise <- vapply(1:3, function(i) {
  f <- stats::lm(stats::rnorm(nrow(o)) ~ 1, weights = o$k)
  ct <- clubSandwich::coef_test(f, vcov = clubSandwich::vcovCR(f, cluster = paper, type = "CR2"),
                                test = "Satterthwaite")
  ct[[grep("^df", names(ct))[1]]][1]
}, numeric(1))
stopifnot(max(abs(df_noise - df_cr2[1])) < 1e-6)
message(sprintf("CR2 Satterthwaite df = %.3f, identical in all 12 cells and under outcomes with no signal (design-only, as it should be)",
        df_cr2[1]))
message(sprintf("  t multiplier %.2f against %.2f for the model-based interval (%.2fx on width alone)",
        stats::qt(0.975, df_cr2[1]), stats::qt(0.975, nrow(o) - 1),
        stats::qt(0.975, df_cr2[1]) / stats::qt(0.975, nrow(o) - 1)))

write_revision(out, "ma_level_uncertainty.csv")

# --- what it shows ------------------------------------------------------------
pct <- function(x) sprintf("%.4g", x)
for (sp in specs$effect_estimator) {
  message(sprintf("\n%s", sp))
  for (mt in METRICS) {
    d <- dplyr::filter(out, effect_estimator == sp, metric == mt)
    message(sprintf("  %-7s point %-10s", mt, pct(d$point_estimate[1])))
    for (i in seq_len(nrow(d)))
      message(sprintf("     %-28s [%9s, %9s]  width(log) %.3f%s%s",
              d$method[i], pct(d$ci_lower[i]), pct(d$ci_upper[i]), d$ci_width_log[i],
              if (!is.na(d$df[i])) sprintf("  df %.1f", d$df[i]) else "",
              if (isTRUE(d$limit_below_bound[i])) "  <-- limit below the metric's bound" else ""))
  }
}

wid <- out |>
  dplyr::select(effect_estimator, metric, method, ci_width_log) |>
  tidyr::pivot_wider(names_from = method, values_from = ci_width_log)
message("\nratio of interval width (log scale) to the model-based interval:")
message(sprintf("  CR2 clustered by paper       : median %.2fx  (range %.2f-%.2f)",
        stats::median(wid$CR2_paper_cluster / wid$model_based),
        min(wid$CR2_paper_cluster / wid$model_based), max(wid$CR2_paper_cluster / wid$model_based)))
message(sprintf("  paper bootstrap, percentile  : median %.2fx  (range %.2f-%.2f)",
        stats::median(wid$bootstrap_paper_percentile / wid$model_based),
        min(wid$bootstrap_paper_percentile / wid$model_based), max(wid$bootstrap_paper_percentile / wid$model_based)))
message(sprintf("  48-model bootstrap (excluded): median %.2fx  (range %.2f-%.2f)",
        stats::median(wid$bootstrap_model_percentile / wid$model_based),
        min(wid$bootstrap_model_percentile / wid$model_based), max(wid$bootstrap_model_percentile / wid$model_based)))
nb <- sum(out$limit_below_bound, na.rm = TRUE)
message(sprintf("\nlimits falling below the metric's lower bound: %d of %d rows (all Type S; NOT censored here)",
        nb, nrow(out)))
