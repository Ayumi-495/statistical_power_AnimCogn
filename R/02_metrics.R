# 02_metrics.R ---------------------------------------------------------------
# Power, Type M, and Type S at both levels, under any assumed effect (mu).
#
# The assumed effect is a parameter here rather than hard-coded, because four of
# the required revision analyses differ only in what mu is: the uncorrected mean,
# the corrected mean, the farther-from-zero confidence limit, and externally
# specified benchmark values. One implementation, several mu.
#
# Two structural facts about these estimators, established in the audit and
# worth keeping in view when reading the outputs:
#
#   All three metrics are functions of t = |mu|/se alone. Power and Type S are
#   so by inspection; in Type M, writing est = se*(t + Z) makes se cancel. So
#   the three are deterministic re-encodings of one scalar, not three
#   independent indicators.
#
#   At the meta-analysis level SE is se_beta0, the standard error of the pooled
#   estimate -- NOT the "average sampling variance across all effect sizes" the
#   manuscript Methods describes. That wording is a documentation error; see
#   docs/03_manuscript_corrections.md. Consequently the uncorrected
#   meta-analysis-level quantity is a monotone transformation of the pooled Wald
#   statistic.
#
# Corrected effects are paired with the UNcorrected se_beta0, matching S2_v2.R
# and Yang et al. That is deliberate: precision is held fixed so that only the
# assumed effect varies between the corrected and uncorrected comparison.

source(here::here("R", "00_setup.R"))

# meta-analysis level --------------------------------------------------------
# mu_col: a column of `est` holding the assumed effect for each meta-analysis.
ma_metrics <- function(est, mu_col, label, N_mc = 10000) {
  mu <- est[[mu_col]]
  se <- est$se_beta0
  # values computed before the tibble() call: tibble evaluates arguments
  # sequentially with data masking, so a column named `mu` would shadow the
  # vector `mu` in any later argument.
  pw <- power_two_tailed(mu, se)
  tM <- vapply(seq_along(mu), function(j) error_M(mu[j], se[j], N = N_mc), numeric(1))
  tS <- error_S(mu, se)
  tibble::tibble(
    case = est$case, es_type = est$es_type, paper = est$paper, k = est$k,
    scenario_label = label, mu = mu, se = se,
    power = pw, type_M = tM, type_S = tS
  )
}

# primary-study level --------------------------------------------------------
# Within one meta-analysis every row shares the same mu, so all within-
# meta-analysis variation in these metrics is a monotone function of the row's
# own sei. Each row keeps its observed sei in every scenario.
primary_metrics <- function(dat, est, mu_col, label, N_mc = 10000) {
  L <- all_datasets(dat)
  stopifnot(identical(names(L), est$case))
  mu <- est[[mu_col]]

  purrr::list_rbind(lapply(seq_along(L), function(i) {
    d <- L[[i]]
    # mu_i is bound outside the tibble() call on purpose: tibble evaluates its
    # arguments sequentially with data masking, so a column named `mu` shadows
    # the vector `mu`, and `mu[i]` would silently index the new column instead
    # (returning NA for every dataset with fewer rows than i).
    mu_i <- mu[[i]]
    se_i <- d$sei
    pw <- power_two_tailed(mu_i, se_i)
    tM <- error_M_vec(mu_i, se_i, N = N_mc)
    tS <- error_S(mu_i, se_i)
    tibble::tibble(
      case = est$case[i], es_type = est$es_type[i], paper = est$paper[i],
      study_ID = as.character(d$study_ID), sei = se_i,
      scenario_label = label, mu = mu_i,
      power = pw, type_M = tM, type_S = tS
    )
  }))
}

# summaries ------------------------------------------------------------------
# TYPE_S_OFFSET: Type S can be exactly 0, so the log-scale model is fitted on
# x + 0.025 and the offset subtracted afterwards. This is inherited from
# S2_v2.R. Note it can push the CI's lower limit below zero (the uncorrected
# meta-analysis-level Type S lower limit is -0.0007), which S2_v2.R floors to 0
# via a code comment rather than a stated method. The floor is applied here
# explicitly and flagged in the output.
TYPE_S_OFFSET <- 0.025

offset_for <- function(metric) if (metric == "type_S") TYPE_S_OFFSET else 0

# Three candidate summaries, reported side by side so the choice is inspectable.
#
#   model_median  the back-transformed log-scale intercept with its 95% CI.
#                 The principal reported summary: it is the fitted model's own
#                 estimate, it carries a valid interval, and it uses the model's
#                 own weighting.
#   legacy_mean   exp(intercept + 0.5*var(log x)), the lognormal mean used in
#                 the submitted manuscript. Reproduced for diagnosis only: it
#                 has no upper bound and returns 1.137 for meta-analysis-level
#                 power, which is impossible for a probability. Never report.
#   arith_mean    the original-scale arithmetic mean. Bounded and assumption-
#                 free but purely descriptive; at the meta-analysis level it is
#                 given both unweighted and k-weighted, because the model
#                 weights by k and the two differ materially (0.710 vs 0.886 for
#                 uncorrected power).
summarise_ma <- function(m, metric) {
  off <- offset_for(metric)
  x <- m[[metric]]
  y <- log(x + off)
  fit <- stats::lm(y ~ 1, weights = m$k)
  ci <- stats::confint(fit)
  tibble::tibble(
    level = "meta-analysis", metric = metric, scenario_label = m$scenario_label[1],
    n_units = length(x),
    model_median = exp(stats::coef(fit)[[1]]) - off,
    ci_lower_raw = exp(ci[1]) - off,
    ci_upper     = exp(ci[2]) - off,
    ci_lower_floored = pmax(ci_lower_raw, 0),
    legacy_lognormal_mean = exp(stats::coef(fit)[[1]] + 0.5 * stats::var(y)) - off,
    arith_mean_unweighted = mean(x),
    arith_mean_kweighted  = sum(x * m$k) / sum(m$k),
    min = min(x), q1 = stats::quantile(x, .25), median_raw = stats::median(x),
    q3 = stats::quantile(x, .75), max = max(x),
    n_at_bound = sum(x >= 0.999)
  )
}

# grouping: "study" is the exact Yang et al. port, (1 | study_ID). "study+model"
# adds a meta-analysis random effect and is reported as a comparison only -- it
# is not a settled decision (pending PI discussion), so it never overwrites the
# primary result.
summarise_primary <- function(m, metric, grouping = c("study", "study+model")) {
  grouping <- match.arg(grouping)
  off <- offset_for(metric)
  x <- m[[metric]]
  d <- data.frame(y = log(x + off), study_ID = m$study_ID, case = m$case)
  fit <- if (grouping == "study") {
    lme4::lmer(y ~ 1 + (1 | study_ID), data = d)
  } else {
    lme4::lmer(y ~ 1 + (1 | study_ID) + (1 | case), data = d)
  }
  ci <- stats::confint(fit, method = "Wald")
  ci <- ci[rownames(ci) == "(Intercept)", ]
  b0 <- lme4::fixef(fit)[[1]]
  tibble::tibble(
    level = "primary-study", metric = metric, scenario_label = m$scenario_label[1],
    grouping = grouping, n_units = length(x),
    model_median = exp(b0) - off,
    ci_lower_raw = exp(ci[[1]]) - off,
    ci_upper     = exp(ci[[2]]) - off,
    ci_lower_floored = pmax(exp(ci[[1]]) - off, 0),
    legacy_lognormal_mean = exp(b0 + 0.5 * stats::var(d$y)) - off,
    arith_mean_unweighted = mean(x),
    arith_mean_kweighted  = NA_real_,   # the primary-level model is unweighted
    min = min(x), q1 = stats::quantile(x, .25), median_raw = stats::median(x),
    q3 = stats::quantile(x, .75), max = max(x),
    n_at_bound = sum(x >= 0.999)
  )
}

summarise_all <- function(ma = NULL, primary = NULL, grouping = "study") {
  metrics <- c("power", "type_M", "type_S")
  out <- list()
  if (!is.null(ma))      out <- c(out, lapply(metrics, function(v) summarise_ma(ma, v)))
  if (!is.null(primary)) out <- c(out, lapply(metrics, function(v) summarise_primary(primary, v, grouping)))
  purrr::list_rbind(out)
}
