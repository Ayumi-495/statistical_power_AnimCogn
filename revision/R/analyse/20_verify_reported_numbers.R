# revision/R/analyse/20_verify_reported_numbers.R ---------------------------------------
# Step 20: a second derivation, by a different route, for everything that reaches the
# manuscript and was not already covered.
#
# WHAT WAS ALREADY COVERED, and is therefore not repeated here:
#   scenarios (Table S1 parts B and C)  17_verify_scenarios.py, 123/123 rows in Python
#   Yang-2024 CR2                        06_validate_yang2024_reference.R, 9/9 published values
#   leave-one-out, -paper, -cluster      gated against the canonical tables inside 07/14/15
#   reversal counts, rho sensitivity     two derivations recorded in the files themselves
#   FE + VCV and UWLS point estimates    closed form against metafor, 6.7e-16 / 8.9e-16
#
# WHAT THIS CLOSES. Four gaps, in descending order of how much damage each could do.
#
# A. THE `beta0_c3` GATE. This is the highest-value check available, because a defect
#    here already bit once: the missing `~` in `S2_v2.R` silently changed every
#    bias-corrected number in the paper. The gate is re-derived from its inputs -
#    `beta0` and `beta0_var_model` - by an independent expression of the same rule, and
#    compared against the stored value model by model.
#
# B. THE THREE CLOSED-FORM METRICS. Every route in this workflow calls the same three
#    functions, so an error inside one of them would survive every other check: R and
#    Python would agree because they implement the same formula. The only way out is to
#    compare against something that shares no formula at all, so power, Type M and
#    Type S are recomputed by direct Monte Carlo simulation of the sampling
#    distribution.
#
# C. THE TWO MODEL-LEVEL ESTIMANDS at the primary-study level. The Python verifier
#    covers `lmer(y ~ 1 + (1|cluster))` only. The equal-weight and random-effect
#    estimands - which produce the 22.4% and 22.3% figures - go through
#    `lmer(y ~ 0 + case + (1|g))` and `lmer(y ~ 1 + (1|case) + (1|g))` and had no second
#    derivation. Refitted here with `nlme`, a different package with a different
#    optimiser and a different parameterisation of the variance components.
#
# D. THE META-ANALYSIS-LEVEL PART A ROWS. `17_verify_scenarios.py` exercises the
#    weighted-least-squares route on parts B and C but never on the four reported
#    specifications. Closed here by a hand-computed weighted mean of logs, which shares
#    no code with `lm()`.
#
# Nothing is written to any canonical table. This writes one audit file,
# `revision/results/verification_audit.csv`, one row per quantity checked.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))
suppressMessages(library(nlme))

message("== 20: second derivations for everything reaching the manuscript ==")
S  <- readRDS(file.path(REV_TMP, "original_estimates.rds"))
o  <- S$original; L <- all_datasets(S$dat)
BR <- readRDS(file.path(REV_TMP, "bias_robust.rds"))
stopifnot(identical(o$MA_model, BR$MA_model), nrow(o) == 48L)

audit <- list()
record <- function(area, quantity, n_checked, route_a, route_b, max_abs, statistic, tol,
                   note = "", criterion = "max relative difference") {
  audit[[length(audit) + 1L]] <<- tibble::tibble(
    area = area, quantity = quantity, n_checked = n_checked,
    route_a = route_a, route_b = route_b,
    max_abs_diff = max_abs, criterion = criterion, statistic = statistic, tolerance = tol,
    passed = statistic <= tol, note = note)
  message(sprintf("  [%s] %-46s n=%-5d %-22s %.3g  %s",
          if (statistic <= tol) "PASS" else "FAIL", quantity, n_checked,
          criterion, statistic, if (nzchar(note)) note else ""))
}
relmax <- function(a, b) max(abs(a - b) / pmax(abs(a), abs(b), 1e-300))

# --- A. the beta0_c3 gate -----------------------------------------------------
# The rule, stated independently of how 01 implements it: use the variance-model
# corrected mean when its absolute value is strictly smaller than the uncorrected
# pooled mean's, otherwise keep the uncorrected mean. Magnitude only; the sign of the
# corrected mean is not consulted, which is why 20 of the 48 selected estimates have
# the opposite sign to their uncorrected counterpart.
message("\nA. the bias-correction gate")
gate_recomputed <- ifelse(abs(o$beta0_var_model) < abs(o$beta0), o$beta0_var_model, o$beta0)
record("upstream", "beta0_c3 reproduced from its own inputs", nrow(o),
       "01_reproduce_original_analysis.R", "gate rule re-expressed here",
       max(abs(gate_recomputed - o$beta0_c3)), relmax(gate_recomputed, o$beta0_c3), 1e-12)

n_sel <- sum(abs(o$beta0_var_model) < abs(o$beta0))
n_rev <- sum(sign(o$beta0_c3) != sign(o$beta0) & o$beta0_c3 != o$beta0)
message(sprintf("       gate selects the corrected mean in %d of 48 models; %d of those reverse sign",
        n_sel, n_rev))
rc <- readr::read_csv(file.path(REV_OUT, "reversal_counts.csv"), show_col_types = FALSE)
stopifnot(rc$n_reversal[rc$effect_estimator == "yang2023_gated_beta0_c3"] == n_rev)
message("       agrees with reversal_counts.csv")

# also check the flag column that records which models the gate switched
stopifnot(identical(o$gate_selects_var_model, abs(o$beta0_var_model) < abs(o$beta0)))
message("       the stored gate_selects_var_model flag agrees model by model")

# --- B. the three metrics, against Monte Carlo --------------------------------
# The definitions, simulated directly rather than evaluated in closed form:
#   power  = P(|X| > crit*se)                       for X ~ N(mu, se)
#   Type S = P(sign(X) != sign(mu) | |X| > crit*se)
#   Type M = E[|X| / |mu|  |  |X| > crit*se]
# A grid of t = |mu|/se covering the whole observed range is used rather than the
# observed values themselves, so the check does not depend on this corpus.
message("\nB. the closed-form metrics, against direct simulation")
set.seed(20260815)
N_SIM <- 4e6
t_grid <- c(0.05, 0.1, 0.25, 0.5, 1, 1.5, 1.96, 2.5, 3, 4, 6)
sim <- vapply(t_grid, function(tt) {
  x   <- stats::rnorm(N_SIM, mean = tt, sd = 1)      # se = 1 without loss of generality
  sig <- abs(x) > CRIT
  ns  <- sum(sig); pw <- ns / N_SIM
  ts  <- sum(sig & x < 0) / ns
  tm  <- mean(abs(x[sig])) / tt
  c(power = pw, type_S = ts, type_M = tm,
    # Monte Carlo standard errors, which are what a simulated value can be held to
    se_power  = sqrt(pw * (1 - pw) / N_SIM),
    se_type_S = sqrt(max(ts * (1 - ts), .Machine$double.eps) / ns),
    se_type_M = stats::sd(abs(x[sig])) / sqrt(ns) / tt)
}, numeric(6))
cf <- rbind(power  = power_two_tailed_cf(t_grid, 1),
            type_S = type_S_cf(t_grid, 1),
            type_M = type_M_cf(t_grid, 1))
# Monte Carlo standard error sets the tolerance: with 4e6 draws the tightest cell is
# the Type S ratio at small t, where the numerator is a small count.
#
# TYPE S IS CHECKED ON A SUBSET OF THE GRID, and the reason is a limit of the
# simulation rather than of the formula. Type S is the probability of a significant
# result in the WRONG direction, which is pnorm(-crit - t): at t = 4 that is 1.4e-6 of
# the significant draws, so 4e6 draws yield zero wrong-sign events and the simulated
# value is exactly 0. Comparing against 0 would report a relative error of 1 for a
# formula that is right. The grid is therefore restricted to the range where the
# simulation can resolve the quantity at all - at least 200 expected wrong-sign
# significant draws - which is t <= 1.96 here. The behaviour above that range is
# covered instead by the analytic limit check below and by the closed form's agreement
# with power and Type M, which share the same two tail probabilities.
#
# THE TOLERANCE IS THE SIMULATION'S OWN ERROR, not a fixed relative difference. A flat
# tolerance is the wrong test here: the Monte Carlo error differs by three orders of
# magnitude across the grid, so any single threshold is either vacuous at one end or
# spuriously failing at the other. Each closed-form value is instead expressed as a
# z-score against the Monte Carlo standard error of the corresponding simulated value,
# and the check is that no cell deviates by more than 4 standard errors.
resolves <- N_SIM * stats::pnorm(-CRIT - t_grid) >= 200
for (mt in rownames(cf)) {
  keep <- if (mt == "type_S") resolves else rep(TRUE, length(t_grid))
  z <- abs(cf[mt, keep] - sim[mt, keep]) / sim[paste0("se_", mt), keep]
  record("metric definitions", sprintf("%s, closed form vs %.0e-draw simulation", mt, N_SIM),
         sum(keep), "00_revision_functions.R closed form", "Monte Carlo simulation",
         max(abs(cf[mt, keep] - sim[mt, keep])), max(z), 4,
         sprintf("t from %.2f to %.2f", min(t_grid[keep]), max(t_grid[keep])),
         criterion = "max |deviation| / MC SE")
}

# the documented limiting behaviour, which is why a near-zero corrected mean gives a
# very large Type M: Type M -> dnorm(crit)/pnorm(-crit) / |t| = 2.3378/|t| as t -> 0
lim <- stats::dnorm(CRIT) / stats::pnorm(-CRIT)
record("metric definitions", "Type M limit constant 2.3378/|t| as t -> 0", 1,
       "type_M_cf at t = 1e-6", "analytic limit dnorm(c)/pnorm(-c)",
       abs(type_M_cf(1e-6, 1) * 1e-6 - lim), abs(type_M_cf(1e-6, 1) * 1e-6 / lim - 1), 1e-6,
       sprintf("constant = %.4f", lim))

# --- C. the primary-study level, lme4 against nlme ----------------------------
message("\nC. primary-study-level summaries, lme4 against nlme")
sei_all <- unlist(lapply(L, function(x) x$sei), use.names = FALSE)
sid_all <- unlist(lapply(L, function(x) as.character(x$study_ID)), use.names = FALSE)
ma_all  <- rep(o$MA_model, vapply(L, nrow, integer(1)))
stopifnot(length(sei_all) == 5740L)
cl_all  <- namespaced_study_id(ma_all, sid_all)
idx     <- rep(seq_along(L), vapply(L, nrow, integer(1)))

specs <- tibble::tribble(
  ~effect_estimator,          ~mu,                ~role,
  "uncorrected_beta0",        o$beta0,            "reference_uncorrected",
  "yang2023_gated_beta0_c3",  o$beta0_c3,         "primary",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, "reported_sensitivity",
  "yang2024_UWLS",            BR$UWLS_estimate,   "supplementary")
METRICS <- c("power", "type_M", "type_S")
metric_fun <- list(power = power_two_tailed_cf, type_M = type_M_cf, type_S = type_S_cf)

ctrl <- nlme::lmeControl(opt = "optim", maxIter = 500, msMaxIter = 500,
                         msMaxEval = 2000, returnObject = TRUE)
cmp <- list()
for (i in seq_len(nrow(specs))) for (mt in METRICS) {
  off <- offset_for(mt)
  v   <- metric_fun[[mt]](specs$mu[[i]][idx], sei_all)
  d   <- data.frame(y = log(v + off), case = factor(ma_all), g = factor(cl_all))

  a_unw <- aggregate_primary(v, cl_all, mt)$geometric_mean
  a_eq  <- aggregate_primary_equal(v, ma_all, cl_all, mt)$geometric_mean
  a_rnd <- aggregate_primary_random(v, ma_all, cl_all, mt)$geometric_mean

  b_unw <- exp(unname(nlme::fixef(nlme::lme(y ~ 1, random = ~ 1 | g, data = d,
                                            method = "REML", control = ctrl)))[1]) - off
  fe    <- nlme::fixef(nlme::lme(y ~ 0 + case, random = ~ 1 | g, data = d,
                                 method = "REML", control = ctrl))
  b_eq  <- exp(mean(fe)) - off
  b_rnd <- exp(unname(nlme::fixef(nlme::lme(y ~ 1, random = ~ 1 | case / g, data = d,
                                            method = "REML", control = ctrl)))[1]) - off

  cmp[[length(cmp) + 1L]] <- tibble::tibble(
    effect_estimator = specs$effect_estimator[i], metric = mt,
    estimand = c("study_cluster_random_intercept", "equal_per_meta_analysis",
                 "meta_analysis_random_effect"),
    lme4 = c(a_unw, a_eq, a_rnd), nlme = c(b_unw, b_eq, b_rnd))
}
cmp <- dplyr::bind_rows(cmp)
for (es in unique(cmp$estimand)) {
  x <- dplyr::filter(cmp, estimand == es)
  record("primary-study level", sprintf("%s estimand", es), nrow(x),
         "lme4::lmer", "nlme::lme",
         max(abs(x$lme4 - x$nlme)), relmax(x$lme4, x$nlme), 1e-5)
}

# and those lme4 values must be the ones in the canonical table
# (effect_estimator, metric, weighting) identifies a row uniquely: 4 x 3 x 3 = 36, the
# whole file. `role` is not part of the key because the two model-level estimands are
# filed as `secondary_descriptive` regardless of which assumed effect produced them.
pl <- readr::read_csv(file.path(REV_OUT, "primary_level_sensitivity.csv"), show_col_types = FALSE) |>
  dplyr::select(effect_estimator, metric, estimand = weighting, canonical = geometric_mean)
stopifnot(nrow(pl) == 36L, !any(duplicated(pl[, 1:3])))
j <- dplyr::inner_join(cmp, pl, by = c("effect_estimator", "metric", "estimand"))
stopifnot(nrow(j) == 36L)
record("primary-study level", "canonical table equals the refitted values", nrow(j),
       "primary_level_sensitivity.csv", "refitted here",
       max(abs(j$lme4 - j$canonical)), relmax(j$lme4, j$canonical), 1e-10)

# --- C2. WHICH estimand the reported summary actually is ----------------------
# Added 2026-08-15 after an external review found the reported primary-study-level rows
# labelled `unweighted_per_effect_size`, which they are not. Everything else in this
# workflow checks that a number is right; nothing checked that the DESCRIPTION of the
# number was right, and that is the gap the mislabel fell through. This closes it by
# measuring the estimand instead of asserting it.
#
# The check is a claim about direction, not a tolerance on a value: the fitted intercept
# must sit closer to an equal-per-study-cluster summary than to an equal-per-effect-size
# one. If that ever ceases to hold, the label has to change again.
message("\nC2. which estimand the reported primary-study-level summary is")
v_unc <- power_two_tailed_cf(specs$mu[[1]][idx], sei_all)
d_unc <- data.frame(y = log(v_unc), g = factor(cl_all))
f_unc <- lme4::lmer(y ~ 1 + (1 | g), data = d_unc)
vc    <- as.data.frame(lme4::VarCorr(f_unc))
lambda <- vc$vcov[1] / vc$vcov[2]
n_cl  <- as.numeric(table(cl_all))
w_cl  <- n_cl / (1 + n_cl * lambda)
fitted_int  <- exp(unname(lme4::fixef(f_unc))[1])
per_cluster <- exp(mean(tapply(d_unc$y, cl_all, mean)))
per_effect  <- exp(mean(d_unc$y))
message(sprintf("       lambda = tau^2/sigma^2 = %.3f | cluster sizes %d-%d over %d clusters",
        lambda, min(n_cl), max(n_cl), length(n_cl)))
message(sprintf("       implied cluster weights span %.2fx while effect-size counts span %.0fx",
        max(w_cl) / min(w_cl), max(n_cl) / min(n_cl)))
message(sprintf("       fitted %.5f | equal per cluster %.5f (%.1f%%) | equal per effect size %.5f (%.1f%%)",
        fitted_int, per_cluster, 100 * abs(fitted_int / per_cluster - 1),
        per_effect, 100 * abs(fitted_int / per_effect - 1)))
record("estimand labelling",
       "reported summary is nearer an equal-per-cluster than an equal-per-effect-size summary",
       1, "lmer intercept", "the two candidate estimands",
       abs(fitted_int - per_cluster),
       abs(fitted_int / per_cluster - 1) / abs(fitted_int / per_effect - 1), 1,
       "so the row is labelled study_cluster_random_intercept",
       criterion = "distance ratio, cluster vs effect size")

# the label must actually be the one in the canonical table
lbl <- unique(readr::read_csv(file.path(REV_OUT, "primary_level_sensitivity.csv"),
                              show_col_types = FALSE)$weighting)
if ("unweighted_per_effect_size" %in% lbl)
  stop("primary_level_sensitivity.csv still carries the retired label ",
       "`unweighted_per_effect_size`; the fitted model does not estimate that quantity")
message("       the retired label `unweighted_per_effect_size` appears nowhere in the table")

# --- D. the meta-analysis level, lm() against a hand-computed weighted mean ----
message("\nD. meta-analysis-level summaries, lm() against a hand-computed weighted mean")
se_for <- list(uncorrected_beta0 = o$se_beta0, yang2023_gated_beta0_c3 = o$se_beta0,
               yang2024_FE_VCV = BR$FE_VCV_CRVE_SE_CR2,
               yang2024_UWLS = BR$UWLS_CRVE_SE_CR2_naive_t)
ma_cmp <- list()
for (i in seq_len(nrow(specs))) for (mt in METRICS) {
  off <- offset_for(mt)
  v   <- metric_fun[[mt]](specs$mu[[i]], se_for[[specs$effect_estimator[i]]])
  for (wn in c("k_effect_sizes", "equal_per_meta_analysis")) {
    w <- if (wn == "k_effect_sizes") o$k else rep(1, length(v))
    ma_cmp[[length(ma_cmp) + 1L]] <- tibble::tibble(
      effect_estimator = specs$effect_estimator[i], role = specs$role[i],
      metric = mt, weighting = wn,
      via_lm = aggregate_ma(v, w, mt)$geometric_mean,
      by_hand = exp(sum(w * log(v + off)) / sum(w)) - off)
  }
}
ma_cmp <- dplyr::bind_rows(ma_cmp)
record("meta-analysis level", "weighted geometric mean, both weightings", nrow(ma_cmp),
       "stats::lm(weights = k)", "hand-computed sum(w*log(x))/sum(w)",
       max(abs(ma_cmp$via_lm - ma_cmp$by_hand)), relmax(ma_cmp$via_lm, ma_cmp$by_hand), 1e-12)

# Dropping the `diagnostic_*` rows leaves exactly the 24 reportable ones, on which
# (effect_estimator, metric, weighting) is unique. `role` is not part of the key: the
# equal-weight rows are filed as `secondary_descriptive` whatever produced them, and
# the diagnostics deliberately share an estimator with a reportable row.
ml <- readr::read_csv(file.path(REV_OUT, "meta_analysis_level_sensitivity.csv"),
                      show_col_types = FALSE) |>
  dplyr::filter(!startsWith(role, "diagnostic")) |>
  dplyr::select(effect_estimator, metric, weighting, canonical = geometric_mean)
stopifnot(nrow(ml) == 24L, !any(duplicated(ml[, 1:3])))
jm <- dplyr::inner_join(ma_cmp, ml, by = c("effect_estimator", "metric", "weighting"))
stopifnot(nrow(jm) == 24L)
record("meta-analysis level", "canonical table equals the recomputed values", nrow(jm),
       "meta_analysis_level_sensitivity.csv", "recomputed here",
       max(abs(jm$via_lm - jm$canonical)), relmax(jm$via_lm, jm$canonical), 1e-12)

# --- summary ------------------------------------------------------------------
audit <- dplyr::bind_rows(audit) |>
  dplyr::mutate(checked_on = "2026-08-15", .before = 1)
write_revision(audit, "verification_audit.csv")
message(sprintf("\n%d checks, %d passed, %d failed",
        nrow(audit), sum(audit$passed), sum(!audit$passed)))
if (any(!audit$passed)) stop("verification failed: see revision/results/verification_audit.csv")
