# 05_loo_study.R -------------------------------------------------------------
# B5: leave-one-primary-study-out assumed effect (Reviewer 2, §3.6).
#
# Reviewer 2's objection is that each primary study contributes to the pooled
# effect used to evaluate it, so the assumed effect is not independent of the
# study under assessment. This recomputes the pooled mean with the focal study
# excluded and re-evaluates that study's metrics against it, holding the study's
# own observed sei fixed.
#
# EXCLUSION UNIT: the (study_ID x dataset) cluster -- every row sharing a
# study_ID within one dataset, not a single row. 903 of the 1,415 clusters
# (63.8%) contain more than one effect size, up to 115, so dropping one row
# would leave most of a focal study's contribution in the proxy.
#
# Note the limit established in the audit: study_ID is not a harmonised primary-
# study identifier. Four schemes coexist across the 48 datasets (author-year with
# format drift, numeric PMID-style in SMD/MA31, opaque codes CD001-CD126 in
# lnRR/MA09, and author-only with no year in MA14/MA17), and in 9 source papers
# study_ID is FINER than "primary study" -- MA41 has 106 labels against 54
# referenced studies. So this analysis removes a study *as the source dataset
# defines it*. Where study_ID is finer than a paper, self-inclusion is only
# partially removed, which makes the result a LOWER bound on the self-inclusion
# effect. It cannot be tightened without the identifier work in §2 of the plan.
#
# UNCORRECTED ONLY. The corrected analogue would need the detection model, the
# scenario assignment, the scenario-specific reduced model, and the gate re-run
# per exclusion (~4 fits x 1,415), and a dropped study can flip a scenario or
# the gate. The resulting spread would then mix self-inclusion with scenario
# instability and would not be interpretable as "the same analysis minus one
# study". Deferred pending the uncorrected result (plan §4b, decision C1).

source(here::here("R", "02_metrics.R"))

e <- readRDS(out_dir("estimates.rds"))
dat <- e$dat; est <- e$est
L <- all_datasets(dat)
stopifnot(identical(names(L), est$case))

clusters <- purrr::list_rbind(lapply(seq_along(L), function(i) {
  tibble::tibble(dataset_index = i, case = est$case[i], es_type = est$es_type[i],
                 study_ID = unique(as.character(L[[i]]$study_ID)))
}))
message(sprintf("leave-one-study-out: %d (study x dataset) clusters across %d datasets",
                nrow(clusters), length(L)))

# For each cluster: refit the intercept-only model without it, then recompute
# the focal cluster's own rows against the leave-one-out pooled mean.
loo_rows <- vector("list", nrow(clusters))
n_failed <- 0L

for (j in seq_len(nrow(clusters))) {
  i <- clusters$dataset_index[j]
  sid <- clusters$study_ID[j]
  d <- L[[i]]
  keep <- as.character(d$study_ID) != sid
  focal <- d[!keep, , drop = FALSE]

  # a dataset whose only study is the focal one cannot yield a LOO proxy
  if (sum(keep) < 2L || dplyr::n_distinct(d$study_ID[keep]) < 2L) {
    n_failed <- n_failed + 1L
    next
  }

  fit <- fit_rma(d[keep, , drop = FALSE],
                 optimizer = optimizer_for(est$es_type[i], "null"))
  if (is.null(fit)) { n_failed <- n_failed + 1L; next }

  mu_loo <- as.numeric(fit$beta)
  mu_self <- est$beta0[i]

  loo_rows[[j]] <- tibble::tibble(
    case = est$case[i], es_type = est$es_type[i], paper = est$paper[i],
    study_ID = sid, sei = focal$sei,
    mu_self_inclusive = mu_self, mu_leave_one_out = mu_loo,
    power_self = power_two_tailed(mu_self, focal$sei),
    power_loo  = power_two_tailed(mu_loo,  focal$sei),
    type_M_self = error_M_vec(mu_self, focal$sei),
    type_M_loo  = error_M_vec(mu_loo,  focal$sei),
    type_S_self = error_S(mu_self, focal$sei),
    type_S_loo  = error_S(mu_loo,  focal$sei)
  )

  if (j %% 200L == 0L) message(sprintf("  %d / %d clusters", j, nrow(clusters)))
}

loo <- purrr::list_rbind(loo_rows)
message(sprintf("clusters with a usable LOO proxy: %d of %d (%d skipped)",
                nrow(clusters) - n_failed, nrow(clusters), n_failed))
readr::write_csv(loo, out_dir("B5_leave_one_study_out_rows.csv"))

# aggregate both ways on the same rows, so the comparison is like-for-like
agg <- function(col_self, col_loo, metric) {
  off <- offset_for(metric)
  a <- summarise_primary(
    dplyr::mutate(loo, !!metric := .data[[col_self]], scenario_label = "self-inclusive"),
    metric)
  b <- summarise_primary(
    dplyr::mutate(loo, !!metric := .data[[col_loo]], scenario_label = "leave-one-out"),
    metric)
  dplyr::bind_rows(a, b)
}

b5 <- purrr::list_rbind(list(
  agg("power_self",  "power_loo",  "power"),
  agg("type_M_self", "type_M_loo", "type_M"),
  agg("type_S_self", "type_S_loo", "type_S")
)) |> dplyr::mutate(analysis = "B5 leave-one-study-out", .before = 1)
readr::write_csv(b5, out_dir("B5_leave_one_study_out_summary.csv"))

b5_delta <- b5 |>
  dplyr::select(metric, scenario_label, model_median) |>
  tidyr::pivot_wider(names_from = scenario_label, values_from = model_median) |>
  dplyr::mutate(abs_change = `leave-one-out` - `self-inclusive`,
                pct_change = 100 * (`leave-one-out` / `self-inclusive` - 1))
readr::write_csv(b5_delta, out_dir("B5_delta.csv"))

b5_shift <- tibble::tibble(
  n_rows = nrow(loo),
  n_clusters = dplyr::n_distinct(paste(loo$case, loo$study_ID)),
  median_abs_mu_shift = stats::median(abs(loo$mu_leave_one_out - loo$mu_self_inclusive)),
  max_abs_mu_shift = max(abs(loo$mu_leave_one_out - loo$mu_self_inclusive)),
  median_pct_mu_shift = stats::median(100 * abs(loo$mu_leave_one_out / loo$mu_self_inclusive - 1)),
  n_sign_flips = sum(sign(loo$mu_leave_one_out) != sign(loo$mu_self_inclusive)),
  median_power_change_pp = stats::median(100 * (loo$power_loo - loo$power_self))
)
readr::write_csv(b5_shift, out_dir("B5_proxy_shift.csv"))

cat("\n--- B5: leave-one-study-out vs self-inclusive ---\n")
print(as.data.frame(b5_delta), digits = 4)
cat("\n--- B5: how much the assumed effect moves ---\n")
print(as.data.frame(b5_shift), digits = 4)
