# R/revision/14_leave_one_cluster_out.R -----------------------------------------
# Step 14: leave-one-cluster-out at the primary-study level, under the adopted
# clustering definition.
#
# Reviewer 2 (R2C9) objects that each primary study contributes to the pooled mean used
# to evaluate it, so the assumed effect is not independent of the study being assessed.
# The remedy is to recompute the pooled mean with the focal study's cluster removed and
# re-evaluate that cluster's own effect sizes against it, holding each effect size's own
# sampling error fixed.
#
# WHAT NEEDED REDOING, AND WHAT DID NOT. The earlier run already removed the correct
# unit: every row sharing a `study_ID` WITHIN one meta-analysis, which is exactly the
# clustering definition later adopted. So the leave-one-out means themselves stand and
# no model needs refitting. What was wrong was the second stage: the summary across the
# 5,740 values was taken with a random effect on the raw `study_ID`, which merges rows
# from different meta-analyses into one cluster. That aggregation is redone here on
# `MA_model x study_ID`, so the reported 17.15% -> 17.22% comparison is superseded.
#
# Metrics are recomputed in closed form from the stored `sei` and the two assumed
# effects rather than carried over, because the earlier run estimated Type M by
# simulation and everything in results/revision/ is analytic.
#
# UNCORRECTED ONLY, as before. A corrected analogue would have to re-run the detection
# model, the scenario assignment and the magnitude comparison for each of the ~1,415
# exclusions, and a dropped cluster can change the scenario or the selection. The
# resulting spread would mix self-inclusion with selection instability and would not be
# interpretable as "the same analysis minus one study".

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 14: leave-one-cluster-out, adopted clustering ==")

SRC <- here::here("outputs", "B5_leave_one_study_out_rows.csv")
if (!file.exists(SRC))
  stop("Row-level leave-one-out output not found: ", SRC,
       "\nRegenerate it with R/05_loo_study.R before running this script.")

b5 <- readr::read_csv(SRC, show_col_types = FALSE)
stopifnot(nrow(b5) == 5740L)

# The adopted clustering unit. `case` is the meta-analysis model.
b5$cluster <- namespaced_study_id(b5$case, as.character(b5$study_ID))
message(sprintf("%d effect sizes, %d clusters under MA_model x study_ID (%d under the raw identifier)",
        nrow(b5), dplyr::n_distinct(b5$cluster), dplyr::n_distinct(b5$study_ID)))

# --- gate: the stored self-inclusive assumed effect must be the pooled mean ----
o <- readRDS(file.path(REV_TMP, "original_estimates.rds"))$original
chk <- dplyr::left_join(dplyr::distinct(b5, case, mu_self_inclusive),
                        dplyr::transmute(o, case = MA_model, beta0), by = "case")
d0 <- max(abs(chk$mu_self_inclusive - chk$beta0))
message(sprintf("stored self-inclusive assumed effect vs beta0: max|diff| = %.3g", d0))
if (d0 > 1e-8)
  stop("The stored leave-one-out file does not correspond to the current estimates.")

n_missing <- sum(is.na(b5$mu_leave_one_out))
message(sprintf("clusters with no leave-one-out mean (the meta-analysis has too few clusters): %d effect sizes",
        n_missing))

# --- re-aggregate -------------------------------------------------------------
METRICS <- c("power", "type_M", "type_S")
metric_fun <- list(power = power_two_tailed_cf, type_M = type_M_cf, type_S = type_S_cf)

use <- b5[!is.na(b5$mu_leave_one_out), ]
rows <- purrr::list_rbind(lapply(c("self_inclusive", "leave_one_cluster_out"), function(sc) {
  mu <- if (sc == "self_inclusive") use$mu_self_inclusive else use$mu_leave_one_out
  purrr::list_rbind(lapply(METRICS, function(mt) {
    v <- metric_fun[[mt]](mu, use$sei)
    dplyr::bind_rows(
      aggregate_primary(v, use$cluster, mt) |>
        dplyr::mutate(weighting = "unweighted_per_effect_size"),
      aggregate_primary_equal(v, use$case, use$cluster, mt) |>
        dplyr::mutate(weighting = "equal_per_meta_analysis")
    ) |> dplyr::mutate(metric = mt, assumed_effect = sc, .before = 1)
  }))
})) |>
  dplyr::mutate(aggregation = "primary_study_level",
                se_source = "own_sampling_error_per_effect_size",
                crit_value_method = "z_1.96",
                clustering_unit = "MA_model_x_study_ID",
                effect_estimator = "uncorrected_beta0",
                role = "influence_check",
                verification_status = "two_derivations", .before = 1)

write_revision(rows, "leave_one_cluster_out.csv")

# --- gate: the self-inclusive row must match the canonical summary -------------
canon <- readr::read_csv(file.path(REV_OUT, "primary_level_sensitivity.csv"),
                         show_col_types = FALSE) |>
  dplyr::filter(effect_estimator == "uncorrected_beta0",
                weighting == "unweighted_per_effect_size")
cmp <- dplyr::inner_join(
  dplyr::filter(rows, assumed_effect == "self_inclusive",
                weighting == "unweighted_per_effect_size") |>
    dplyr::select(metric, loo_file = geometric_mean),
  dplyr::select(canon, metric, canonical = geometric_mean), by = "metric")
message("\nself-inclusive rows against primary_level_sensitivity.csv:")
for (i in seq_len(nrow(cmp)))
  message(sprintf("   %-7s recomputed %.6f | canonical %.6f | diff %.2e",
          cmp$metric[i], cmp$loo_file[i], cmp$canonical[i],
          abs(cmp$loo_file[i] - cmp$canonical[i])))

# --- the comparison the response letter needs ---------------------------------
pw <- function(sc, w) rows$geometric_mean[rows$metric == "power" &
                                          rows$assumed_effect == sc & rows$weighting == w]
message(sprintf("\nPOWER, effect-size-weighted : self-inclusive %.5f -> leave-one-cluster-out %.5f (%+.2f%%)",
        pw("self_inclusive", "unweighted_per_effect_size"),
        pw("leave_one_cluster_out", "unweighted_per_effect_size"),
        100 * (pw("leave_one_cluster_out", "unweighted_per_effect_size") /
               pw("self_inclusive", "unweighted_per_effect_size") - 1)))
message(sprintf("POWER, meta-analysis-weighted: self-inclusive %.5f -> leave-one-cluster-out %.5f",
        pw("self_inclusive", "equal_per_meta_analysis"),
        pw("leave_one_cluster_out", "equal_per_meta_analysis")))
tm <- function(sc) rows$geometric_mean[rows$metric == "type_M" &
                                       rows$assumed_effect == sc &
                                       rows$weighting == "unweighted_per_effect_size"]
message(sprintf("TYPE M, effect-size-weighted : %.4f -> %.4f", tm("self_inclusive"), tm("leave_one_cluster_out")))
