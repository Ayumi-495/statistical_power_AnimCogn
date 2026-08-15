# R/revision/15_leave_one_paper_out.R -------------------------------------------
# Step 15: leave-one-SOURCE-PAPER-out influence on the meta-analysis-level summaries.
#
# A COMPANION TO 07, NOT A REPLACEMENT. `07_influence_loo.R` removes one of the 48
# meta-analysis MODELS at a time. This removes one of the 28 source PAPERS at a time,
# and therefore drops every model that paper contributed. Twelve papers contribute more
# than one model, so the two analyses ask different questions:
#
#   07  how much does the summary depend on any single fitted model?
#   15  how much does the summary depend on any single published paper?
#
# The second is the one a reader asks about a second-order meta-analysis, because a
# paper is the unit of publication and of any shared authorship, laboratory, or
# analytical convention. Both are kept.
#
# The earlier audit-phase output `outputs/B2_leave_one_paper_out.csv` is NOT used and
# must not be quoted alongside these numbers. It predates the adopted clustering, used a
# Monte Carlo Type M, and covered only the uncorrected and Yang-2023 specifications.
#
# Everything here is aligned with the current pipeline:
#   - the current corrected estimates and model definitions (from 01 and 03)
#   - four specifications, each with the standard error it is reported with
#   - closed-form Type M
#   - Type S summarised through `aggregate_ma()`, so the 0.025 offset, the raw quartiles
#     and the `summary_dominated_by_offset` flag travel with every row exactly as they do
#     in the canonical tables
#   - both meta-analysis-level weightings
#   - the aggregate recomputed from ONLY the models that remain after the paper is removed

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 15: leave-one-source-paper-out, meta-analysis level ==")
o  <- readRDS(file.path(REV_TMP, "original_estimates.rds"))$original
BR <- readRDS(file.path(REV_TMP, "bias_robust.rds"))
stopifnot(identical(o$MA_model, BR$MA_model), nrow(o) == 48L)

papers <- sort(unique(o$source_paper))
message(sprintf("%d source papers contributing %d models; %d papers contribute more than one",
        length(papers), nrow(o), sum(table(o$source_paper) > 1)))

METRICS <- c("power", "type_M", "type_S")
metric_fun <- list(power = power_two_tailed_cf, type_M = type_M_cf, type_S = type_S_cf)

# The four reported specifications, each paired with the standard error it is reported
# with, matching 04_revision_sensitivity_summaries.R. Critical value z = 1.96 throughout.
specs <- tibble::tribble(
  ~effect_estimator,          ~mu,                ~se,                          ~se_source,  ~se_method,          ~spec_role,
  "uncorrected_beta0",        o$beta0,            o$se_beta0,                   "se_beta0",  NA_character_,       "reference_uncorrected",
  "yang2023_gated_beta0_c3",  o$beta0_c3,         o$se_beta0,                   "se_beta0",  NA_character_,       "primary",
  "yang2024_FE_VCV",          BR$FE_VCV_estimate, BR$FE_VCV_CRVE_SE_CR2,        "own_CRVE",  "CR2_Satterthwaite", "reported_sensitivity",
  "yang2024_UWLS",            BR$UWLS_estimate,   BR$UWLS_CRVE_SE_CR2_naive_t,  "own_CRVE",  "CR2_naive_t",       "supplementary"
)

summarise <- function(v, k, mt) {
  dplyr::bind_rows(
    aggregate_ma(v, k, mt)                 |> dplyr::mutate(weighting = "k_effect_sizes"),
    aggregate_ma(v, rep(1, length(v)), mt) |> dplyr::mutate(weighting = "equal_per_meta_analysis"))
}

lopo <- purrr::list_rbind(lapply(seq_len(nrow(specs)), function(i) {
  sp <- specs[i, ]; mu <- sp$mu[[1]]; se <- sp$se[[1]]
  purrr::list_rbind(lapply(METRICS, function(mt) {
    v <- metric_fun[[mt]](mu, se)
    base <- summarise(v, o$k, mt) |>
      dplyr::select(weighting, base = geometric_mean)
    purrr::list_rbind(lapply(papers, function(pp) {
      keep <- o$source_paper != pp        # every model this paper contributed is removed
      summarise(v[keep], o$k[keep], mt) |>
        dplyr::mutate(
          dropped_source_paper = pp,
          n_models_dropped = sum(!keep), n_models_remaining = sum(keep),
          n_effect_sizes_dropped = sum(o$k[!keep]),
          n_effect_sizes_remaining = sum(o$k[keep]),
          pct_of_effect_sizes_dropped = 100 * sum(o$k[!keep]) / sum(o$k),
          metric = mt, effect_estimator = sp$effect_estimator,
          se_source = sp$se_source, se_method = sp$se_method,
          spec_role = sp$spec_role, .before = 1)
    })) |>
      dplyr::left_join(base, by = "weighting") |>
      dplyr::mutate(summary_all_28_papers = base,
                    abs_change = geometric_mean - base,
                    pct_change = 100 * (geometric_mean / base - 1)) |>
      dplyr::select(-base) |>
      dplyr::group_by(weighting) |>
      dplyr::mutate(influence_rank = rank(-abs(pct_change), ties.method = "min")) |>
      dplyr::ungroup()
  }))
})) |>
  dplyr::mutate(aggregation = "meta_analysis_level", crit_value_method = "z_1.96",
                exclusion_unit = "source_paper", role = "influence_check",
                verification_status = "two_derivations", .before = 1)

# --- verification gates -------------------------------------------------------
expected <- length(papers) * nrow(specs) * length(METRICS) * 2L
stopifnot(nrow(lopo) == expected)
message(sprintf("rows: %d = %d papers x %d specifications x %d metrics x 2 weightings",
        nrow(lopo), length(papers), nrow(specs), length(METRICS)))

combo <- table(lopo$dropped_source_paper, lopo$effect_estimator,
               lopo$metric, lopo$weighting)
if (!all(combo == 1L))
  stop("incomplete or duplicated combinations in the leave-one-paper-out grid")
message("complete combination check: every (paper x specification x metric x weighting) appears exactly once")

# every paper must remove at least one model and leave at least two behind
stopifnot(all(lopo$n_models_dropped >= 1L),
          all(lopo$n_models_remaining + lopo$n_models_dropped == 48L),
          all(lopo$n_effect_sizes_remaining + lopo$n_effect_sizes_dropped == 5740L))
message(sprintf("models remaining after a paper is dropped: min %d, max %d (of 48)",
        min(lopo$n_models_remaining), max(lopo$n_models_remaining)))

# the baseline must equal the canonical summary
canon <- readr::read_csv(file.path(REV_OUT, "meta_analysis_level_sensitivity.csv"),
                         show_col_types = FALSE) |>
  dplyr::filter(role %in% c("reference_uncorrected", "primary", "reported_sensitivity",
                            "supplementary", "secondary_descriptive"),
                weighting %in% c("k_effect_sizes", "equal_per_meta_analysis")) |>
  dplyr::select(effect_estimator, metric, weighting, canonical = geometric_mean)
cmp <- lopo |>
  dplyr::distinct(effect_estimator, metric, weighting, summary_all_28_papers) |>
  dplyr::inner_join(canon, by = c("effect_estimator", "metric", "weighting"))
d <- max(abs(cmp$summary_all_28_papers - cmp$canonical))
message(sprintf("baseline agreement with meta_analysis_level_sensitivity.csv over %d cells: max|diff| = %.3g",
        nrow(cmp), d))
if (nrow(cmp) < 24L || d > 1e-10)
  stop("leave-one-paper-out baselines do not match the canonical summaries")

write_revision(lopo, "leave_one_paper_out.csv")

# --- what it shows ------------------------------------------------------------
per_paper <- o |>
  dplyr::count(source_paper, name = "n_models") |>
  dplyr::left_join(dplyr::summarise(dplyr::group_by(o, source_paper),
                                    k = sum(k), .groups = "drop"), by = "source_paper")
message("\nmodels contributed per paper:")
print(as.data.frame(dplyr::arrange(per_paper, dplyr::desc(n_models), dplyr::desc(k))[1:6, ]),
      row.names = FALSE)

for (r in c("primary", "reported_sensitivity")) {
  x <- lopo |>
    dplyr::filter(metric == "power", spec_role == r, weighting == "k_effect_sizes") |>
    dplyr::arrange(dplyr::desc(abs(pct_change)))
  message(sprintf("\n%s, power, k-weighted: summary over all 28 papers %.5f", r, x$summary_all_28_papers[1]))
  for (j in 1:3)
    message(sprintf("   drop %-6s (%d model(s), %4.1f%% of effect sizes) -> %.5f (%+.1f%%)",
            x$dropped_source_paper[j], x$n_models_dropped[j],
            x$pct_of_effect_sizes_dropped[j], x$geometric_mean[j], x$pct_change[j]))
  message(sprintf("   median |change| %.1f%% | above 10%%: %d of 28 | above 20%%: %d",
          stats::median(abs(x$pct_change)), sum(abs(x$pct_change) > 10), sum(abs(x$pct_change) > 20)))
  y <- lopo |>
    dplyr::filter(metric == "power", spec_role == r, weighting == "equal_per_meta_analysis") |>
    dplyr::arrange(dplyr::desc(abs(pct_change)))
  message(sprintf("   equally weighted: baseline %.5f | largest %s %+.1f%% | median |change| %.1f%%",
          y$summary_all_28_papers[1], y$dropped_source_paper[1], y$pct_change[1],
          stats::median(abs(y$pct_change))))
}
