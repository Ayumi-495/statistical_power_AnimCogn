# revision/R/analyse/12_supplementary_tables.R ------------------------------------------
# Step 12: the two supplementary tables that go into the manuscript.
#
# Only two. Everything else stays as CSV in the repository and in the archived,
# DOI-bearing deposit, which is also the answer to the FAIR / DOI / licence request:
# large tables belong in a citable archive rather than compressed into a PDF.
#
#   Table S1  the characteristics of the 28 included meta-analytical papers. Four
#             reviewer comments depend on it (R1C17, R2C11, R2C15, R2C17).
#   Table S2  the reported metrics, both levels, every assumed effect, with the
#             summary, its interval, and the raw quantiles. This is what a reader
#             needs in order to check the claims in the text.
#   Table S3  leave-one-MODEL-out influence, all 48 models.
#   Table S4  leave-one-PAPER-out influence, all 28 source papers.
#
# S3 AND S4 EXIST BECAUSE THE MAIN TEXT QUOTES NUMBERS OUT OF THEM. The text gives the
# largest single-model influence (MA09, +44.7% on the bias-robust summary) and the
# largest single-paper influence under equal weighting. A reader who wants to check
# either, or to see where it sits among the other 47 or 27 values, should not have to
# clone the repository to do it. The underlying CSVs stay in the archived deposit; these
# are the display versions, rounded and labelled like S1 and S2.
#
# They are sorted by influence within each block, so the value the text quotes is the
# FIRST ROW of its specification x metric x weighting group and the reader can see at a
# glance how far it stands above the rest.
#
# All four are written as CSV, for pasting into the manuscript.
#
# TYPE S IS REPORTED BOTH WAYS AND THE TABLE SAYS SO. The summary is fitted on
# log(Type S + 0.025) and the offset dominates when the values are far below it, so
# every row carries the raw median and quartiles beside the model-based value. Which
# of the two becomes the headline is still with the PI; the table does not prejudge it.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))

message("== 12: supplementary tables ==")
SUP <- file.path(REV_OUT, "supplementary"); dir.create(SUP, showWarnings = FALSE)

# --- constraining out-of-range confidence limits ------------------------------
# The intervals are back-transformed from models fitted on the log scale and are not
# bounded, so a limit can fall outside the range the metric can actually take. Following
# the convention used in the submitted supplement, such limits are constrained to the
# bound and marked with an asterisk. No POINT estimate falls outside its range; only
# interval limits do (13 upper limits for power, 24 lower limits for Type S, and 13
# lower limits for Type M).
#
#   power   in [0, 1]
#   Type S  in [0, 0.5]; 0.5 is the value taken when the assumed effect is zero
#   Type M  at least 1; a significant estimate cannot be smaller than the assumed effect
#           in expectation
BOUNDS <- list("Statistical power" = c(0, 1), "Type S error" = c(0, 0.5),
               "Type M error" = c(1, Inf))
constrain <- function(v, metric) {
  b <- BOUNDS[[metric]]
  pmin(pmax(v, b[1]), b[2])
}
was_constrained <- function(v, metric) {
  b <- BOUNDS[[metric]]
  !is.na(v) & (v < b[1] | v > b[2])
}

fmt <- function(v, metric, digits = 3) {
  ifelse(metric == "type_M", sprintf("%.2f", v),
         ifelse(abs(v) < 0.001 & v != 0, sprintf("%.2g", v), sprintf("%.4f", v)))
}
pct <- function(v) sprintf("%.2f", 100 * v)

# --- Table S2: reported metrics -----------------------------------------------
core <- dplyr::bind_rows(
  readr::read_csv(file.path(REV_OUT, "primary_level_sensitivity.csv"), show_col_types = FALSE),
  readr::read_csv(file.path(REV_OUT, "meta_analysis_level_sensitivity.csv"), show_col_types = FALSE)) |>
  dplyr::filter(role %in% c("reference_uncorrected", "primary", "reported_sensitivity",
                            "supplementary", "secondary_descriptive")) |>
  dplyr::transmute(
    part = "A. Reported results",
    level = ifelse(aggregation == "primary_study_level",
                   "Primary-study level", "Meta-analysis level"),
    assumed_effect = dplyr::recode(effect_estimator,
      uncorrected_beta0 = "Uncorrected pooled mean",
      yang2023_gated_beta0_c3 = "Yang 2023 bias-corrected",
      yang2024_FE_VCV = "Yang 2024 bias-robust (FE + VCV)",
      yang2024_UWLS = "Yang 2024 bias-robust (UWLS)"),
    # Display names. `.default = NA` rather than passing the input through, so that an
    # unmapped internal label fails the gate below instead of leaking into a
    # supplementary table.
    weighting = dplyr::recode(weighting,
      unweighted = "-",
      # NOT "each effect-size estimate weighted equally". That was the label until
      # 2026-08-15 and it was wrong: the random intercept makes this close to an
      # equal-per-study-cluster summary. See the note above `aggregate_primary()`.
      study_cluster_random_intercept = "Study cluster as a random effect",
      k_effect_sizes              = "By effect-size count",
      equal_per_meta_analysis     = "Each meta-analysis weighted equally",
      meta_analysis_random_effect = "Meta-analysis as a random effect",
      .default = NA_character_),
    metric, n_unit, geometric_mean, ci_lower, ci_upper,
    raw_median, raw_q1, raw_q3, summary_dominated_by_offset)

sc_all <- readr::read_csv(file.path(REV_OUT, "assumed_effect_scenarios.csv"),
                          show_col_types = FALSE)

scen <- sc_all |>
  dplyr::filter(aggregation == "primary_study_level",
                scenario_family %in% c("optimistic", "external"),
                grouping != "all metrics" | scenario_family == "optimistic") |>
  dplyr::transmute(
    part = "B. Sensitivity to the assumed effect (primary-study level)",
    level = ifelse(grouping == "all metrics", "All metrics", grouping),
    assumed_effect = scenario, weighting = "-",
    metric, n_unit, geometric_mean, ci_lower, ci_upper,
    raw_median, raw_q1, raw_q3, summary_dominated_by_offset)

# Part C. The externally specified assumed effects at the META-ANALYSIS level. This is
# the comparison Reviewer 2 asks for: what these same models could detect when the
# assumed effect does not come from the models themselves. It was previously computed
# but filtered out of the table.
scen_ma <- sc_all |>
  dplyr::filter(aggregation == "meta_analysis_level", scenario_family == "external") |>
  dplyr::transmute(
    part = "C. Externally specified assumed effects (meta-analysis level)",
    level = grouping, assumed_effect = scenario,
    weighting = dplyr::recode(weighting,
      k_effect_sizes          = "By effect-size count",
      equal_per_meta_analysis = "Each meta-analysis weighted equally",
      .default = NA_character_),
    metric, n_unit, geometric_mean, ci_lower, ci_upper,
    raw_median, raw_q1, raw_q3, summary_dominated_by_offset)

# --- part D: leave-one-cluster-out --------------------------------------------
# The Methods say the leave-one-cluster-out analysis recomputed power, Type M error AND
# Type S error, and it did - `14_leave_one_cluster_out.R` writes all three. But until
# 2026-08-17 none of it reached this table, so a reader following the Methods to the
# supplement found the optimistic and external scenarios and nothing for the dependence
# check. All twelve rows are included here: two conditions (the assumed effect estimated
# with the cluster in, and with it held out) x three metrics x the two primary-study-level
# estimands. Only the uncorrected pooled mean is used, because the question is whether an
# effect size helps estimate the mean it is judged against, and that dependence is the
# same whichever correction is applied afterwards.
loco <- readr::read_csv(file.path(REV_OUT, "leave_one_cluster_out.csv"),
                        show_col_types = FALSE) |>
  dplyr::transmute(
    part = "D. Leave-one-cluster-out (primary-study level)",
    level = "Primary-study level",
    assumed_effect = dplyr::recode(assumed_effect,
      self_inclusive        = "Uncorrected pooled mean, cluster included",
      leave_one_cluster_out = "Uncorrected pooled mean, cluster held out",
      .default = NA_character_),
    weighting = dplyr::recode(weighting,
      study_cluster_random_intercept = "Study cluster as a random effect",
      equal_per_meta_analysis        = "Each meta-analysis weighted equally",
      .default = NA_character_),
    metric, n_unit, geometric_mean, ci_lower, ci_upper,
    raw_median, raw_q1, raw_q3, summary_dominated_by_offset)
stopifnot(nrow(loco) == 12L, !any(is.na(loco$assumed_effect)), !any(is.na(loco$weighting)))

S2 <- dplyr::bind_rows(core, scen, scen_ma, loco) |>
  dplyr::mutate(
    metric = dplyr::recode(metric, power = "Statistical power",
                           type_M = "Type M error", type_S = "Type S error"),
    summary_estimate = fmt(geometric_mean, metric),
    ci = sprintf("[%s%s, %s%s]",
                 fmt(mapply(constrain, ci_lower, metric), metric),
                 ifelse(mapply(was_constrained, ci_lower, metric), "*", ""),
                 fmt(mapply(constrain, ci_upper, metric), metric),
                 ifelse(mapply(was_constrained, ci_upper, metric), "*", "")),
    n_constrained_limits = mapply(was_constrained, ci_lower, metric) +
                           mapply(was_constrained, ci_upper, metric),
    raw_median_iqr = sprintf("%s [%s, %s]", fmt(raw_median, metric),
                             fmt(raw_q1, metric), fmt(raw_q3, metric)),
    offset_note = ifelse(summary_dominated_by_offset,
                         "summary sensitive to the 0.025 offset", "")) |>
  dplyr::select(part, level, assumed_effect, weighting, metric, n_unit,
                summary_estimate, ci, n_constrained_limits, raw_median_iqr, offset_note)

# No point estimate should ever need constraining; if one does, that is a modelling
# problem rather than a display problem and must not be hidden by an asterisk.
bad <- dplyr::bind_rows(core, scen, scen_ma, loco) |>
  dplyr::mutate(metric = dplyr::recode(metric, power = "Statistical power",
                                       type_M = "Type M error", type_S = "Type S error")) |>
  dplyr::filter(mapply(was_constrained, geometric_mean, metric))
if (nrow(bad)) stop(sprintf("%d point estimates fall outside the metric's range", nrow(bad)))

# No internal label may reach the table.
if (any(is.na(S2$weighting)))
  stop("unmapped `weighting` value reached Table S2; add it to the recode above")
leak <- grep("_", c(S2$weighting, S2$level, S2$assumed_effect, S2$metric), value = TRUE)
if (length(leak)) stop("internal shorthand in Table S2: ", paste(unique(leak), collapse = ", "))

unlink(file.path(SUP, "TableS1_reported_metrics.csv"))
readr::write_csv(S2, file.path(SUP, "TableS2_reported_metrics.csv"))
message(sprintf("Table S2: %d rows -> supplementary/TableS2_reported_metrics.csv", nrow(S2)))

# --- Table S1: evidence-base characteristics ----------------------------------
short_species <- function(x) {
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    v <- trimws(unlist(strsplit(s, ";")))
    v <- v[nzchar(v)]
    if (length(v) <= 3) paste(v, collapse = "; ")
    else sprintf("%s and %d other species", paste(v[1:2], collapse = "; "), length(v) - 2)
  }, character(1), USE.NAMES = FALSE)
}
clip <- function(x, n = 80) ifelse(is.na(x) | nchar(x) <= n, x,
                                   paste0(substr(x, 1, n - 1), "…"))

S1 <- readr::read_csv(file.path(REV_OUT, "evidence_base_characteristics.csv"),
                      show_col_types = FALSE) |>
  dplyr::transmute(
    paper = meta_id, study = authors_id, year,
    species = short_species(species_focus),
    cognitive_domain = cognition_domain,
    task = clip(task),
    manipulation = clip(manipulation_category, 40),
    life_stage = ifelse(tolower(trimws(life_stage_reported)) %in% c("yes","y","true"),
                        ifelse(tolower(trimws(life_stage_included)) %in% c("yes","y","true"),
                               "reported, modelled", "reported"), "not reported"),
    sex = ifelse(tolower(trimws(sex_reported)) %in% c("yes","y","true"),
                 ifelse(tolower(trimws(sex_included)) %in% c("yes","y","true"),
                        "reported, modelled", "reported"), "not reported"),
    effect_size_metric = effect_size_metrics,
    n_models = n_models_included,
    n_effect_sizes = k_effect_sizes,
    n_study_clusters = n_study_clusters,
    n_sign_reversals = n_sign_reversals_reported)

stopifnot(nrow(S1) == 28L, sum(S1$n_effect_sizes) == 5740L, sum(S1$n_models) == 48L)
unlink(file.path(SUP, "TableS2_evidence_base.csv"))
readr::write_csv(S1, file.path(SUP, "TableS1_evidence_base.csv"))
message(sprintf("Table S1: %d rows, %d models, %d effect sizes -> supplementary/TableS1_evidence_base.csv",
        nrow(S1), sum(S1$n_models), sum(S1$n_effect_sizes)))

# --- Tables S3 and S4: the two influence analyses -----------------------------
# Shared display vocabulary, so a reader moving between S1, S3 and S4 meets the same
# words for the same things. `.default = NA_character_` for the same reason as above: an
# unmapped internal label must fail the gate rather than leak into a supplementary table.
spec_label <- function(x) dplyr::recode(x,
  uncorrected_beta0       = "Uncorrected pooled mean",
  yang2023_gated_beta0_c3 = "Yang 2023 bias-corrected",
  yang2024_FE_VCV         = "Yang 2024 bias-robust (FE + VCV)",
  yang2024_UWLS           = "Yang 2024 bias-robust (UWLS)",
  .default = NA_character_)
weight_label <- function(x) dplyr::recode(x,
  k_effect_sizes          = "By effect-size count",
  equal_per_meta_analysis = "Each meta-analysis weighted equally",
  .default = NA_character_)
metric_label <- function(x) dplyr::recode(x,
  power = "Statistical power", type_M = "Type M error", type_S = "Type S error",
  .default = NA_character_)

S3 <- readr::read_csv(file.path(REV_OUT, "loo_influence.csv"), show_col_types = FALSE) |>
  dplyr::transmute(
    specification = spec_label(effect_estimator),
    weighting     = weight_label(weighting),
    metric        = metric_label(metric),
    dropped_model = sub("[.]csv$", "", dropped_MA_model),
    n_effect_sizes_contributed  = dropped_k,
    pct_of_effect_sizes         = sprintf("%.1f", dropped_pct_of_k),
    summary_all_48_models = fmt(summary_all_48, metric),
    summary_without       = fmt(summary_without, metric),
    change_pct            = sprintf("%+.1f", pct_change),
    influence_rank) |>
  dplyr::arrange(weighting, specification, metric, influence_rank)

S4 <- readr::read_csv(file.path(REV_OUT, "leave_one_paper_out.csv"), show_col_types = FALSE) |>
  dplyr::transmute(
    specification = spec_label(effect_estimator),
    weighting     = weight_label(weighting),
    metric        = metric_label(metric),
    dropped_paper = dropped_source_paper,
    n_models_dropped, n_effect_sizes_dropped,
    pct_of_effect_sizes = sprintf("%.1f", pct_of_effect_sizes_dropped),
    summary_all_28_papers = fmt(summary_all_28_papers, metric),
    summary_without       = fmt(geometric_mean, metric),
    change_pct            = sprintf("%+.1f", pct_change),
    influence_rank) |>
  dplyr::arrange(weighting, specification, metric, influence_rank)

# gates: no unmapped label may reach a supplementary table, and the row counts must be
# the full grids rather than a silently filtered subset
for (nm in c("S3", "S4")) {
  tb <- get(nm)
  leaked <- vapply(tb, function(col) any(is.na(col)), logical(1))
  if (any(leaked)) stop(nm, ": unmapped or missing values in ",
                        paste(names(tb)[leaked], collapse = ", "))
}
stopifnot(nrow(S3) == 48L * 4L * 3L * 2L, nrow(S4) == 28L * 4L * 3L * 2L)

readr::write_csv(S3, file.path(SUP, "TableS3_leave_one_model_out.csv"))
readr::write_csv(S4, file.path(SUP, "TableS4_leave_one_paper_out.csv"))
message(sprintf("Table S3: %d rows -> supplementary/TableS3_leave_one_model_out.csv", nrow(S3)))
message(sprintf("Table S4: %d rows -> supplementary/TableS4_leave_one_paper_out.csv", nrow(S4)))

# the number the main text quotes must be the top row of its block, so that a reader
# following the sentence lands on it
top <- S3 |>
  dplyr::filter(specification == "Yang 2024 bias-robust (FE + VCV)",
                metric == "Statistical power",
                weighting == "By effect-size count", influence_rank == 1L)
stopifnot(nrow(top) == 1L, top$dropped_model == "MA09", top$change_pct == "+44.7")
message(sprintf("  main-text check: largest single-model influence is %s at %s%% (rank 1 of 48)",
        top$dropped_model, top$change_pct))

# --- what stays as CSV rather than becoming a supplementary table -------------
keep <- setdiff(basename(list.files(REV_OUT, pattern = "[.]csv$")), character(0))
message("\nremaining as archived CSV (not supplementary tables): ",
        paste(sub("[.]csv$", "", keep), collapse = ", "))
