# R/revision/12_supplementary_tables.R ------------------------------------------
# Step 12: the two supplementary tables that go into the manuscript.
#
# Only two. Everything else stays as CSV in the repository and in the archived,
# DOI-bearing deposit, which is also the answer to the FAIR / DOI / licence request:
# large tables belong in a citable archive rather than compressed into a PDF.
#
#   Table S1  the reported metrics, both levels, every assumed effect, with the
#             summary, its interval, and the raw quantiles. This is what a reader
#             needs in order to check the claims in the text.
#   Table S2  the characteristics of the 28 included meta-analytical papers. Four
#             reviewer comments depend on it (R1C17, R2C11, R2C15, R2C17).
#
# Both are written as CSV, for pasting into the manuscript, and echoed as markdown
# so they can be read without opening a spreadsheet.
#
# TYPE S IS REPORTED BOTH WAYS AND THE TABLE SAYS SO. The summary is fitted on
# log(Type S + 0.025) and the offset dominates when the values are far below it, so
# every row carries the raw median and quartiles beside the model-based value. Which
# of the two becomes the headline is still with the PI; the table does not prejudge it.

source(here::here("R", "revision", "00_revision_functions.R"))

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

# --- Table S1 -----------------------------------------------------------------
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
      unweighted_per_effect_size  = "Each effect-size estimate weighted equally",
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

S1 <- dplyr::bind_rows(core, scen, scen_ma) |>
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
bad <- dplyr::bind_rows(core, scen, scen_ma) |>
  dplyr::mutate(metric = dplyr::recode(metric, power = "Statistical power",
                                       type_M = "Type M error", type_S = "Type S error")) |>
  dplyr::filter(mapply(was_constrained, geometric_mean, metric))
if (nrow(bad)) stop(sprintf("%d point estimates fall outside the metric's range", nrow(bad)))

# No internal label may reach the table.
if (any(is.na(S1$weighting)))
  stop("unmapped `weighting` value reached Table S1; add it to the recode above")
leak <- grep("_", c(S1$weighting, S1$level, S1$assumed_effect, S1$metric), value = TRUE)
if (length(leak)) stop("internal shorthand in Table S1: ", paste(unique(leak), collapse = ", "))

readr::write_csv(S1, file.path(SUP, "TableS1_reported_metrics.csv"))
message(sprintf("Table S1: %d rows -> supplementary/TableS1_reported_metrics.csv", nrow(S1)))

# --- Table S2 -----------------------------------------------------------------
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

S2 <- readr::read_csv(file.path(REV_OUT, "evidence_base_characteristics.csv"),
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

stopifnot(nrow(S2) == 28L, sum(S2$n_effect_sizes) == 5740L, sum(S2$n_models) == 48L)
readr::write_csv(S2, file.path(SUP, "TableS2_evidence_base.csv"))
message(sprintf("Table S2: %d rows, %d models, %d effect sizes -> supplementary/TableS2_evidence_base.csv",
        nrow(S2), sum(S2$n_models), sum(S2$n_effect_sizes)))

# --- markdown echo, so both can be read without a spreadsheet -----------------
md <- function(d) {
  hdr <- paste0("| ", paste(names(d), collapse = " | "), " |")
  sep <- paste0("|", paste(rep(" --- ", ncol(d)), collapse = "|"), "|")
  rows <- apply(d, 1, function(r) paste0("| ", paste(ifelse(is.na(r), "", r), collapse = " | "), " |"))
  c(hdr, sep, rows)
}
writeLines(c("# Table S1. Reported metrics", "",
             md(dplyr::filter(S1, part == "A. Reported results") |> dplyr::select(-part)), "",
             "# Table S1, part B. Sensitivity to the assumed effect (primary-study level)", "",
             md(dplyr::filter(S1, startsWith(part, "B.")) |> dplyr::select(-part)), "",
             "# Table S1, part C. External assumed effects (meta-analysis level)", "",
             md(dplyr::filter(S1, startsWith(part, "C.")) |> dplyr::select(-part)), "",
             "# Table S2. Characteristics of the 28 included meta-analytical papers", "",
             md(S2)),
           file.path(SUP, "supplementary_tables.md"))
message("markdown echo -> supplementary/supplementary_tables.md")

# --- what stays as CSV rather than becoming a supplementary table -------------
keep <- setdiff(basename(list.files(REV_OUT, pattern = "[.]csv$")), character(0))
message("\nremaining as archived CSV (not supplementary tables): ",
        paste(sub("[.]csv$", "", keep), collapse = ", "))
