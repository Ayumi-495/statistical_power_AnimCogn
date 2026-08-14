# R/revision/13_table_metadata.R ------------------------------------------------
# Step 13: captions and column metadata for the supplementary tables, and a
# file-level index of everything that goes into the archived deposit.
#
# Biology Open accepts no supplementary text other than captions, so each caption has
# to be self-contained: a reader must be able to interpret the table without the
# Methods in front of them.
#
# The column dictionary follows the convention already used by the sibling systematic
# map (`data/metadata.xlsx`, sheet `variable`). It is generated here rather than typed
# by hand so that it cannot drift out of step with the tables: the script reads the
# actual column names and fails if a description is missing or refers to a column that
# no longer exists.
#
# Metadata is not decoration. Reviewer 1 asked for "data, code, and meta-data in a
# repository with a permanent DOI ... and with a licence that describes re-use". The
# dictionary is the meta-data half of that request.

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 13: captions and metadata ==")
SUP <- file.path(REV_OUT, "supplementary")

# --- captions -----------------------------------------------------------------
cap_S1 <- paste0(
  "**Table S1. Statistical power, Type M error and Type S error under each assumed ",
  "underlying effect.** Part A gives the results reported in the main text, at the ",
  "primary-study level (5,740 effect-size estimates) and at the meta-analysis level ",
  "(48 models). Meta-analysis-level summaries are shown both weighted by the number of ",
  "effect-size estimates contributed by each meta-analysis and with each meta-analysis ",
  "weighted equally; the two answer different questions and are not alternative ",
  "estimates of one quantity. Part B gives sensitivity analyses at the primary-study ",
  "level: an optimistic assumed effect, defined as the limit of the 95% confidence ",
  "interval of the uncorrected pooled mean lying farther from zero, and externally ",
  "specified assumed effects. External values are reported separately for each ",
  "effect-size metric because the three sets are not commensurable: Cohen's ",
  "conventions for the standardised mean difference and for Fisher's Zr correspond to ",
  "different underlying magnitudes, and the log response ratio values are reference ",
  "points rather than conventions. Summary estimates are back-transformed intercepts ",
  "from models fitted on the log scale, with 95% confidence intervals; raw medians and ",
  "interquartile ranges are given alongside. Type S is fitted with an offset of 0.025 ",
  "because it can be exactly zero, and rows where the summary is sensitive to that ",
  "offset are flagged. All metrics use a common normal-theory two-sided threshold at ",
  "alpha = 0.05 rather than the reference distribution of each fitted model.")

cap_S2 <- paste0(
  "**Table S2. Characteristics of the 28 meta-analytical papers included in the ",
  "quantitative reanalysis.** Species, cognitive domain, task and manipulation are ",
  "taken from the systematic map that supplied the corpus. Cognitive domain follows the ",
  "categories assigned in that map; a paper spanning more than one category is listed ",
  "under each. 'Life stage' and 'sex' record whether the meta-analysis reported the ",
  "variable and whether it was also included as a moderator. The last four columns ",
  "describe what each paper contributes to the present analysis: the number of ",
  "meta-analysis models, effect-size estimates and study clusters, and the number of ",
  "models whose bias-corrected mean had the opposite sign to the uncorrected mean. ",
  "Task is recorded as 'unclear' where the systematic map could not identify a ",
  "specific paradigm, and species were not recorded for one paper.")

# --- column dictionary --------------------------------------------------------
dict <- tibble::tribble(
  ~table, ~column, ~description, ~units_or_values,
  "Table S1", "part", "Whether the row belongs to the results reported in the main text or to a sensitivity analysis", "A. Reported results; B. Sensitivity to the assumed effect",
  "Table S1", "level", "The unit the metric is computed for. At the primary-study level each effect-size estimate is evaluated against its own sampling standard error; at the meta-analysis level each model is evaluated against the standard error of its pooled estimate. In Part B, the effect-size metric the external value applies to", "Primary-study level; Meta-analysis level; All metrics; SMD; Zr; lnRR",
  "Table S1", "assumed_effect", "The value assumed to be the underlying true effect when computing the metric", "text",
  "Table S1", "weighting", "How meta-analyses are weighted when summarising across them. Not applicable at the primary-study level", "By effect-size count; Equal per meta-analysis; -",
  "Table S1", "metric", "Which design-analysis quantity the row reports", "Statistical power; Type M error; Type S error",
  "Table S1", "n_unit", "Number of units the summary is computed over", "count",
  "Table S1", "summary_estimate", "Back-transformed intercept of a model fitted on the log scale: a weighted linear model at the meta-analysis level, a linear mixed model with a study random effect at the primary-study level", "probability (power, Type S) or ratio (Type M)",
  "Table S1", "ci", "95% confidence interval for summary_estimate, back-transformed on the same scale", "interval",
  "Table S1", "raw_median_iqr", "Median of the underlying values with the first and third quartiles, computed directly and without any model or offset", "median [Q1, Q3]",
  "Table S1", "offset_note", "Flag for rows where the Type S values are mostly below the 0.025 offset used to allow the log transformation, so that summary_estimate reflects the offset more than the data and raw_median_iqr should be preferred", "text or empty",

  "Table S2", "paper", "Identifier of the meta-analytical paper in the systematic map that supplied the corpus", "MA01-MA47",
  "Table S2", "study", "First author and year of the meta-analytical paper", "text",
  "Table S2", "year", "Year of publication of the meta-analytical paper", "year",
  "Table S2", "species", "Species studied, from the systematic map. Papers covering more than three species are summarised as the first two followed by a count", "text",
  "Table S2", "cognitive_domain", "Cognitive domain assigned in the systematic map. A paper spanning more than one domain is listed under each", "Learning; Memory; Decision-making; Other",
  "Table S2", "task", "Behavioural task or paradigm, from the systematic map, truncated for display. 'unclear' where no specific paradigm could be identified", "text",
  "Table S2", "manipulation", "Category of experimental manipulation, from the systematic map", "Drug; Environmental; Nutrient; Others; Not applicable; Unclear",
  "Table S2", "life_stage", "Whether the meta-analysis reported life stage, and whether it also included life stage as a moderator", "not reported; reported; reported, modelled",
  "Table S2", "sex", "Whether the meta-analysis reported sex, and whether it also included sex as a moderator", "not reported; reported; reported, modelled",
  "Table S2", "effect_size_metric", "Effect-size metric or metrics used by the models taken from this paper", "SMD; lnRR; Zr",
  "Table S2", "n_models", "Number of meta-analysis models this paper contributes to the present analysis", "count",
  "Table S2", "n_effect_sizes", "Number of effect-size estimates contributed, summed over that paper's models", "count",
  "Table S2", "n_study_clusters", "Number of study clusters contributed, summed over that paper's models. Study identifiers are defined within each meta-analysis and are not harmonised across papers, so this is not a count of distinct primary publications", "count",
  "Table S2", "n_sign_reversals", "Number of this paper's models whose bias-corrected mean has the opposite sign to the uncorrected pooled mean", "count"
)

# --- gate: the dictionary must match the tables exactly ------------------------
for (nm in c("Table S1", "Table S2")) {
  f <- if (nm == "Table S1") "TableS1_reported_metrics.csv" else "TableS2_evidence_base.csv"
  actual <- names(readr::read_csv(file.path(SUP, f), n_max = 0, show_col_types = FALSE))
  described <- dict$column[dict$table == nm]
  missing <- setdiff(actual, described); extra <- setdiff(described, actual)
  if (length(missing)) stop(nm, ": columns present but undescribed: ", paste(missing, collapse = ", "))
  if (length(extra))   stop(nm, ": described but not present: ", paste(extra, collapse = ", "))
  message(sprintf("  %s: %d columns, all described", nm, length(actual)))
}

readr::write_csv(dict, file.path(SUP, "metadata_columns.csv"))

# --- file-level index for the archived deposit --------------------------------
files <- tibble::tribble(
  ~file, ~contents,
  "TableS1_reported_metrics.csv", "Supplementary Table S1. Summary metrics under each assumed effect.",
  "TableS2_evidence_base.csv", "Supplementary Table S2. Characteristics of the 28 included papers.",
  "per_meta_analysis_estimates.csv", "One row per meta-analysis model: uncorrected and bias-corrected means, standard errors and intervals, the bias-robust estimate with its cluster-robust standard error and Satterthwaite degrees of freedom, weighting diagnostics and sign-reversal flags.",
  "model_level_metrics.csv", "Power, Type M and Type S for every model under every assumed effect, with the assumed effect and standard error that produced each.",
  "assumed_effect_scenarios.csv", "Optimistic and externally specified assumed-effect scenarios at both levels.",
  "primary_level_sensitivity.csv", "Primary-study-level summaries with provenance columns.",
  "meta_analysis_level_sensitivity.csv", "Meta-analysis-level summaries with provenance columns, both weightings.",
  "loo_influence.csv", "Leave-one-model-out influence on the meta-analysis-level summaries: all 48 models by 4 specifications by 3 metrics.",
  "reversal_counts.csv", "Sign-reversal counts under each bias-correction approach.",
  "rho_sensitivity.csv", "Sensitivity of the bias-robust analysis to the assumed within-study sampling correlation.",
  "yang2024_reference_validation.csv", "Validation of the bias-robust implementation against the published worked example of the method.",
  "evidence_base_characteristics.csv", "Full-width source table behind Table S2, before columns were trimmed for display."
)
have <- c(list.files(REV_OUT, pattern = "[.]csv$"), list.files(SUP, pattern = "[.]csv$"))
undocumented <- setdiff(setdiff(have, c("metadata_columns.csv", "metadata_files.csv")), files$file)
if (length(undocumented)) stop("result files with no entry in the file index: ",
                               paste(undocumented, collapse = ", "))
files$n_rows <- vapply(files$file, function(f) {
  p <- if (file.exists(file.path(SUP, f))) file.path(SUP, f) else file.path(REV_OUT, f)
  nrow(readr::read_csv(p, show_col_types = FALSE)) }, numeric(1))
readr::write_csv(files, file.path(SUP, "metadata_files.csv"))
message(sprintf("  file index: %d files, all documented", nrow(files)))

writeLines(c("# Supplementary table captions", "", cap_S1, "", cap_S2, "",
             "Column definitions are in metadata_columns.csv; the file index for the archived",
             "deposit is in metadata_files.csv."),
           file.path(SUP, "captions.md"))
message("wrote supplementary/{captions.md, metadata_columns.csv, metadata_files.csv}")
