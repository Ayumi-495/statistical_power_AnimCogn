# revision/R/analyse/13_table_metadata.R ------------------------------------------------
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

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))

message("== 13: captions and metadata ==")
SUP <- file.path(REV_OUT, "supplementary")

# --- captions -----------------------------------------------------------------
cap_S2 <- paste0(
  "**Table S2. Statistical power, Type M error and Type S error under each assumed ",
  "underlying effect.** Part A gives the results reported in the main text, at the ",
  "primary-study level (5,740 effect-size estimates) and at the meta-analysis level ",
  "(48 models). At the primary-study level the reported result is the row that fits ",
  "study cluster as a random effect. Because that model shrinks large clusters, it ",
  "weights study clusters almost equally rather than weighting each effect-size ",
  "estimate equally; the other two primary-study-level rows give each meta-analysis ",
  "equal influence and are sensitivity analyses. Meta-analysis-level summaries are shown both weighted by the number of ",
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
  "offset are flagged. Confidence limits are back-transformed from models fitted on the ",
  "log scale and are therefore unbounded; limits falling outside the range a metric can ",
  "take are constrained to the bound and marked with an asterisk, as in the original ",
  "supplement. Power lies in [0, 1], Type S error in [0, 0.5], and Type M error is at ",
  "least 1. No point estimate required constraining. All metrics use a common ",
  "normal-theory two-sided threshold at alpha = 0.05 rather than the reference ",
  "distribution of each fitted model.")

cap_S1 <- paste0(
  "**Table S1. Characteristics of the 28 meta-analytical papers included in the ",
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

cap_S3 <- paste0(
  "**Table S3. Leave-one-meta-analysis-model-out influence on the meta-analysis-level ",
  "summaries.** Each of the 48 meta-analytic models is removed in turn and the summary ",
  "recomputed from the remaining 47, for each of the four assumed effects, each of the ",
  "three metrics, and both weightings. `change_pct` is the percentage change in the ",
  "summary caused by removing that model, and `influence_rank` orders the 48 models by ",
  "the absolute size of that change within each specification, metric and weighting; ",
  "rows are sorted so that the most influential model appears first in each block. ",
  "Single-model leverage is largely a property of the weighting rather than of any ",
  "model: weighting by effect-size count gives one model 22.6% of the weight and it ",
  "moves the bias-robust power summary by 44.7%, whereas under equal weighting no model ",
  "moves any power summary by more than 5.1%. The contrast is smaller but in the same ",
  "direction for Type M error (34.3% against 14.3%) and for Type S error (78.0% against ",
  "34.9%), which are more sensitive to individual models because both diverge as the ",
  "assumed effect approaches zero. Nothing is excluded from any reported analysis; ",
  "this is an influence diagnostic. Companion to Table S4, which removes a source paper ",
  "at a time.")

cap_S4 <- paste0(
  "**Table S4. Leave-one-source-paper-out influence on the meta-analysis-level ",
  "summaries.** Each of the 28 source meta-analytical papers is removed in turn, ",
  "together with every meta-analytic model it contributed, and the summary recomputed ",
  "from the models that remain. Twelve papers contribute more than one model, so this ",
  "asks a different question from Table S3: how much does the summary depend on any ",
  "single published paper, rather than on any single fitted model. A paper is the unit ",
  "of publication and of any shared authorship, laboratory or analytical convention. ",
  "Columns and sorting follow Table S3. Nothing is excluded from any reported analysis.")

# --- column dictionary --------------------------------------------------------
dict <- tibble::tribble(
  ~table, ~column, ~description, ~units_or_values,
  "Table S2", "part", "Whether the row belongs to the results reported in the main text or to a sensitivity analysis", "A. Reported results; B. Sensitivity to the assumed effect (primary-study level); C. Externally specified assumed effects (meta-analysis level); D. Leave-one-cluster-out (primary-study level)",
  "Table S2", "level", "The unit the metric is computed for. At the primary-study level each effect-size estimate is evaluated against its own sampling standard error; at the meta-analysis level each model is evaluated against the standard error of its pooled estimate. In Part B, the effect-size metric the external value applies to", "Primary-study level; Meta-analysis level; All metrics; SMD; Zr; lnRR",
  "Table S2", "assumed_effect", "The value assumed to be the underlying true effect when computing the metric. In Part D, whether that value was estimated with the effect size's own study cluster included or with it held out", "text",
  "Table S2", "weighting", "How the individual values are combined into the summary. At the meta-analysis level, whether models are weighted by the number of effect-size estimates they contribute or equally. At the primary-study level, which of three estimands the row reports: a random intercept for study cluster, which is the reported result and which weights study clusters almost equally rather than weighting each effect-size estimate equally; each meta-analysis weighted equally; or meta-analysis identity fitted as a second random effect", "By effect-size count; Equal per meta-analysis; Study cluster as a random effect; Meta-analysis as a random effect; -",
  "Table S2", "metric", "Which design-analysis quantity the row reports", "Statistical power; Type M error; Type S error",
  "Table S2", "n_unit", "Number of units the summary is computed over", "count",
  "Table S2", "summary_estimate", "Back-transformed intercept of a model fitted on the log scale: a weighted linear model at the meta-analysis level, a linear mixed model with a study random effect at the primary-study level", "probability (power, Type S) or ratio (Type M)",
  "Table S2", "ci", "95% confidence interval for summary_estimate, back-transformed on the same scale. Limits falling outside the range the metric can take are constrained to the bound and marked with an asterisk; power lies in [0, 1], Type S error in [0, 0.5] and Type M error is at least 1. No point estimate required constraining", "interval",
  "Table S2", "raw_median_iqr", "Median of the underlying values with the first and third quartiles, computed directly and without any model or offset", "median [Q1, Q3]",
  "Table S2", "n_constrained_limits", "How many of the two interval limits in this row were constrained to the metric's bound", "0, 1 or 2",
  "Table S2", "offset_note", "Flag for rows where the Type S values are mostly below the 0.025 offset used to allow the log transformation, so that summary_estimate reflects the offset more than the data and raw_median_iqr should be preferred", "text or empty",

  "Table S1", "paper", "Identifier of the meta-analytical paper in the systematic map that supplied the corpus", "MA01-MA47",
  "Table S1", "study", "First author and year of the meta-analytical paper", "text",
  "Table S1", "year", "Year of publication of the meta-analytical paper", "year",
  "Table S1", "species", "Species studied, from the systematic map. Papers covering more than three species are summarised as the first two followed by a count", "text",
  "Table S1", "cognitive_domain", "Cognitive domain assigned in the systematic map. A paper spanning more than one domain is listed under each", "Learning; Memory; Decision-making; Other",
  "Table S1", "task", "Behavioural task or paradigm, from the systematic map, truncated for display. 'unclear' where no specific paradigm could be identified", "text",
  "Table S1", "manipulation", "Category of experimental manipulation, from the systematic map", "Drug; Environmental; Nutrient; Others; Not applicable; Unclear",
  "Table S1", "life_stage", "Whether the meta-analysis reported life stage, and whether it also included life stage as a moderator", "not reported; reported; reported, modelled",
  "Table S1", "sex", "Whether the meta-analysis reported sex, and whether it also included sex as a moderator", "not reported; reported; reported, modelled",
  "Table S1", "effect_size_metric", "Effect-size metric or metrics used by the models taken from this paper", "SMD; lnRR; Zr",
  "Table S1", "n_models", "Number of meta-analysis models this paper contributes to the present analysis", "count",
  "Table S1", "n_effect_sizes", "Number of effect-size estimates contributed, summed over that paper's models", "count",
  "Table S1", "n_study_clusters", "Number of study clusters contributed, summed over that paper's models. Study identifiers are defined within each meta-analysis and are not harmonised across papers, so this is not a count of distinct primary publications", "count",
  "Table S1", "n_sign_reversals", "Number of this paper's models whose bias-corrected mean has the opposite sign to the uncorrected pooled mean", "count",

  "Table S3", "specification", "Which assumed underlying effect the summary was computed against", "Uncorrected pooled mean; Yang 2023 bias-corrected; Yang 2024 bias-robust (FE + VCV); Yang 2024 bias-robust (UWLS)",
  "Table S3", "weighting", "How the remaining models are combined into the summary", "By effect-size count; Each meta-analysis weighted equally",
  "Table S3", "metric", "Which design-analysis quantity the row reports", "Statistical power; Type M error; Type S error",
  "Table S3", "dropped_model", "The meta-analytic model removed in this row", "MA01_01-MA47",
  "Table S3", "n_effect_sizes_contributed", "Number of effect-size estimates that model contributes", "count",
  "Table S3", "pct_of_effect_sizes", "That model's share of all 5,740 effect-size estimates, which is its share of the weight under weighting by effect-size count", "per cent",
  "Table S3", "summary_all_48_models", "The summary computed from all 48 models, repeated on every row of the block for comparison", "probability (power, Type S) or ratio (Type M)",
  "Table S3", "summary_without", "The summary recomputed from the 47 models that remain", "probability (power, Type S) or ratio (Type M)",
  "Table S3", "change_pct", "Percentage change from summary_all_48_models to summary_without", "per cent, signed",
  "Table S3", "influence_rank", "Rank of this model by the absolute size of change_pct, within this specification, metric and weighting; 1 is the most influential", "1-48",

  "Table S4", "specification", "Which assumed underlying effect the summary was computed against", "Uncorrected pooled mean; Yang 2023 bias-corrected; Yang 2024 bias-robust (FE + VCV); Yang 2024 bias-robust (UWLS)",
  "Table S4", "weighting", "How the remaining models are combined into the summary", "By effect-size count; Each meta-analysis weighted equally",
  "Table S4", "metric", "Which design-analysis quantity the row reports", "Statistical power; Type M error; Type S error",
  "Table S4", "dropped_paper", "The source meta-analytical paper removed in this row, together with every model it contributed", "MA01-MA47",
  "Table S4", "n_models_dropped", "How many of the 48 meta-analytic models that paper contributed, and therefore how many were removed", "count",
  "Table S4", "n_effect_sizes_dropped", "How many effect-size estimates were removed with it", "count",
  "Table S4", "pct_of_effect_sizes", "Those estimates as a share of all 5,740", "per cent",
  "Table S4", "summary_all_28_papers", "The summary computed from all 48 models, repeated on every row of the block for comparison", "probability (power, Type S) or ratio (Type M)",
  "Table S4", "summary_without", "The summary recomputed from the models that remain after the paper is removed", "probability (power, Type S) or ratio (Type M)",
  "Table S4", "change_pct", "Percentage change from summary_all_28_papers to summary_without", "per cent, signed",
  "Table S4", "influence_rank", "Rank of this paper by the absolute size of change_pct, within this specification, metric and weighting; 1 is the most influential", "1-28"
)

# --- gate: the dictionary must match the tables exactly ------------------------
TABLE_FILE <- c("Table S2" = "TableS2_reported_metrics.csv",
                "Table S1" = "TableS1_evidence_base.csv",
                "Table S3" = "TableS3_leave_one_model_out.csv",
                "Table S4" = "TableS4_leave_one_paper_out.csv")
for (nm in names(TABLE_FILE)) {
  f <- TABLE_FILE[[nm]]
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
  "TableS2_reported_metrics.csv", "Supplementary Table S2. Summary metrics under each assumed effect.",
  "TableS1_evidence_base.csv", "Supplementary Table S1. Characteristics of the 28 included papers.",
  "TableS3_leave_one_model_out.csv", "Supplementary Table S3. Leave-one-model-out influence, display version of loo_influence.csv.",
  "TableS4_leave_one_paper_out.csv", "Supplementary Table S4. Leave-one-source-paper-out influence, display version of leave_one_paper_out.csv.",
  "per_meta_analysis_estimates.csv", "One row per meta-analysis model: uncorrected and bias-corrected means, standard errors and intervals, the bias-robust estimate with its cluster-robust standard error and Satterthwaite degrees of freedom, weighting diagnostics and sign-reversal flags.",
  "model_level_metrics.csv", "Power, Type M and Type S for every model under every assumed effect, with the assumed effect and standard error that produced each.",
  "assumed_effect_scenarios.csv", "Optimistic and externally specified assumed-effect scenarios at both levels.",
  "primary_level_sensitivity.csv", "Primary-study-level summaries with provenance columns.",
  "meta_analysis_level_sensitivity.csv", "Meta-analysis-level summaries with provenance columns, both weightings.",
  "loo_influence.csv", "Leave-one-model-out influence on the meta-analysis-level summaries: all 48 models by 4 specifications by 3 metrics by both weightings. Companion to leave_one_paper_out.csv, which removes one source paper at a time together with every model it contributed. Single-model leverage is largely a property of the weighting: equal weighting roughly halves the largest single-model influence for every metric, and for statistical power removes it entirely (44.7% by effect-size count against 5.1% equally weighted). Type M and Type S remain sensitive under either weighting because both diverge as the assumed effect approaches zero.",
  "leave_one_paper_out.csv", "Leave-one-source-paper-out influence on the meta-analysis-level summaries: each of the 28 source papers removed in turn together with every model it contributed, for 4 specifications, 3 metrics and both weightings. Companion to loo_influence.csv, which removes one model at a time.",
  "leave_one_cluster_out.csv", "Leave-one-cluster-out at the primary-study level: the assumed effect recomputed with each study cluster removed, aggregated under the adopted clustering definition and under both primary-study-level estimands.",
  "ma_level_uncertainty.csv", "Comparison of three ways of estimating the uncertainty of the meta-analysis-level summary, holding the point estimate fixed: the model-based interval currently reported, a cluster-robust CR2 interval clustered by source paper with Satterthwaite degrees of freedom, and a nonparametric cluster bootstrap over the 28 source papers (percentile and BCa). A 48-model bootstrap is included as a labelled comparison and is not a candidate.",
  "ma_level_paired_contrasts.csv", "Paired paper-level bootstrap contrasts between every pair of the four assumed-effect specifications, for all three metrics. Both terms of each contrast are recomputed inside the same resampled set of source papers, so the pairing between them is preserved; comparing the marginal intervals instead would discard it. Same 20,000 resamples and seed as ma_level_uncertainty.csv.",
  "verification_audit.csv", "One row per verification check: what was compared, by which two routes, the criterion it was held to and whether it passed. Includes the check that the reported primary-study-level summary is the estimand it is labelled as.",
  "reversal_counts.csv", "Sign-reversal counts under each bias-correction approach.",
  "rho_sensitivity.csv", "Sensitivity of the bias-robust analysis to the assumed within-study sampling correlation.",
  "yang2024_reference_validation.csv", "Validation of the bias-robust implementation against the published worked example of the method.",
  "evidence_base_characteristics.csv", "Full-width source table behind Table S1, before columns were trimmed for display."
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
             cap_S3, "", cap_S4, "",
             "Column definitions are in metadata_columns.csv; the file index for the archived",
             "deposit is in metadata_files.csv."),
           file.path(SUP, "captions.md"))
message("wrote supplementary/{captions.md, metadata_columns.csv, metadata_files.csv}")
