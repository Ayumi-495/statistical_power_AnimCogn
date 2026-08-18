# revision/R/analyse/16_export_scenario_inputs.R ----------------------------------------
# Step 16a: export the raw inputs the assumed-effect scenarios are built from, so
# that `17_verify_scenarios.py` can recompute the whole of
# `revision/results/assumed_effect_scenarios.csv` from scratch in another language.
#
# WHY. Every headline number in this paper carries two derivations by different
# routes. The scenario rows were migrated from the audit pipeline and were the only
# canonical table still marked `single_derivation_...`, and about ten of them are
# about to be quoted in the Results. This closes that gap.
#
# The export is deliberately minimal: assumed effects, standard errors, cluster
# labels and effect-size metrics. Nothing derived. The Python side re-derives the
# scenario definitions, the three closed-form metrics, the REML random-intercept
# summary and the weighted least-squares summary independently, so an error shared
# between the two would have to be an error in these inputs, which come from the
# fitted models that `01_reproduce_original_analysis.R` already verifies.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))

message("== 16: export scenario inputs for independent verification ==")
S <- readRDS(file.path(REV_TMP, "original_estimates.rds"))
o <- S$original; L <- all_datasets(S$dat)
stopifnot(identical(names(L), o$MA_model), nrow(o) == 48L)

# meta-analysis level: one row per model
ma <- tibble::tibble(
  MA_model = o$MA_model, effect_size_type = o$effect_size_type, k = o$k,
  beta0 = o$beta0, se_beta0 = o$se_beta0,
  ci_lb_beta0 = o$ci_lb_beta0, ci_ub_beta0 = o$ci_ub_beta0)

# primary-study level: one row per effect-size estimate
prim <- purrr::list_rbind(lapply(seq_along(L), function(i) {
  x <- L[[i]]
  tibble::tibble(MA_model = o$MA_model[i], effect_size_type = o$effect_size_type[i],
                 study_ID = as.character(x$study_ID), sei = x$sei)
}))
stopifnot(nrow(prim) == 5740L)
prim$cluster <- namespaced_study_id(prim$MA_model, prim$study_ID)

readr::write_csv(ma,   file.path(REV_TMP, "verify_inputs_ma.csv"))
readr::write_csv(prim, file.path(REV_TMP, "verify_inputs_primary.csv"))
message(sprintf("wrote verify_inputs_ma.csv (%d models) and verify_inputs_primary.csv (%d estimates, %d clusters)",
        nrow(ma), nrow(prim), dplyr::n_distinct(prim$cluster)))
message("now run: python3 revision/R/analyse/17_verify_scenarios.py")
