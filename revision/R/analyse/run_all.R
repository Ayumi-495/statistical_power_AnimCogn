# revision/R/analyse/run_all.R ---------------------------------------------------------
# Runs the revision workflow in order. Each script is also runnable on its own,
# provided the earlier ones have been run (they exchange fitted objects through
# revision/results/intermediate/, which is gitignored and regenerable).
#
# Requires the `clubSandwich` package: the Yang-2024 sensitivity analysis is
# specified as CR2 with Satterthwaite degrees of freedom, and 00_revision_functions.R
# stops rather than falling back to the hand-written CR1 sandwich, which is a
# different specification.
#
# NOT included here: 06_validate_yang2024_reference.R, which validates the CR2
# implementation against the authors' published worked example. It needs to download
# a data file from another project's repository, and this script must work offline.
# Run it separately after any change to the CR2 code path:
#   Rscript revision/R/analyse/06_validate_yang2024_reference.R
#
# ALSO NOT included here: 17_verify_scenarios.py, the independent Python
# re-derivation of every assumed-effect scenario row. Step 16 below writes the
# inputs it needs; run it after any change to 09:
#   python3 revision/R/analyse/17_verify_scenarios.py
#
# Steps 07 and 08 must follow 05: 07 checks its baselines against the canonical
# summaries, and 08 reads them to draw the reference rules.
#
# Total runtime is a few minutes, dominated by the rma.mv fits in steps 01 and 03.
# TWO COMMITTED INPUTS THIS SCRIPT CONSUMES BUT CANNOT REGENERATE. Both come from
# scripts that need something outside this repository, so they are committed as results
# and `run_all.R` reads them:
#
#   evidence_base_characteristics.csv   from 10_evidence_base_table.R, which reads the
#                                       sibling systematic map. Step 12 needs it.
#   yang2024_reference_validation.csv   from 06_validate_yang2024_reference.R, which
#                                       downloads the method authors' example data.
#
# Deleting either and running this script offline will fail at step 12 or 13. That is
# deliberate - a missing input should stop the run, not silently produce a shorter table.
for (f in c("01_reproduce_original_analysis.R", "02_overshoot_diagnostics.R",
            "03_yang2024_bias_robust.R", "04_revision_sensitivity_summaries.R",
            "05_make_revision_tables.R", "07_influence_loo.R",
            "08_model_level_figure.R", "09_assumed_effect_scenarios.R",
            "11_main_figure.R", "14_leave_one_cluster_out.R",
            # 12 builds Tables S3 and S4 out of 07's and 15's output, so it has to
            # follow both. It used to run before 14 and 15, when it made only S1 and S2.
            "15_leave_one_paper_out.R", "12_supplementary_tables.R",
            "16_export_scenario_inputs.R",
            "18_ma_level_uncertainty.R", "19_paired_bootstrap_contrasts.R",
            "20_verify_reported_numbers.R",
            "13_table_metadata.R"))  {
  message("\n>>> ", f)
  source(here::here("revision", "R", "analyse", f))
}
