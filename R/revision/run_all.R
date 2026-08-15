# R/revision/run_all.R ---------------------------------------------------------
# Runs the revision workflow in order. Each script is also runnable on its own,
# provided the earlier ones have been run (they exchange fitted objects through
# results/revision/intermediate/, which is gitignored and regenerable).
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
#   Rscript R/revision/06_validate_yang2024_reference.R
#
# Steps 07 and 08 must follow 05: 07 checks its baselines against the canonical
# summaries, and 08 reads them to draw the reference rules.
#
# Total runtime is a few minutes, dominated by the rma.mv fits in steps 01 and 03.
for (f in c("01_reproduce_original_analysis.R", "02_overshoot_diagnostics.R",
            "03_yang2024_bias_robust.R", "04_revision_sensitivity_summaries.R",
            "05_make_revision_tables.R", "07_influence_loo.R",
            "08_model_level_figure.R", "09_assumed_effect_scenarios.R",
            "11_main_figure.R", "12_supplementary_tables.R",
            "13_table_metadata.R"))  {
  message("\n>>> ", f)
  source(here::here("R", "revision", f))
}
