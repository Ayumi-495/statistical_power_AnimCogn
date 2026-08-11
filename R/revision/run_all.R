# R/revision/run_all.R ---------------------------------------------------------
# Runs the revision workflow in order. Each script is also runnable on its own,
# provided the earlier ones have been run (they exchange fitted objects through
# results/revision/intermediate/, which is gitignored and regenerable).
#
# Total runtime is a few minutes, dominated by the 144 rma.mv fits in step 01.
for (f in c("01_reproduce_original_analysis.R", "02_overshoot_diagnostics.R",
            "03_yang2024_bias_robust.R", "04_revision_sensitivity_summaries.R",
            "05_make_revision_tables.R")) {
  message("\n>>> ", f)
  source(here::here("R", "revision", f))
}
