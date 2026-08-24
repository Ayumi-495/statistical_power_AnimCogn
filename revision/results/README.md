# Revision analysis outputs

This directory contains the machine-readable outputs of the revised analysis. The
analysis estimates statistical power, Type M error and Type S error at the
primary-study and meta-analysis levels for 48 meta-analytic models from 28 source
papers, comprising 5,740 effect-size estimates.

The current workflow is implemented in `revision/R/`. The analysis submitted before
revision is retained separately as `S2_v2.R` in the repository root.

## Output index

| File | Contents |
|---|---|
| `per_meta_analysis_estimates.csv` | Uncorrected, Yang-2023-style bias-corrected and Yang-2024 bias-robust estimates for each meta-analytic model, with standard errors and provenance fields. |
| `model_level_metrics.csv` | Power, Type M error and Type S error for each model under each assumed-effect specification. |
| `primary_level_sensitivity.csv` | Primary-study-level summaries for the retained assumed-effect specifications and aggregation methods. |
| `meta_analysis_level_sensitivity.csv` | Meta-analysis-level summaries using effect-size-count and equal-model weighting. |
| `assumed_effect_scenarios.csv` | Optimistic and externally specified assumed-effect scenarios at both analysis levels. |
| `leave_one_cluster_out.csv` | Primary-study-level sensitivity to recomputing the assumed effect after removing each study cluster. |
| `loo_influence.csv` | Leave-one-meta-analysis-model-out influence results. |
| `leave_one_paper_out.csv` | Leave-one-source-paper-out influence results. |
| `ma_level_uncertainty.csv` | Alternative uncertainty estimates for the meta-analysis-level summaries. |
| `ma_level_paired_contrasts.csv` | Paired paper-level bootstrap contrasts between assumed-effect specifications. |
| `reversal_counts.csv` | Counts of sign reversals under the bias-correction approaches. |
| `rho_sensitivity.csv` | Sensitivity of the Yang-2024 bias-robust estimates to the assumed within-study sampling correlation. |
| `evidence_base_characteristics.csv` | Characteristics of the 28 source papers and their contributions to the analysed evidence base. |
| `verification_audit.csv` | Independent checks of metric calculations, estimands and summary statistics. |
| `yang2024_reference_validation.csv` | Validation against the published Yang et al. worked example. |

`supplementary/metadata_files.csv` provides a file-level index with row counts.

## Supplementary tables

The header-first supplementary files are:

- `supplementary/TableS1_evidence_base.csv`: characteristics of the 28 source papers;
- `supplementary/TableS2_reported_metrics.csv`: reported power, Type M error and Type S error results and assumed-effect sensitivities;
- `supplementary/TableS3_leave_one_model_out.csv`: leave-one-model-out influence results; and
- `supplementary/TableS4_leave_one_paper_out.csv`: leave-one-paper-out influence results.

Captions are in `supplementary/captions.md`. Column definitions and the complete file
index are in `supplementary/metadata_columns.csv` and
`supplementary/metadata_files.csv`, respectively. CSV exports preserving the layout of
the final Google Sheet are available in the repository-level `supplementary_tables/`
directory.

## Reproduction

Restore the package versions recorded in `renv.lock`, then run the workflow from the
repository root:

```bash
Rscript revision/R/reproduce/01_estimates.R
Rscript revision/R/reproduce/03_baseline.R
Rscript revision/R/reproduce/05_loo_study.R
Rscript revision/R/analyse/run_all.R
```

The leave-one-study-out step takes approximately 30 minutes. The following independent
checks are run separately:

```bash
Rscript revision/R/analyse/06_validate_yang2024_reference.R
python3 revision/R/analyse/17_verify_scenarios.py
python3 revision/R/analyse/21_audit_manuscript_claims.py
```

Script 06 downloads a checksum-pinned worked-example dataset and therefore requires
network access. Script 10 reads evidence-base metadata from a sibling project through
a local path; its output, `evidence_base_characteristics.csv`, is included here so the
remaining workflow can run offline.

## Analysis levels and provenance

At the primary-study level, each effect-size estimate is evaluated using its sampling
standard error and the corresponding meta-analytic assumed effect. At the
meta-analysis level, each of the 48 models is evaluated using the uncertainty of its
pooled estimate. Results are summarised both by the number of effect-size estimates
contributed by each model and with equal weight per model where indicated.

The output tables retain fields identifying the assumed-effect estimator, standard
error source, critical-value convention, analysis level, weighting and analytical role.
Power and Type S error are bounded probabilities. Type M error is a ratio conditional
on statistical significance and can become large when the assumed effect approaches
zero.

## Verification

The archived verification records confirm the expected hierarchy of 48 models, 28
source papers and 5,740 effect-size estimates. `verification_audit.csv` records the
independent R checks, and `yang2024_reference_validation.csv` records validation against
the external worked example. The Python scenario verification and manuscript-number
audit provide independent checks of the exported results.

## Licences

The analysis code is licensed under the MIT License. Project-created metadata,
supplementary-table exports and derived result tables are licensed under CC BY 4.0.
See `LICENSE-CODE` and `LICENSE-DATA` in the repository root for scope and third-party
data exclusions.
