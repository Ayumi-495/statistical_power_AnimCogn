# Statistical fragility in the meta-analytic evidence base of animal cognition

Data, code and metadata for Mizuno et al., *Biology Open*. The study estimates statistical power, Type M (magnitude) error and Type S (sign) error in a defined evidence base of 48 meta-analytic models from 28 published meta-analytical papers, comprising 5,740 effect-size estimates. These results describe the included evidence base and are not field-wide estimates for all animal cognition research.

## Which analysis to use

The current analysis is under `revision/R/`. Run the scripts in `revision/R/reproduce/` first, followed by `revision/R/analyse/run_all.R`.

`S2_v2.R` is the analysis submitted before revision. It is retained unchanged as an audit record. The current workflow does not edit, source or run it.

## Repository contents

| Path | Contents |
|---|---|
| `SMD/`, `lnRR/`, `zr/` | Current effect-size input datasets, one CSV per meta-analytic model. |
| `old/` | Superseded versions of two datasets. No current script reads this directory. |
| `S2_v2.R` | Frozen analysis submitted before revision. |
| `scatter_plot.R`, `figure/` | Figure code and outputs submitted before revision. |
| `revision/R/reproduce/` | Independent reimplementation of the submitted analysis and the fitted objects used by the revised workflow. |
| `revision/R/analyse/` | Current revision analyses and verification scripts. |
| `revision/results/` | Machine-readable results and technical documentation. |
| `supplementary_tables/` | Portable CSV exports of the five tabs in the final Google Sheet: metadata and Tables S1-S4. |
| `renv.lock`, `session-info.txt` | R package lockfile and the environment used to prepare this repository snapshot. |

An output index and additional reproduction notes are in `revision/results/README.md`.

## Reproducing the current analysis

Use R 4.6.0. From the repository root, restore the recorded package environment:

```r
install.packages("renv")
renv::restore()
```

Run the standard workflow:

```bash
Rscript revision/R/reproduce/01_estimates.R
Rscript revision/R/reproduce/03_baseline.R
Rscript revision/R/reproduce/05_loo_study.R
Rscript revision/R/analyse/run_all.R
```

The leave-one-study-out step takes about 30 minutes. `run_all.R` takes a few minutes
and recreates the CSV outputs and revision figures.

### Additional checks and external inputs

```bash
Rscript revision/R/analyse/06_validate_yang2024_reference.R
python3 revision/R/analyse/17_verify_scenarios.py
python3 revision/R/analyse/21_audit_manuscript_claims.py
```

Script 06 downloads a pinned worked-example dataset and checks its MD5 checksum, so it requires network access. Script 10 reads the sibling `systematic_mapping_AnimCogn` project through a local absolute path and is not portable. Its committed output, `revision/results/evidence_base_characteristics.csv`, is an input to the standard offline workflow. The Python scripts use the standard library only.

## Supplementary tables

`supplementary_tables/` contains CSV exports of the exact tabs in the final Google Sheet. The four table files retain the caption in the first row, a blank separator row and the table header in the third row. Header-first, script-generated versions and their column definitions remain under `revision/results/supplementary/`.

Both the exported tabs and the script-generated files use Table S1 for evidence-base characteristics and Table S2 for the reported metrics.

## Verification

The verification records distributed with this repository are in `revision/results/verification_audit.csv` and `revision/results/README.md`. 
They cover:

- the bias-correction gate and metric definitions;
- independent `lme4` and `nlme` derivations of the model-level estimands;
- an independent Python derivation of all assumed-effect scenarios;
- checks of the numbers reported in the manuscript; and
- external validation against the Yang et al. worked example.

The expected hierarchy is asserted in code: 48 meta-analytic models, 28 source papers and 5,740 effect-size rows.

For this pre-deposit preparation, the lockfile was restored into an empty R library. The baseline reproduction passed 13/13 checks, the independent scenario verification passed 123/123 checks, the manuscript-claim audit passed 93/93 checks, and the pinned Yang et al. worked example passed 9/9 checks after its MD5 checksum was confirmed. The approximately 30-minute leave-one-study-out step was not rerun because its code was unchanged; its existing output is retained.

## Use of AI-assisted tools

Claude Code (Anthropic) and OpenAI Codex were used during revision to assist with code audit, implementation, testing and technical documentation. The authors made the methodological decisions, reviewed the code and outputs, and remain responsible for the analysis and manuscript. The AI tools were not treated as authors or as evidence sources.

## Data provenance

Effect sizes and sampling variances were extracted or reconstructed from published articles, supplementary materials, public repositories and data shared by study authors. Study identification, inclusion criteria and reconstruction procedures are described in the manuscript and supplementary material.

## Licences

Analysis code is available under the MIT License (`LICENSE-CODE`). Metadata, supplementary-table exports and derived result tables created for this project are available under CC BY 4.0 (`LICENSE-DATA`). The CC BY 4.0 licence does not extend to third-party input data under `SMD/`, `lnRR/`, `zr/` or `old/`; consult the original sources before reusing those files.

## Related records

- Preregistration: <https://osf.io/8hnkb>
- Preprint: <https://doi.org/10.32942/X2Z36M>

## Contact

Ayumi Mizuno, University of Alberta<br>
ayumi.mizuno5[at]gmail.com or amizuno[at]ualberta.ca
