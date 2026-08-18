# Statistical fragility in the meta-analytic evidence base of animal cognition

Data and code for Mizuno et al., *Biology Open* (bio.062841): a second-order
meta-analysis of statistical power, Type M (magnitude) error and Type S (sign) error
across **48 meta-analytic models from 28 published meta-analytical papers**, covering
**5,740 effect-size estimates** in log response ratio, standardised mean difference and
Fisher's *Zr*.

## Layout

| | |
|---|---|
| `SMD/`, `lnRR/`, `zr/` | the effect-size datasets, one CSV per meta-analytic model |
| `old/` | superseded versions of two datasets, kept as a record of what data cleaning changed: `old/MA20/` was merged into `lnRR/MA20_cleaned.csv`, and `old/MA41_0*.csv` are earlier versions of `zr/MA41_0*.csv` with one extra row and a different column set. Not read by any script. |
| `S2_v2.R` | **the analysis as submitted.** Frozen: nothing in this repository edits, sources or runs it. |
| `scatter_plot.R`, `figure/` | the figures as submitted |
| `revision/` | **everything added during revision.** See below. |

## `revision/`

```
revision/
  R/
    reproduce/   re-implements S2_v2.R from scratch, to check that it reproduces and to
                 find defects in it. Also produces the model fits the revision analyses
                 build on. Run: 01 -> 03 -> 05.
    analyse/     the revision analyses. Run: run_all.R
  results/
    reproduce/   what R/reproduce/ produces
    *.csv        the revision results
    supplementary/  Supplementary Tables S1 and S2, their captions and column metadata
    figures/     the revision figures
    README.md    the technical note: what was done, what was decided, and why
```

Start with `revision/results/README.md`. It documents every methodological decision, the
two defects found in the submitted analysis, and what is verified against what.

## Reproducing

```bash
Rscript revision/R/reproduce/01_estimates.R
Rscript revision/R/reproduce/03_baseline.R
Rscript revision/R/reproduce/05_loo_study.R      # ~30 min
Rscript revision/R/analyse/run_all.R
```

Then the checks that are not part of the standard run, because they need the network or
Python:

```bash
Rscript revision/R/analyse/06_validate_yang2024_reference.R
python3  revision/R/analyse/17_verify_scenarios.py
python3  revision/R/analyse/21_audit_manuscript_claims.py
```

`run_all.R` is deterministic: a rerun reproduces every CSV byte for byte.

Requires R with `metafor`, `lme4`, `nlme`, `clubSandwich`, `dplyr`, `ggplot2`, `here`,
and Python 3 (standard library only).

## Verification

Every number reported in the paper has two derivations by different routes. The current
state is recorded in `revision/results/verification_audit.csv` and summarised in
`revision/results/README.md`.

| check | result |
|---|---|
| `20_verify_reported_numbers.R` — the bias-correction gate, the three metric definitions against Monte Carlo, the model-level estimands (`lme4` against `nlme`), the meta-analysis-level summaries | 12 of 12 pass |
| `17_verify_scenarios.py` — every assumed-effect scenario row re-derived from scratch in Python | 123 of 123, within 1e-6 |
| `21_audit_manuscript_claims.py` — every number written into the manuscript, against these files | 90 of 90 |
| `06_validate_yang2024_reference.R` — the CR2 implementation against the method authors' published worked example | 9 of 9 published values |
