This repository contains the data-processing and analysis code for the manuscript:
**The statistical fragility of animal cognition findings: a meta-meta-analytic reappraisal**

The study evaluates estimated statistical power, Type M magnitude error, and Type S sign error in a defined subset of published animal-cognition meta-analyses. The analysis includes 28 meta-analytical papers, represented by 48 meta-analytic models and 5,740 primary studies. The estimates apply to this meta-analytic evidence base and should not be interpreted as field-wide estimates for all animal cognition research.

## Repository structure

- `S2_v2.R`: Main analysis script. It imports the effect-size datasets, fits the meta-analytic and publication-bias models, calculates estimated power and Type M and Type S errors, aggregates results across meta-analytic models, and produces the main summaries and figures.
- `lnRR/`: Input datasets using log response ratios.
- `SMD/`: Input datasets using standardised mean differences.
- `SMD/des_stat/`: Datasets supplied as descriptive statistics, from which effect sizes and sampling variances are calculated in the script.
- `Zr/`: Input datasets using Fisher's z-transformed correlations.

The script uses paths relative to the repository root. Run it with the repository root as the working directory.

## Software requirements

The analyses reported in the manuscript were conducted in R 4.4.2. The main script uses the following R packages:

`pacman`, `tidyverse`, `dplyr`, `tidyr`, `purrr`, `stringr`, `readr`, `janitor`, `here`, `ggplot2`, `gghighlight`, `ggrepel`, `patchwork`, `metafor`, `esc`, `orchaRd`, `lme4`, `kableExtra`, and `gt`.

Install any missing packages before running the analysis.

## Reproducing the analysis

1. Clone or download this repository.
2. Open R or RStudio and set the working directory to the repository root.
3. Confirm that the `lnRR/`, `SMD/`, and `Zr/` directories are present.
4. Run the main script:

```r
source("S2_v2.R")
```

The script reads the input files, reconstructs effect sizes where required, fits uncorrected and publication-bias-adjusted models, and calculates the statistical reliability metrics used in the manuscript.

## Interpretation of the main quantities

- **Estimated power** is calculated retrospectively using the meta-analytic mean, bias-corrected when indicated, as the assumed underlying effect. It is not a prospective power calculation for a planned study.
- **Type M error** is the expected exaggeration in the absolute effect-size estimate, conditional on statistical significance.
- **Type S error** is the probability that a statistically significant estimate has the opposite sign from the assumed underlying effect.
- The assumed underlying or “true” effect is not directly observed. It is represented by the uncorrected or bias-corrected meta-analytic mean, so all reported metrics are conditional on this proxy and the model assumptions.

## Data provenance

The effect sizes and sampling variances were extracted or reconstructed from published articles, supplementary materials, public repositories, and data shared by study authors. Details of study identification, inclusion criteria, and reconstruction procedures are provided in the manuscript and its supplementary material.

## Preregistration

The study was preregistered on OSF: [https://osf.io/8hnkb](https://osf.io/8hnkb)
Also, here is the preprint: [https://doi.org/10.32942/X2Z36M](https://doi.org/10.32942/X2Z36M)

## Contact
Ayumi Mizuno  
University of Alberta  
ayumi.mizuno5[at]gmail.com or amizuno[at]alberta.ca
