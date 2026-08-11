# Revision sensitivity analyses — technical note

Scope: the analyses added during the Biology Open revision that have been
independently verified. This is technical documentation, not manuscript prose.

**`S2_v2.R` is frozen.** Nothing in `R/revision/` edits, sources or runs it, and no
existing output is overwritten. The submitted analysis and its results remain exactly
as at commit `ac1e5cd` (`S2_v2.R`, md5 `a910b4cc1ac134b9792f2da5d0558ef9`).

**Yang et al. (2024) is being evaluated as a sensitivity analysis. It has not been
adopted as a replacement for the Yang-2023-style primary correction.** The 2024
paper's own framing supports that: *"Our proposed approach should be used as an
effective sensitivity analysis to understand the effects of publication bias on the
inferences drawn from a meta-analysis"* (Section 5.1, p.1605).

---

## 1. What the submitted analysis did

Per meta-analysis model (48 models, 28 source papers, 5,740 effect-size estimates,
in lnRR, SMD and Fisher's Zr):

1. Fit a multilevel model `rma.mv(yi = es, V = var, random = ~1|study_ID/obs_ID)`.
   The intercept is the uncorrected pooled mean, `beta0`.
2. Fit a bias-detection model adding a sampling-error moderator and a
   latest-year-centred publication year. Assign one of four scenarios from the
   **sign** of `beta0*beta1` and `beta0*beta2`.
3. Re-fit a scenario-specific correction model, once with a sampling-**error**
   moderator (`beta0_se_model`, Yang 2023 Eq. 2) and once with a sampling-**variance**
   moderator (`beta0_var_model`, Eq. 3).
4. Report `beta0_c3` = `beta0_var_model` when `|beta0| > |beta0_var_model|`, else
   `beta0`. A magnitude-only gate.
5. Compute power, Type M and Type S at two levels, and aggregate: weighted `lm` on
   the log scale over the 48 models, and `lmer` with a study random effect over the
   5,740 effect sizes.

Two properties of step 2 that the manuscript does not currently state, and which the
Methods text will need to reflect:

- **Scenario assignment uses sign only, with no significance test.** This is Yang
  2023's stated method, not an error in our port. Their paper defines the four
  scenarios purely by "expected sign" / "unexpected sign", and notes that the
  criterion captures a "(statistically non-significant) tendency", with the
  probability being 50% if there is no real effect. Our Methods currently describe a
  significance-based rule ("only the significant bias term"), which is a
  documentation error.
- Corrected means are deliberately paired with the **uncorrected** standard error
  when computing the metrics, holding precision fixed so only the assumed effect
  varies. This is Yang 2023's explicit choice.

## 2. Why these sensitivity analyses were added

Reviewer 1 asked whether the correction procedure itself drives the power, Type M and
Type S results. Reviewer 2 asked about self-reference in the assumed effect. Auditing
the correction step surfaced two mechanisms that make the corrected-side results
fragile:

- **Sign reversal.** The corrected mean can come out with the opposite sign to the
  uncorrected mean. The gate compares magnitudes only and is blind to a zero
  crossing, so small reversals pass it.
- **Divergence as the assumed effect approaches zero.** Type M tends to
  `2.3378 / |t|` as `t = |mu|/se` tends to 0, so a corrected mean near zero produces
  an arbitrarily large Type M. The relationship is inverted relative to intuition:
  the *smallest* reversals produce the *largest* Type M.

The Yang et al. (2024) estimators were added because they obtain a bias-robust point
estimate by **re-weighting** rather than by extrapolating an intercept to zero
sampling error. The relevant author claim is narrow and worth quoting accurately:

> "extrapolation involved in marginalizing predictor effects may yield poor estimates
> of publication-biased-adjusted effect size … **In contrast, the proposed two-step
> method does not rely on extrapolation** and relaxes the assumption of infinite
> precision" (Section 5.1, p.1605).

The paper does **not** claim that the estimator cannot overshoot or cannot cross
zero, and it does not use the words overshoot, over-correction or sign reversal
anywhere. Do not describe it as "never overshooting". What it delivers is the absence
of extrapolation; empirically the estimator still reverses sign in 6 of our 48 models.

## 3. Which script generates which table

| script | what it does |
|---|---|
| `R/revision/00_revision_functions.R` | shared functions: closed-form metrics, VCV construction, FE + VCV and UWLS, the CRVE sandwich, the two aggregation schemes. Reuses `load_datasets()` from `R/00_setup.R` so data handling cannot drift. |
| `01_reproduce_original_analysis.R` | reproduces the Yang-2023-style estimates (`beta0`, `beta0_se_model`, `beta0_var_model`, `beta0_c3`) as the fixed baseline |
| `02_overshoot_diagnostics.R` | sign reversals, shrinkage, the change in the assumed effect, `t` on both the fixed and the own standard error, and the meta-analysis-level metrics under `beta0_c3` |
| `03_yang2024_bias_robust.R` | FE + VCV and UWLS with CRVE; writes **`rho_sensitivity.csv`** |
| `04_revision_sensitivity_summaries.R` | aggregates power / Type M / Type S under each retained assumed effect, with provenance |
| `05_make_revision_tables.R` | writes **`per_meta_analysis_estimates.csv`**, **`reversal_counts.csv`**, **`primary_level_sensitivity.csv`**, **`meta_analysis_level_sensitivity.csv`** |
| `run_all.R` | runs 01–05 in order |

`results/revision/intermediate/` holds fitted objects passed between scripts. It is
gitignored and fully regenerable; delete it freely.

Type M is computed in **closed form** throughout, not by Monte Carlo. `S2_v2.R` draws
10,000 deviates per unit with no seed, which makes its Type M irreproducible across
runs and quantised to three decimals. The closed form agrees with the simulation to
0.34% and removes the dependence entirely.

### Provenance columns

Every aggregated row carries four columns that fix what the number means:

- `effect_estimator` — which assumed effect `mu` was used
- `se_source` — which standard error was paired with it (`se_beta0` or `own_CRVE`)
- `se_method` — how that standard error was computed (`CR1`; `NA` for a plain model
  standard error)
- `aggregation` — `meta_analysis_level` (48 units) or `primary_study_level` (5,740)
- `weighting` — `k_effect_sizes` or `unweighted`

**The finding that matters here is qualitative: pairing a bias-robust estimator with
its own cluster-robust standard error materially changes the meta-analysis-level
summary, because the robust standard error is larger than `se_beta0`.** Neither
figure below should be treated as a headline value yet.

The reproducible value in this repository uses the **CR1** sandwich implemented in
`00_revision_functions.R`: meta-analysis-level power is **0.708** with `se_beta0` and
**0.529** with the estimator's own CR1 standard error. Rows of the latter kind carry
`se_source = own_CRVE`, `se_method = CR1`, `verification_status = single_derivation`.

**Provisional alternative, not in the canonical tables.** An earlier implementation
used a **CR2 / Satterthwaite** correction (clubSandwich, and
`metafor::robust(..., clubSandwich = TRUE)`) and gave **0.479** for the same quantity.
That variant is **not implemented in this repository and has not been independently
reproduced**, so it is deliberately excluded from
`meta_analysis_level_sensitivity.csv`. It is recorded here only to document that the
CRVE variant itself moves the answer by about 0.05 on a quantity of order 0.5 — which
is a further reason no `own_CRVE` figure should be quoted without naming both the
convention and the method.

At the primary-study level the standard error is always the individual effect size's
own sampling error, so `se_source` does not vary there.

## 4. Yang et al. (2024) implementation, in full

So that another person can check we followed the method:

- **VCV construction** (paper Eq. 4), block-diagonal by clustering unit:
  diagonal `V[i,i] = v_i`; within-study off-diagonal `V[i,j] = rho * sqrt(v_i * v_j)`;
  zero between studies.
- **Clustering unit**: `study_ID` — the same grouping the submitted analysis uses as
  the outer random effect in `~1|study_ID/obs_ID`.
- **rho**: 0.5 by default, matching the authors' tutorial, with a sensitivity
  analysis over {0, 0.25, 0.5, 0.75} in `rho_sensitivity.csv`.
- **FE + VCV** (Eq. 3): fixed-effect GLS intercept, `beta = (1'V^-1 1)^-1 1'V^-1 y`,
  i.e. a weighted mean with weights `a = V^-1 1`. Implemented in closed form so the
  weights are inspectable and no optimizer is involved; agreement with a `metafor`
  fit was checked at 1.1e-04. **The weights are not constrained to be positive** —
  the paper says so directly (Section 2.2.1, p.1596: "the off-diagonal elements
  become negative values") — and `n_negative_weight` records this per meta-analysis:
  25 of 48 models have at least one negative weight, up to 48% of weights in one
  model. FE + VCV therefore carries no convex-combination guarantee.
- **UWLS**: the same estimator with a diagonal VCV, i.e. `rho = 0`, so `a = 1/v_i`.
  The paper notes UWLS, Henmi-Copas and IVhet are special cases of its framework with
  zero sampling correlation. Being a convex combination, UWLS is bounded by the
  observed effect sizes.
- **CRVE** (Eq. 5): sandwich estimator clustered by `study_ID`,
  `Var(beta) = (sum a)^-2 * sum_over_clusters ( sum_{i in cluster} a_i e_i )^2`, with
  `e = y - beta`. **CR0** unadjusted and **CR1** with the small-sample factor
  `J/(J-1)`; intervals use `t` on `J-1` degrees of freedom.
  **The CR2 / Satterthwaite correction is not implemented here.** An earlier pass used
  it and its standard errors differ from CR1 by up to 0.036 on standard errors of
  order 0.1–0.2. It has not been independently replicated, so all CRVE columns carry
  `verification_status = "single_derivation"`.
- **Which SE is used downstream**: both conventions are reported —
  `se_source = "se_beta0"` (the submitted analysis's convention, so the rows are
  comparable with the published numbers) and `se_source = "own_CRVE_SE_CR1"`.

The empirical property the method does deliver: **both estimators fell inside
`[min(observed effect), max(observed effect)]` in 48 of 48 models**. But the observed
effects straddle zero in **41 of 48**, so an estimate inside that range can still
fall either side of zero — which is why reversals remain. Bounded by the data is not
the same as unable to cross zero.

## 5. What is currently considered verified

`two_derivations` means the value was obtained twice by different routes — typically
a `metafor` fit and an independent closed-form implementation, or two separate
reviewers.

Verified:

- the Yang-2023-style estimates and the reversal counts 33 / 26 / 20 of 48;
- FE + VCV and UWLS point estimates, and their reversal counts 6 and 5 of 48;
- median `|estimate| / |beta0|` of 0.923 (FE + VCV) and 0.839 (UWLS);
- negative-weight counts; the 48/48 within-range result; 41/48 straddling zero;
- primary-study-level and meta-analysis-level summaries under `se_beta0` for all
  four assumed effects;
- the closed-form Type M against the Monte Carlo version, and against the
  `2.3378/|t|` limit;
- **external corpus**, documented here rather than tabulated because reproducing it
  needs another project's data: in the Yang et al. (2024) corpus, rebuilt from
  `Yefeng0920/WLS_RVE`'s raw `all effect size data edited 24-08.csv`, FE + VCV
  reverses sign in **66 of 456 (14.5%)** and UWLS in **50 of 456 (11.0%)**; on the
  paper's own 448-model analysis set, 65/447 and 48/447. So ~10–15% appears to be the
  method's baseline rate, and ours (12.5%) sits on it. Verification script:
  `/private/tmp/claude-501/meeverify/mee.R` (not committed; the raw data is 15 MB and
  belongs to another project).

Single-derivation, marked as such in the tables:

- all CRVE standard errors and the counts derived from them (CIs including zero:
  FE + VCV 19 of 48, UWLS 18 of 48). Under a CR2 variant the UWLS count was 20. The
  qualitative statement — roughly 40% of models have a bias-robust estimate
  indistinguishable from zero — holds at 17–20 of 48 under every variant tried;
- the meta-analysis-level rows with `se_source = "own_CRVE_SE_CR1"`, since they
  inherit the CRVE standard error.

Not carried into these tables at all: reversal counts for the Yang et al. **2023**
87-model corpus, and meta-analysis-level Type M maxima from the earlier audit passes.
Both are single-sourced and are recorded only in `docs/09_overshoot_investigation.md`.

## 6. Methodological decisions still open

These are recorded in `docs/07_final_audit_report.md` and are not resolved by
anything in this directory:

1. **The magnitude-only gate.** Whether to retain it, floor `|mu|`, restrict which
   Type M values are interpreted, or require sign preservation. Adopting the 2024
   estimators reduces sign reversals from 20 to 6 of 48 but does **not** resolve the
   `mu -> 0` divergence: the maximum primary-level Type M is still 1,029 under
   FE + VCV and 609 under UWLS, against 26,326 as submitted. A reporting rule for
   near-zero assumed effects is needed regardless of which estimator is chosen.
2. **The dependence structure and estimand** at the primary-study level: whether to
   weight study clusters or meta-analyses equally, and whether the 48 models are the
   population or a sample.
3. **Whether to build an optimistic-effect scenario on the variance-model estimate**,
   given that 37 of 48 of those estimates have a confidence interval including zero.
4. **Type S reporting.** The log-scale model needs an offset because Type S can be
   exactly zero, and at the meta-analysis level the offset (0.025) exceeds the
   estimate itself. The lower confidence limit is negative under the uncorrected
   estimate, under a sign-preserving gate and under FE + VCV alike, so this is not a
   consequence of the correction choice.
5. **Whether Yang et al. (2024) becomes a reported sensitivity analysis, stays a
   response-to-reviewers calculation, or is not used.** Note it is not implemented in
   `orchaRd` (version 2.2.1 inspected): `pub_bias_plot()` takes already-fitted models
   and only adds plot layers, and no exported function implements UWLS, FE + VCV or
   CRVE. Step 1 and Step 2 must be fitted directly, as done here and in the authors'
   tutorial.
