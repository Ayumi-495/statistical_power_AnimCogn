# Revision sensitivity analyses — technical note

Scope: the analyses added during the Biology Open revision that have been
independently verified. This is technical documentation, not manuscript prose.

**`S2_v2.R` is frozen.** Nothing in `revision/R/analyse/` edits, sources or runs it, and no
existing output is overwritten. The submitted analysis and its results remain exactly
as at commit `ac1e5cd` (`S2_v2.R`, md5 `a910b4cc1ac134b9792f2da5d0558ef9`).

**Yang et al. (2024) is a sensitivity analysis, not a replacement for the
Yang-2023-style primary correction.** The 2024 paper's own framing supports that:
*"Our proposed approach should be used as an effective sensitivity analysis to
understand the effects of publication bias on the inferences drawn from a
meta-analysis"* (Section 5.1, p.1605).

## Decisions this directory now implements

Settled, and encoded in the `role` column of the summary tables so no reader has to
infer them:

1. **Yang et al. (2023) remains the primary analysis** — `role = "primary"`. It keeps
   the paper consistent with its framework of origin and directly comparable to it.
2. **FE + VCV with its own cluster-robust standard error is the reported sensitivity
   analysis** — `role = "reported_sensitivity"`.
3. **At the meta-analysis level the FE + VCV estimate is paired with its own
   cluster-robust standard error, not with the original `se_beta0`.** They are outputs
   of the same fitted analysis; combining the new estimate with the old standard error
   would create a hybrid quantity that neither model estimates. The hybrid is retained
   only as `role = "diagnostic_hybrid_not_reported"`.
4. **Type M and Type S are handled by interpretation, not by arithmetic.** No floor on
   `|mu|`, no exclusion rule and no sign-preserving correction is applied anywhere in
   this directory. Very large Type M values are expected when the assumed true effect
   is close to zero, because the assumed effect is the denominator, and they are not
   evidence of an effect in the opposite direction — particularly since the confidence
   intervals of the reversed estimates include zero.
5. **UWLS is supplementary** — `role = "supplementary"` — not co-primary with FE + VCV.
6. **The critical value stays `z = 1.96` for power, Type M and Type S, in both the
   Yang-2023 and the FE + VCV analyses.** This is Yang et al.'s own choice — their
   published script defines `power.ma_Shinichi` with the identical `qnorm(1 - alpha/2)`
   (`EcoEvo_PB_script.Rmd:46-48`), and `error_S` and `error_M` likewise — and it gives
   both approaches the same decision rule, which is what makes them comparable. The
   Methods must state it explicitly: these design metrics use a common normal-theory
   two-sided 5% threshold rather than each fitted model's own inferential reference
   distribution. The Satterthwaite alternative remains a labelled diagnostic only.
7. **The `k`-weighted meta-analysis-level aggregate is a secondary descriptive
   summary, not the headline.** The paper centres on the primary-study-level findings.
   `MA09` is **not** excluded from anything; its influence is reported as a diagnostic
   (`07_influence_loo.R`), alongside an equally weighted model-level summary and a
   figure of all 48 model values (`08_model_level_figure.R`).

Open items are in section 6.

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
| `revision/R/analyse/00_revision_functions.R` | shared functions: closed-form metrics, VCV construction, FE + VCV and UWLS, the CRVE sandwich, the two aggregation schemes. Reuses `load_datasets()` from `revision/R/reproduce/00_setup.R` so data handling cannot drift. |
| `01_reproduce_original_analysis.R` | reproduces the Yang-2023-style estimates (`beta0`, `beta0_se_model`, `beta0_var_model`, `beta0_c3`) as the fixed baseline |
| `02_overshoot_diagnostics.R` | sign reversals, shrinkage, the change in the assumed effect, `t` on both the fixed and the own standard error, and the meta-analysis-level metrics under `beta0_c3` |
| `03_yang2024_bias_robust.R` | FE + VCV and UWLS with CRVE; writes **`rho_sensitivity.csv`** |
| `04_revision_sensitivity_summaries.R` | aggregates power / Type M / Type S under each retained assumed effect, with provenance |
| `05_make_revision_tables.R` | writes **`per_meta_analysis_estimates.csv`**, **`reversal_counts.csv`**, **`primary_level_sensitivity.csv`**, **`meta_analysis_level_sensitivity.csv`** |
| `06_validate_yang2024_reference.R` | validates the CR2 implementation against the authors' published worked example; writes **`yang2024_reference_validation.csv`**. **Not in `run_all.R`** — it downloads a data file from another project's repository, and `run_all.R` must work offline. Re-run it after any change to the CR2 code path. |
| `07_influence_loo.R` | leave-one-model-out influence on the meta-analysis-level summaries, all 48 models × 4 specifications × 3 metrics × 2 weightings; writes **`loo_influence.csv`** |
| `08_model_level_figure.R` | the per-model values behind those summaries; writes **`model_level_metrics.csv`** and **`figures/model_level_metrics.{pdf,png}`** |
| `09_assumed_effect_scenarios.R` | optimistic and externally specified assumed effects at both levels; writes **`assumed_effect_scenarios.csv`** |
| `10_evidence_base_table.R` | the 28-paper characteristics table; writes **`evidence_base_characteristics.csv`**. **Not in `run_all.R`** — it reads the sibling systematic map. |
| `11_main_figure.R` | the replacement for manuscript Figure 3; writes **`figures/main_metrics.{pdf,png}`** |
| `12_supplementary_tables.R` | writes **`supplementary/TableS1_reported_metrics.csv`** and **`supplementary/TableS2_evidence_base.csv`** |
| `13_table_metadata.R` | captions, the column dictionary and the file index; writes **`supplementary/{captions.md, metadata_columns.csv, metadata_files.csv}`**. Runs last, because its gates check every other output. |
| `14_leave_one_cluster_out.R` | leave-one-cluster-out at the primary-study level; writes **`leave_one_cluster_out.csv`** |
| `15_leave_one_paper_out.R` | leave-one-source-paper-out at the meta-analysis level; writes **`leave_one_paper_out.csv`** |
| `16_export_scenario_inputs.R` | exports the fitted inputs that `17` needs |
| `17_verify_scenarios.py` | **verification.** Independent Python re-derivation of every row of `assumed_effect_scenarios.csv`. **Not in `run_all.R`** (Python). |
| `18_ma_level_uncertainty.R` | three ways of estimating the meta-analysis-level interval; writes **`ma_level_uncertainty.csv`** |
| `19_paired_bootstrap_contrasts.R` | paired paper-level bootstrap contrasts; writes **`ma_level_paired_contrasts.csv`**. **Computed but not currently used** — see the note on that file below. |
| `20_verify_reported_numbers.R` | **verification.** Second derivations for the gate, the metric definitions, the two model-level estimands and the meta-analysis-level summaries; writes **`verification_audit.csv`** |
| `21_audit_manuscript_claims.py` | **verification.** Checks every number written into the manuscript text against these files. **Not in `run_all.R`** (Python). |
| `run_all.R` | runs 01–05, 07, 08, 09, 11, 12, 14, 15, 16, 18, 19, then 13 last |

Scripts 06, 10, 17 and 21 are run separately: 06 needs the network, 10 needs the sibling
systematic map, and 17 and 21 are Python. Full verification is therefore:

```
Rscript revision/R/analyse/run_all.R
Rscript revision/R/analyse/06_validate_yang2024_reference.R
Rscript revision/R/analyse/20_verify_reported_numbers.R
python3  revision/R/analyse/17_verify_scenarios.py
python3  revision/R/analyse/21_audit_manuscript_claims.py
```

`clubSandwich` is a hard requirement. `00_revision_functions.R` stops if it is absent
rather than falling back to the hand-written CR1 sandwich, which is a different
specification. Versions used: R 4.6.0, `metafor` 5.0.1, `clubSandwich` 0.7.0, `lme4`
2.0.6; the first three are also stamped into `per_meta_analysis_estimates.csv`.

`revision/results/intermediate/` holds fitted objects passed between scripts. It is
gitignored and fully regenerable; delete it freely.

Type M is computed in **closed form** throughout, not by Monte Carlo. `S2_v2.R` draws
10,000 deviates per unit with no seed, which makes its Type M irreproducible across
runs and quantised to three decimals. The closed form agrees with the simulation to
0.34% and removes the dependence entirely.

### Provenance columns

Every aggregated row carries the columns that fix what the number means:

- `role` — what the row is for: `primary`, `reported_sensitivity`, `supplementary`,
  `reference_uncorrected`, or `diagnostic_*` (**not** for reporting; retained so that a
  choice is visible rather than assumed)
- `effect_estimator` — which assumed effect `mu` was used
- `se_source` — which standard error was paired with it (`se_beta0`, `own_CRVE`, or
  `own_sampling_error_per_effect_size`)
- `se_method` — how that standard error was computed (`CR2_Satterthwaite`,
  `CR2_naive_t`, `CR1`; `NA` for a plain model standard error)
- `crit_value_method` — which critical value defines the metric
- `aggregation` — `meta_analysis_level` (48 units) or `primary_study_level` (5,740)
- `weighting` — `k_effect_sizes` or `unweighted`

At the primary-study level the standard error is always the individual effect size's
own sampling error, so `se_source` does not vary there.

### The critical value is an explicit choice

Adopting the 2024 estimator introduces a second candidate critical value: the
Satterthwaite `t` that clubSandwich reports for the CR2 test. **The canonical choice
throughout is `z = 1.96`**, the design-analysis convention of Gelman & Carlin (2014)
and the one the submitted analysis uses, so that a sensitivity analysis varies the
assumed effect and its standard error and *not* the test that defines the metric.

This is not cosmetic. The Satterthwaite degrees of freedom here have median 10.4 and
range 1.8–170.6, so `qt(0.975, df)` often exceeds 2.4 (11 of 48 models have df < 5).
Substituting it lowers meta-analysis-level power from 0.479 to **0.306** and raises
Type M from 1.53 to 1.65. That row is in the table as
`role = "diagnostic_critical_value"` and must not be quoted as a result.

**The submitted analysis already declines to use its fitting model's own reference
distribution.** All 48 intercept-only models are `rma.mv(..., test = "t")`, and their
degrees of freedom are exactly `k − 1` (verified for all 48; range 3–1296), giving
`qt(0.975, df)` a median of 2.028 and a maximum of 3.182. The metric uses 1.96 anyway.
So "use the model's own test" is not an obligation created by the 2024 estimator — it
is one the original metric already declines, and applying it consistently would move
every row:

| specification | z = 1.96 | each model's own t |
|---|---|---|
| uncorrected | 0.82207 | 0.81138 |
| Yang-2023 gated `beta0_c3` | 0.39038 | 0.37959 |
| FE + VCV, own CR2 SE | 0.47895 | 0.30578 |

Two derivations: R, and an independent re-implementation in pure Python (normal CDF
from `erf`, t quantiles from a hand-coded incomplete beta plus bisection, weighted
geometric mean without `lm()`); all values agree to five decimals.

**Using z for the Yang-2023 rows and Satterthwaite for the FE + VCV row is the one
combination that is not defensible**, and it is what a casual reading of "pair the
estimator with its own uncertainty" produces. A further reason to prefer z: the
Satterthwaite df is a *dependence* quantity — it tracks the effective number of
independent study clusters under the FE + VCV weighting (Spearman 0.9997 across the 48
models, sitting about 9% below it), whereas `k − 1` ignores clustering entirely. The two
df notions are not commensurable, so switching one row alone would make that row carry
a pseudo-replication penalty no other row carries.

**This choice is not settled** and is open item 6 below.

### Meta-analysis-level power under each pairing

| pairing | `role` | power | Type M | Type S |
|---|---|---|---|---|
| uncorrected `beta0`, `se_beta0` | `reference_uncorrected` | 0.822 | 1.109 | 0.00058 |
| **submitted `beta0_c3`, `se_beta0`** | **`primary`** | **0.390** | **2.177** | **0.0122** |
| **FE + VCV, own CR2 SE** | **`reported_sensitivity`** | **0.479** | **1.526** | **0.0056** |
| UWLS, own CR2 SE (naive-t) | `supplementary` | 0.383 | 2.011 | 0.0230 |
| FE + VCV, own CR1 SE | `diagnostic_crve_variant` | 0.529 | 1.433 | 0.0031 |
| FE + VCV, hybrid `se_beta0` | `diagnostic_hybrid_not_reported` | 0.708 | 1.239 | 0.0022 |
| FE + VCV, own CR2 SE, Satterthwaite crit. | `diagnostic_critical_value` | 0.306 | 1.646 | 0.0030 |

**The reported figure is 0.479 [0.368, 0.624].** Note it lies *between* the submitted
0.390 and the uncorrected 0.822: the sensitivity analysis raises power relative to the
submitted analysis while leaving it far below the uncorrected value.

**Why CR2 gives 0.479 where CR1 gave 0.529.** The unweighted median CR2/CR1 standard
error ratio is only 1.008 (range 0.968–1.344; CR2 is larger in 42 of 48). But the
aggregation is `k`-weighted, and the `k`-weighted mean log ratio is **+7.65%**, driven
almost entirely by `MA09` — 1,297 of the 5,740 rows, ratio 1.317, Satterthwaite df 4.2.
Its power falls 0.200 → 0.135. The `k`-weighted mean log power ratio is −9.39%, which
reproduces 0.529 → 0.479 exactly. So this difference is one large dataset, not a
general property of the correction, and it is a concrete instance of the open question
in section 6 about how much weight `k`-weighting places on individual models.

### How much of the meta-analysis-level summary is one dataset

Separate question from the CRVE and critical-value choices, and larger than either.
`MA09` holds **1,297 of 5,740 effect sizes (22.6%)** and therefore 22.6% of the
`k`-weight. Dropping it:

| specification | all 48 | dropping MA09 | change |
|---|---|---|---|
| uncorrected, z | 0.82207 | 0.77639 | −5.6% |
| Yang-2023 gated `beta0_c3`, z | 0.39038 | 0.38060 | −2.5% |
| **FE + VCV, own CR2 SE, z** | **0.47895** | **0.69313** | **+44.7%** |
| FE + VCV, own CR2 SE, Satterthwaite | 0.30578 | 0.59967 | +96.1% |

So the Yang-2023 summary is barely sensitive to `MA09` and the FE + VCV summary is
dominated by it. **This is not an argument for excluding `MA09`.** It is a property of
`lm(log(metric) ~ 1, weights = k)` over 48 units with very unequal `k`: the top 1, 5 and
10 models hold 22.6%, 62.3% and 75.0% of all effect sizes. Weighting each
meta-analysis equally instead changes the Yang-2023 summary by −36% and the FE + VCV
summary by +11%, so **the weighting convention moves this quantity more than the choice
of estimator does.** For the submitted specification the most influential single model
is `MA26` (+0.076), not `MA09`.

Why `MA09` behaves as it does under FE + VCV: its estimate is **−0.127** against
`beta0` = **+0.219** — one of the six reversals — with a CR2 interval of
[−0.538, +0.283]. Its CR2 standard error is **0.151** against a working standard error
of **0.0025**, a 60-fold inflation, because **4.6 of its 126 nominal study clusters
carry the weight** (one holds 44%, the top three 62%). Its FE + VCV power is 0.135. It
is also the largest dataset, the only one affected by the `~` defect (section 4b), and
the one whose `study_ID` values are opaque codes `CD001`–`CD126`, so the clustering
cannot easily be checked against the source papers.

Two derivations (R and the independent Python re-implementation described above).

### Full leave-one-out, and the equally weighted summary

Reporting the MA09 figure alone would invite the question of why that model was
singled out, and the honest answer is that it is the largest. `loo_influence.csv` gives
the **full leave-one-out** — all 48 models × 4 specifications × 3 metrics × 2 weightings,
1,152 rows — so MA09 appears as the largest of 48 influence values rather than as a
special case. Baselines are checked against `meta_analysis_level_sensitivity.csv` over
all 24 cells to 8.9e-16 (a weighted mean of logs against the canonical `lm()` route, so a
second derivation).

For power, the three largest influences on each summary:

| | weighting | most influential | second | third | median \|change\| | n > 10% | n > 20% |
|---|---|---|---|---|---|---|---|
| Yang-2023 (primary) | by effect-size count | MA26 **+19.5%** | MA31 −12.9% | MA08 −6.9% | 0.6% | 2 | 0 |
| Yang-2023 (primary) | equal per model | MA39_1 +3.5% | MA26 +3.5% | MA13_03 +3.4% | 2.1% | 0 | 0 |
| FE + VCV, own CR2 SE | by effect-size count | MA09 **+44.7%** | MA31 −10.2% | MA08 −7.2% | 0.4% | 2 | 1 |
| FE + VCV, own CR2 SE | equal per model | MA22_02 +5.1% | MA39_2 +4.5% | MA39_1 +3.6% | 1.3% | 0 | 0 |

Two things follow. Single-model leverage is not unique to the FE + VCV summary — the
Yang-2023 summary has its own influential model in MA26 — and **MA09 is the only model
anywhere in the table that moves a summary by more than 20%**.

But that leverage is **largely a property of the weighting, not of any model**. For
power, under equal weighting no model moves any summary by more than **5.1%**, against
44.7% under weighting by effect-size count, and the most influential model is a different
one in every specification.

The same contrast holds for the other two metrics but is less complete, and the table
below is the honest version rather than the headline one:

| metric | weighting | largest \|change\| | most influential | n > 10% | n > 20% | median \|change\| |
|---|---|---|---|---|---|---|
| power | by effect-size count | 72.7% | MA09 | 6 | 2 | 0.5% |
| power | equal per model | **5.1%** | MA15_04 | 0 | 0 | 1.3% |
| Type M | by effect-size count | 34.3% | MA09 | 6 | 2 | 0.3% |
| Type M | equal per model | 14.3% | MA39_1 | 1 | 0 | 0.8% |
| Type S | by effect-size count | 78.0% | MA09 | 22 | 9 | 0.9% |
| Type S | equal per model | 34.9% | MA15_04 | 8 | 4 | 2.3% |

Counts are over 4 specifications x 48 models = 192 values per cell. Equal weighting
roughly halves the largest influence for every metric, but only for power does it remove
single-model leverage entirely. Type M and Type S remain sensitive under either weighting
because both diverge as the assumed effect approaches zero, so a model whose corrected
mean sits near zero moves them a long way whatever weight it carries. Note also that the
72.7% and 78.0% figures are the UWLS rows, which are supplementary; the reported
bias-robust specification is FE + VCV, where the corresponding figure is 44.7%.

This is the same conclusion `leave_one_paper_out.csv` reaches at the level of source
papers, reached independently at the level of models, and it is part of the reason the
meta-analysis-level aggregate is reported as a descriptive summary rather than as a
principal result.

`meta_analysis_level_sensitivity.csv` also carries an **equally weighted** summary
(`weighting = "equal_per_meta_analysis"`, `role = "secondary_descriptive"`), giving each
of the 48 meta-analyses equal say rather than each effect size. For power:

| | `k`-weighted | equally weighted |
|---|---|---|
| uncorrected | 0.82207 | 0.56934 |
| Yang-2023 (primary) | 0.39038 | **0.25047** |
| FE + VCV, own CR2 SE | 0.47895 | **0.53083** |
| UWLS | 0.38269 | 0.48799 |

**This is a second descriptive summary, not a robustness check on the first**, and it
does not simply confirm the `k`-weighted picture: the two weightings move the primary
and the sensitivity analyses in **opposite** directions, so the gap between them goes
from **+0.089** to **+0.280**. Equal weighting answers "what is the typical
meta-analysis like?"; `k`-weighting answers "what is the typical effect-size estimate
like?". Neither is the more correct one, which is itself part of the case for treating
the meta-analysis level as descriptive.

### The 48 model values

`model_level_metrics.csv` (192 rows: 48 models × 4 specifications) holds the per-model
power, Type M and Type S, with the assumed effect and standard error that produced
each. `figures/model_level_metrics.pdf` plots them, ordered by `k` so that leverage
reads top to bottom, with the `k`-weighted summaries drawn as rules so the gap between
a summary and the models it summarises is visible rather than asserted. Design
constraints are recorded at the top of `08_model_level_figure.R`; in particular the
figure deliberately avoids the two defects of the submitted Figure 2 (a
`scale_fill_gradient(limits = c(0, 20))` with no `oob`, censoring 11 tiles to grey that
read as missing data, and a ramp running to white that hides the 19 models at
power ≥ 0.99), and the script asserts that no plotted point falls outside its axis
limits.

**CR1 is retained as a diagnostic only** (`role = "diagnostic_crve_variant"`,
`verification_status = "single_derivation"`). It is not the specification Yang et al.
use. Note that the *qualitative* conclusions are unchanged by the variant: the CI
includes zero for **19 of 48** models under CR2 and CR1 alike, and it is the **same 19
models**.

## 4. Yang et al. (2024) implementation, in full

### What the authors actually do — verified from the primary sources

Verified 2026-08-12 by fetching both locations named in the paper's data-availability
statement, not from any earlier notes:

| source | file | md5 |
|---|---|---|
| tutorial (`yefeng0920.github.io/BiasRobustMA_tutorial/`) | `R/hands_on_R.Rmd` | `24f70634e4157e63321faac86e38f3e7` |
| analysis repo (`github.com/Yefeng0920/WLS_RVE`) | `R/final_GLS_V5.Rmd` | `ddfdf14a1840070e0956ae7c794c283e` |

The tutorial's two steps, verbatim:

```r
VCV <- vcalc(vi = var.eff.size, cluster = study, rho = 0.5, obs = obs, data = dat)
mod_MLFE     <- rma.mv(yi = eff.size, V = VCV, method = "REML", test = "t",
                       dfs = "contain", data = dat)
mod_MLFE_RVE <- robust(mod_MLFE, cluster = study, adjust = TRUE, clubSandwich = TRUE)
```

The paper's own re-analysis of the 448 meta-analyses uses the same two calls
(`final_GLS_V5.Rmd:271` and `:281`; the CRVE call there omits `adjust`, which
`robust()` ignores on the clubSandwich path).

**What `clubSandwich = TRUE` means** was read out of metafor 5.0.1's
`robust.rma.mv` source rather than inferred: the defaults are `vcov = "CR2"` and
`coef_test = "Satterthwaite"`, with `conf_test` inheriting `coef_test`. **So the
specification is CR2 with Satterthwaite degrees of freedom.** The tutorial gives the
reason in a technical note: *"CR2 correction performs better than CR1. However, CR2 is
not applicable to models with non-nested random effects … In our case, the model in the
first step does not include random effects."*

**The authors' UWLS specification is different from their FE + VCV specification**, and
this implementation follows each source call rather than copying one across. In
`final_GLS_V5.Rmd:402` and `:451` the objects carried forward are
`lm(eff.size ~ 1, weights = 1/var.eff.size)` and
`coef_test(mod, vcov = "CR2", cluster = study, test = "naive-t", target = var.eff.size)`
— CR2 with **naive-t** degrees of freedom and an explicit `target` working variance,
not Satterthwaite. (Satterthwaite variants appear at `:491` and `:497` as exploratory
code, one annotated "too conservative", and are not used downstream.) Hence
`se_method = "CR2_naive_t"` on the UWLS rows.

### External validation

The CR2 standard error has one implementation — clubSandwich's. A hand-rolled second
CR2 would not be an independent derivation, only a re-implementation tuned to agree, so
the second check is an **external anchor** instead: `06_validate_yang2024_reference.R`
runs this code path on the authors' own example data (`bird.et.al.2019.ecoletts.csv`,
md5 `1bcad390a96bdc6a8e07a81a2a31347e`) and reproduces **all nine** values in the
rendered tutorial, each to the precision it is published at:

| | ours | published |
|---|---|---|
| step one estimate | 0.073806 | 0.074 |
| step one SE | 0.017649 | 0.018 |
| step one 95% CI | [0.039190, 0.108423] | [0.039, 0.108] |
| step two CR2 SE | 0.052827 | 0.053 |
| step two Satterthwaite df | 51.96 | 52 |
| step two p | 0.168312 | 0.168 |
| step two 95% CI | [−0.032201, 0.179814] | [−0.032, 0.18] |

The prose of the paper quotes slightly different figures for the same example
(estimate 0.075, SE 0.054, t 1.375, p 0.175, CI [−0.034, 0.184]). Those come from the
448-model corpus pipeline, which applies its own exclusions, not from the tutorial's
data preparation. Both are CR2 with Satterthwaite df and they agree to two decimals;
the tutorial is the reproducible one and is therefore the anchor.

Three further cross-checks, all run on every one of the 48 models and all enforced as
hard stops in `03_yang2024_bias_robust.R`:

- **FE + VCV point estimate, two independent derivations**: closed-form `a = V^-1 1`
  weighted mean against `metafor::rma.mv`, max |difference| **6.7e-16**. UWLS closed
  form against `lm`, **8.9e-16**. (An earlier pass reported 1.1e-04 for this
  comparison; that was a data-preparation artefact of comparing two differently
  prepared frames, and it is gone now that both derivations use the identical frame.)
- **VCV construction**: `metafor::vcalc` against our own `build_vcv()`, max
  |difference| **1.4e-14**.
- **Wrapper check**: `metafor::robust(..., clubSandwich = TRUE)` against
  `clubSandwich::coef_test`/`conf_int` called directly, max |difference| **0**. This
  catches argument-passing errors in the wrapper, which is the realistic failure mode;
  it is not a second derivation of the CR2 algebra.

CR2 succeeded for **48 of 48** models (`cr2_status`, `uwls_cr2_status`). Had any failed,
the 48-model aggregation denominator would have changed, so the pipeline stops rather
than dropping a model silently.

### The rest of the implementation

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
  weights are inspectable and no optimizer is involved, and independently against
  `rma.mv` (agreement 6.7e-16). **The weights are not constrained to be positive** —
  the paper says so directly (Section 2.2.1, p.1596: "the off-diagonal elements
  become negative values") — and `n_negative_weight` records this per meta-analysis:
  25 of 48 models have at least one negative weight, up to 48% of weights in one
  model. FE + VCV therefore carries no convex-combination guarantee.
- **UWLS**: the same estimator with a diagonal VCV, i.e. `rho = 0`, so `a = 1/v_i`.
  The paper notes UWLS, Henmi-Copas and IVhet are special cases of its framework with
  zero sampling correlation. Being a convex combination, UWLS is bounded by the
  observed effect sizes.
- **CRVE** (Eq. 5), the reported specification: **CR2 with Satterthwaite degrees of
  freedom**, clustered by `study_ID`, via `metafor::robust(..., clubSandwich = TRUE)`
  exactly as the authors call it. Columns `FE_VCV_CRVE_SE_CR2`,
  `FE_VCV_CRVE_df_Satterthwaite`, `FE_VCV_CR2_ci_lb/ub`, `FE_VCV_CR2_pval`. Intervals
  are taken from `clubSandwich::conf_int` rather than reconstructed by hand.
- **CRVE, diagnostic only**: the hand-written sandwich
  `Var(beta) = (sum a)^-2 * sum_over_clusters ( sum_{i in cluster} a_i e_i )^2` with
  `e = y - beta`, as **CR0** (unadjusted) and **CR1** (small-sample factor `J/(J-1)`,
  `t` on `J-1` df). Retained so that the effect of the small-sample correction is
  visible rather than assumed; marked `single_derivation` and never reported.
- **Which SE is used downstream**: the reported meta-analysis-level rows use
  `se_source = "own_CRVE"`, `se_method = "CR2_Satterthwaite"`. The `se_beta0` pairing
  is kept as `diagnostic_hybrid_not_reported` for comparability with the published
  numbers only. At the primary-study level the SE is always the individual effect
  size's own sampling error, and **the CRVE choice cannot reach that table** — see
  section 4a.

The empirical property the method does deliver: **both estimators fell inside
`[min(observed effect), max(observed effect)]` in 48 of 48 models**. But the observed
effects straddle zero in **41 of 48**, so an estimate inside that range can still
fall either side of zero — which is why reversals remain. Bounded by the data is not
the same as unable to cross zero.

## 4a. The primary-study level, and why the CRVE choice cannot reach it

At this level `mu` is the meta-analytic estimate and `se` is **each primary effect
size's own sampling error**. That pairing *is* the design analysis: it asks what
precision an individual study had against an assumed true effect. A cluster-robust
standard error is a property of the meta-analytic mean, so substituting it here would
answer a different question. The decision to pair the Yang-2024 estimate with its own
robust SE therefore applies to the **meta-analysis-level summary only**.

A consequence worth stating because it is a free regression test: these rows are
**unchanged to the last digit** by the move from CR1 to CR2. Verified — every numeric
column of `primary_level_sensitivity.csv` differs by exactly 0 from the CR1-era file;
only the new `role` and `crit_value_method` columns were added.

`lmer(log(metric) ~ 1 + (1|study_ID))`, unweighted, 5,740 units:

| assumed effect | `role` | power | Type M | Type S | max Type M |
|---|---|---|---|---|---|
| uncorrected `beta0` | `reference_uncorrected` | 0.17354 | 2.891 | 0.02764 | 311 |
| **submitted `beta0_c3`** | **`primary`** | **0.08988** | **7.879** | **0.10207** | **26,326** |
| **FE + VCV** | **`reported_sensitivity`** | **0.13390** | **3.860** | **0.04836** | **1,029** |
| UWLS | `supplementary` | 0.12095 | 4.557 | 0.06441 | 609 |

> **Updated 2026-08-15.** This table previously read 0.17154 / 0.08774 / 0.12941 /
> 0.11659 for power. Those values used the **raw** `study_ID` as the clustering unit and
> predate the decision of 2026-08-15 to prefix it with the meta-analysis. Do not quote
> the old set.

This is the table that carries the paper's conclusion, and it supports the reading that
**the low-power finding is stable while the exact Type M and Type S values are not**:
power spans 9.0–17.4% across every assumed effect — low throughout — while Type M
spans 2.9–7.9, a factor of 2.7, and the maximum spans 311 to 26,326.

## 4b. The missing `~` in the Yang-2023 correction models

Separate from anything above: `S2_v2.R:544` and `:554` write `mod = sei + year_pub.l`
and `mod = var + year_pub.l`, omitting the `~`. `metafor` partial-matches `mod` to
`mods`, but the value is then an expression rather than a formula, so **one composite
moderator equal to the arithmetic sum of the two variables** is fitted instead of two
separate moderators. A sampling variance and a centred publication year are on
unrelated scales, so the composite has no interpretation. There are two occurrences in
`S2_v2.R` and **zero** in Yang et al.'s own script: the defect is ours.

**Yang et al. (2023) remaining the primary framework does not mean preserving this
typo.** It is an implementation error, not a methodological choice, and comparability is
not a reason to keep it.

Handled as follows. `S2_v2.R` is **unchanged**. `revision/R/reproduce/01_estimates.R` produces both
specifications: `legacy = TRUE` reproduces the submitted manuscript, and `legacy =
FALSE` — the corrected two-moderator specification — is the default and is what every
number in `revision/results/` uses.

The typo is confined to the lnRR scenario-1 block, whose only member is `MA09.csv`
(1,297 of the 5,740 rows, 126 studies). Its corrected mean is 0.1060 as written and
**0.0681** under the intended model. `beta0`, `beta1`, `beta2` and the scenario
assignment are bit-identical between the two paths, so it cannot propagate through
scenario assignment. Effect on the reported results
(`revision/results/reproduce/legacy_vs_corrected_specification.csv`):

| quantity | submitted | corrected | change |
|---|---|---|---|
| MA-level power, corrected | 0.4486 | **0.3904** | −13.0% |
| MA-level Type S, corrected | 0.01215 | **0.01223** | +0.6% |
| primary power, corrected | 0.09057 | **0.08774** | −3.1% |
| primary Type S, corrected | 0.09854 | **0.10431** | +5.9% |
| MA-level Type M, corrected | 2.040 | **2.180** | +6.9% |
| primary Type M, corrected | 7.791 | **8.124** | +4.3% |

Uncorrected results are unaffected. The power and Type S changes above are exact — the
uncorrected rows of that file change by exactly 0, which confirms it. **The two Type M
rows carry Monte Carlo noise**: that file was produced by the audit pipeline, which
matches `S2_v2.R` in drawing 10,000 unseeded deviates, and its uncorrected Type M rows
move spuriously by 0.12% (meta-analysis level) and 0.008% (primary). The Type M changes
are an order of magnitude larger than that noise, so their direction and rough size are
sound, but they should not be quoted to four digits. `revision/results/` uses the
closed form throughout and has no such noise.

Every one of these changes moves in the direction that **strengthens** the paper's
conclusion: lower power, higher Type M, higher Type S. It must be reported as a
correction in the response to reviewers, not absorbed silently.

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
- **the CR2 / Satterthwaite implementation**, against the authors' published worked
  example (all nine values), plus the three internal cross-checks in section 4;
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

- the **CR0 / CR1** hand-written sandwich and everything derived from it. These are
  diagnostics and are not reported. Note the qualitative conclusion does not depend on
  the variant: **19 of 48** FE + VCV confidence intervals include zero under CR2 and
  under CR1, and it is the same 19 models. UWLS gives 18 of 48. So roughly 40% of
  models have a bias-robust estimate indistinguishable from zero under every variant
  tried.

One number retracted by this pass. An earlier note recorded the own-robust-SE
meta-analysis-level power as **~0.479** from an unreproduced CR2 pass, and this
repository's CR1 implementation gave **0.529**; the README previously said neither
should be treated as a headline. The reproducible CR2 implementation gives
**0.47895**, derived fresh and only compared with the earlier figure afterwards. The
earlier 0.479 is therefore corroborated and 0.529 is a CR1 diagnostic, not a candidate
result.

Not carried into these tables at all: reversal counts for the Yang et al. **2023**
87-model corpus, and meta-analysis-level Type M maxima from the earlier audit passes.
Both are single-sourced and were not carried forward.

## 6. Methodological decisions still open

These are not resolved by anything in this directory. Items 5 and 6 were open when this file was written and
have since been settled; they are struck through rather than deleted so the reasoning
behind them is not lost:

1. **The magnitude-only gate**, whether to retain it. Note this is now narrower than it
   was: the reporting question above it is settled (interpretation, not arithmetic), but
   the gate is a property of the primary analysis itself. Adopting the 2024 estimators
   reduces sign reversals from 20 to 6 of 48 and does **not** resolve the `mu -> 0`
   divergence — the maximum primary-level Type M is still 1,029 under FE + VCV and 609
   under UWLS, against 26,326 as submitted.
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
~~5. How prominently to report the meta-analysis-level results.~~ **Settled**:
   a secondary descriptive summary, with the paper centred on the primary-study-level
   findings, and MA09's influence reported as a diagnostic rather than excluded. See
   the decisions list at the top and section 3.

~~6. The critical value defining the metric.~~ **Settled**: `z = 1.96` throughout, for
   both the Yang-2023 and the FE + VCV analyses, stated explicitly in the Methods. It
   is Yang et al.'s own choice (`EcoEvo_PB_script.Rmd:46-48`) and it gives both
   approaches the same decision rule. See the decisions list at the top.

Settled since this file was first written, and no longer open: whether Yang et al.
(2024) becomes a reported sensitivity analysis (**yes**, FE + VCV, with UWLS
supplementary); which standard error it is paired with at the meta-analysis level
(**its own**); which CRVE variant and degrees-of-freedom method (**CR2 with
Satterthwaite**, verified from the primary sources); and how near-zero assumed effects
are handled (**interpretation, no arithmetic rule**). Also settled: whether Type S should
be elevated over Type M as the headline for near-zero cases — **not for now**, that
would be a separate decision.

For reference, the method is **not** implemented in `orchaRd` (version 2.2.1
inspected): `pub_bias_plot()` takes already-fitted models and only adds plot layers, and
no exported function implements UWLS, FE + VCV or CRVE. Step one and step two must be
fitted directly, as done here and in the authors' tutorial.
