# Manuscript corrections required by the audit

The manuscript file is not edited by this branch. This is the list of changes to
apply to `docs/AnimalCogn_statistical power_BiologyOpen.docx`, each with the
evidence that motivates it. Line references are to the submitted manuscript's
numbering as it appears in the reviewer comments.

Every number quoted below was re-derived from the data; see
`outputs/reproduction_check.csv` (all 13 checks pass against the submitted
manuscript) and `outputs/baseline_summaries.csv`.

---

## 1. Methods: the variance used at the meta-analysis level

**Current text (Methods, twice):** "we calculated these metrics using the average
sampling variance across all effect sizes within each meta-analysis"; and "based
on the corresponding average sampling variance".

**What the code does:** `S2_v2.R:1167` calls
`power.ma_Shinichi(mu = beta0, SE = se_beta0)`, where `se_beta0` is the standard
error of the pooled estimate returned by `rma.mv` — not an average of the primary
studies' sampling variances. Yang et al. (`EcoEvo_PB_script.Rmd:1760`) use the
same quantity.

**Why it matters:** these are different quantities, and Reviewer 2 (§3.7) asked
precisely this. The implementation is defensible; the description is not.

**Replacement:** "At the meta-analysis level we used the standard error of the
pooled estimate from each multilevel model as the precision term."

**Consequential addition (Reviewer 2 §3.6, and PI discussion pending):** because
this SE is the pooled estimate's own standard error, the resulting quantity is a
monotone transformation of the pooled effect-to-SE ratio, i.e. of the Wald
statistic. The Methods should say so, and should frame the meta-analysis-level
metrics as retrospective and conditional rather than as prospective design
properties. The same is true of meta-analysis-level Type M and Type S: writing
`t = |mu|/se`, power is `2 − Φ(1.96 − t) − Φ(1.96 + t)`, Type S is a function of
`t` alone, and in Type M the `se` cancels because `est = se·(t + Z)`.

---

## 2. Methods: which aggregation model is used

**Current text:** Methods line 45 says the meta-analysis-level estimates "were
then summarised across all meta-analyses using linear mixed-effects models";
line 46 says "weighted linear regression (lm function in base R), with the number
of effect sizes per meta-analysis as weights".

**What the code does:** `S2_v2.R:1397, 2529` use `lm(log(x) ~ 1, weights = k)`.
The `lmer` alternative is commented out at `S2_v2.R:1384` under an unresolved
`TODO` at `:1381`. Line 46 is correct; line 45 is not.

**Replacement:** delete "using linear mixed-effects models" from line 45, or make
it "using weighted linear regression (see below)". Keep line 46 as written. The
primary-study-level description ("linear mixed-effects models … including study
identity as a random effect") is accurate.

---

## 3. Abstract, Methods, and Figure 1: the unit of "5,740"

**Current text:** "48 meta-analysis models based on 5,740 primary studies"
(Abstract); "48 meta-analyses encompassing 5,740 primary studies" (Figure 1
legend).

**What the data contains:** 5,740 is the number of **effect-size rows** after the
NA / zero-variance / missing-year filter (5,845 raw − 105 dropped). It is
identical to `sum(k)`, the vector used as the aggregation weights. Verified by
`R/00_setup.R`'s `check_hierarchy()`, which asserts 48 models, 28 papers, and
5,740 rows.

**Replacement:** "48 meta-analysis models comprising 5,740 effect-size
estimates".

**Study count.** A de-duplicated count of the contributing primary studies could
not be established from the repository: `study_ID` is not a harmonised
identifier (see §6 below). What *can* be stated exactly, from
`Ayumi-495/systematic_mapping_AnimCogn`, `data/papers/primary_paper.xlsx` (all 28
papers matched, zero missing DOIs):

- **1,187 unique primary studies referenced** by the 28 meta-analytical papers;
- **74** of those referenced by more than one of the 28 (one by 13);
- **18** spanning more than one effect-size metric.

Recommended Abstract phrasing: "5,740 effect-size estimates from 48
meta-analysis models in 28 meta-analytical papers, which collectively reference
1,187 unique primary studies". This also answers Reviewer 1's question about how
the prior systematic map's 1,824 primary studies became 5,740 — the two numbers
count different things.

---

## 4. Results: "average" where the quantity is a median, and the removal of the >100% mean

**Current text:** "the average statistical power to detect the uncorrected mean
effect was 82.2% (95% CI: 71.7%–94.2%; mean: 114%)", and eleven further
parentheticals of the same form; plus the Abstract's "Average primary-study power
declined from 17.2% to 9.1%" and "average meta-analysis-level power declined from
82.2% to 44.9%".

**What the code does:** the headline value is `exp(intercept)` from a log-scale
model — a back-transformed conditional mean of `log(x)`, i.e. an approximate
median, not an average. The parenthetical "mean" is
`exp(intercept + 0.5·var(log x))`, the lognormal mean. Reproduced exactly:
`exp(−0.1959 + 0.5 × 0.649) = 1.1372` → "114%", with `var(log(MA.power)) = 0.649`
and 16 of the 48 meta-analyses at power ≥ 0.999.

**Why it must go:** power is bounded in [0,1] and this estimator has no upper
bound. It also combines a `k`-weighted location with an unweighted dispersion,
and it is the mean of a fitted distribution rather than an estimate of a mean —
which is why it falls outside the CI, as Reviewer 1 noticed. The same formula in
Yang et al. never breached the bound because their power values (55% → 36%) were
far from the ceiling, so precedent does not carry over.

**Changes:**
1. Delete all twelve "(…; mean: X)" parentheticals.
2. Replace "average" with "median" throughout Results and Abstract.
3. Report **median with 95% CI** as the principal summary.
4. Name Type M's units at first use — "a dimensionless exaggeration ratio;
   1.11 means significant estimates exceeded the assumed effect by 11% on
   average" — answering Reviewer 1's "1.11 what?".

**Comparison of the candidate summaries** (uncorrected, from
`outputs/summary_statistic_comparison_power.csv`):

| Summary | MA-level power | Bounded? | Same weighting as the model? | Status |
|---|---|---|---|---|
| Legacy lognormal mean | **1.1372** | **No** | No (k-weighted location, unweighted spread) | Removed |
| Model-based median + 95% CI | **0.8221 [0.7172, 0.9423]** | Yes | Yes | **Principal summary** |
| Arithmetic mean, unweighted | 0.7100 | Yes | No | Supplementary, descriptive |
| Arithmetic mean, k-weighted | 0.8857 | Yes | Yes | Supplementary, descriptive |

At the primary-study level: median 0.1715, legacy 0.2310, arithmetic 0.2266. The
substantive conclusion — low power — is identical under all three, so nothing in
the paper's argument depends on this choice.

**Recommended supplementary addition:** an original-scale descriptive table
(arithmetic mean plus min/quartiles/max, per metric and level, with the
meta-analysis-level version given both unweighted and k-weighted). It is bounded,
needs no back-transformation assumption, shows the heavy tails the log-scale
median deliberately down-weights, and answers Reviewer 1's mean-vs-median
confusion by showing both rather than asserting one.

---

## 5. Results: Type S confidence intervals

**Current text:** "Type S error was near zero before publication bias correction
(0.06%, 95% CI: 0–0.19%)".

**What the code does:** the log-scale model is fitted on `x + 0.025` and the
offset subtracted afterwards, which puts the interval's lower limit at
**−0.0007**. `S2_v2.R:1468` floors it to 0 in a code comment
("if the lower boundary was negative … we used 0 to replace it") rather than as a
stated method.

**Change:** either state the flooring explicitly in the Methods, or move Type S
to a scale where the bound holds. The offset is 2.5%, not the "25%" the code
comment claims. This is decision D8 in the plan.

---

## 6. Methods: what "study" means

**Current text:** "including study identity as a random effect".

**What the data contains:** `study_ID` runs four incompatible schemes across the
48 datasets — author-year strings with inconsistent punctuation and mixed
encodings; numeric PMID-style identifiers in `SMD/MA31.csv` (212 values);
opaque codes `CD001`–`CD126` in `lnRR/MA09.csv` (126 of the 164 lnRR labels); and
author-only labels with no year in `SMD/des_stat/MA14_01.csv`, `MA14_02.csv`, and
`zr/MA17.csv`. In 9 of the 28 source papers there are **more** `study_ID` labels
than the paper references primary studies (`MA41`: 106 against 54), so there
`study_ID` is finer than "primary study".

**Why it matters:** the random effect is well defined within a dataset, so
`beta0` is unaffected. But the pooled primary-study model groups labels across
datasets, where matching is correct for author-year, inert for MA31/MA09, and
demonstrably wrong in at least 20 cases (studies appearing under two metrics).

**Change:** state that study identity is as defined by each source
meta-analysis, that identifier conventions differ between source datasets, and
that a de-duplicated cross-dataset study count could not be established. Add to
the limitations.

---

## 7. Methods: the bias-correction procedure (Reviewer 2 §3.4)

Currently described only as "the full bias-correction model … a reduced version …
otherwise the uncorrected mean". Three specifics need stating:

1. **Scenario assignment uses sign, not significance.** The branch is chosen by
   the sign of `beta0 × beta1` and `beta0 × beta2`; the `p < 0.05` tests are used
   only to report which meta-analyses show a small-study or decline effect.
   Observed: scenario 1 (both slopes in the expected direction) 21 models,
   scenario 2 (drop the error moderator) 5, scenario 3 (drop the year moderator)
   20, scenario 4 (intercept only) 2.
2. **The reported corrected mean comes from the sampling-variance moderator
   model**, not the sampling-error one. `beta0_c` is the intercept of the
   error-moderator model, `beta0_c2` of the variance-moderator model, and the
   reported value is `beta0_c2`.
3. **The selection gate is one-directional.** `beta0_c3 = beta0_c2` only when
   `|beta0| > |beta0_c2|`, otherwise it reverts to `beta0`. Correction can
   shrink a magnitude but never grow one, and "no bias detected" and "correction
   moved the estimate away from zero" both collapse to `beta0`. Observed: the
   gate selects `beta0_c2` for 34 of 48 models and reverts for 14. Whether to
   retain the gate is decision D3 in the plan.

Also worth stating (Reviewer 2 §3.2, §3.3): 43 of the 48 datasets supply effect
sizes and sampling variances directly, and 5 SMD datasets are computed from raw
descriptive statistics with `escalc()` — which is why "retained in original
metrics" and "calculated from raw descriptive statistics" are both true and
neither implies conversion between metrics. And the effective-sample-size
moderator is used for those 5 SMD datasets only; lnRR has the same artefactual
effect-size/sampling-error correlation but its datasets carry no group sample
sizes, so `ess` could not be computed there — a data-limited deviation from Yang
et al. affecting 5 of 48 models.

---

## 8. A coding defect in the submitted analysis (`S2_v2.R:544, 554`)

`S2_v2.R` writes `mod = sei + year_pub.l` and `mod = var + year_pub.l`, omitting
the `~`. `metafor` partially matches `mod` to `mods`, but the value is then an
expression rather than a formula, so a **single composite moderator equal to the
arithmetic sum of the two variables** is fitted instead of two separate
moderators. A sampling variance and a latest-year-centred publication year are on
unrelated scales, so the composite has no interpretation.

The typo is confined to the lnRR scenario-1 block, whose only member is
`MA09.csv` — the largest lnRR dataset, contributing 1,297 of the 5,740
effect-size rows and 126 studies. Its corrected mean is 0.1060 under the
composite moderator and **0.0681** under the intended two-moderator model
(verified directly: the composite fit reports one moderator, the formula fit
two).

`S2_v2.R` has not been modified. `R/01_estimates.R` produces both: `legacy =
TRUE` reproduces the submitted manuscript, and the default corrected
specification is used for all revision analyses.

**Effect on the reported results** (`outputs/legacy_vs_corrected_specification.csv`).
Uncorrected results are unchanged. Corrected results move as follows:

| Quantity | Submitted | Corrected specification | Change |
|---|---|---|---|
| MA-level power, corrected | 0.4486 | **0.3904** | −13.0% |
| MA-level Type M, corrected | 2.040 | **2.180** | +6.9% |
| MA-level Type S, corrected | 1.215% | **1.223%** | +0.6% |
| Primary power, corrected | 0.0906 | **0.0877** | −3.1% |
| Primary Type M, corrected | 7.79 | **8.12** | +4.3% |
| Primary Type S, corrected | 9.85% | **10.43%** | +5.9% |

Every corrected value moves in the direction that **strengthens** the paper's
conclusion: lower power, higher Type M and Type S. This must be reported as a
correction in the response to reviewers, not absorbed silently.

---

## 9. Reproducibility items (Reviewer 1; rubric 2c)

- **No seed.** `error_M()` is a Monte Carlo estimator with `N = 10000` and
  neither `S2_v2.R` nor Yang et al. set a seed (`grep -c set.seed` returns 0 in
  both). Two runs changed individual Type M values by up to 2.625, and the
  manuscript highlights individual values "exceeding 20" and colours individual
  cells in Figure 2b. `R/00_setup.R` sets `SEED <- 20260810`.
- **Case-sensitive path.** `S2_v2.R:98,103` read `path = "Zr"` while the
  directory is `zr/`. This resolves only on a case-insensitive filesystem; on
  Linux `list.files("Zr")` returns `character(0)` and the pipeline silently
  yields 37 models instead of 48. `R/00_setup.R` reads `zr` and asserts the 48 /
  28 / 5,740 hierarchy so a silent load failure cannot pass as a valid run.
- **DOI archive.** Reviewer 1 asked for a DOI-issuing deposit with a re-use
  licence. The seed and the count assertions are prerequisites; deposit after the
  revision analyses are final.

---

## 10. Figure 2

- The violin summary line is `stat_summary(fun = mean)` over the plotted values
  (`scatter_plot.R:44-67, 134-150, 195-214`) — the arithmetic mean of 48 bounded
  values, which is **a different quantity from the text's "mean: 114%"**. Replace
  it with the model-based median and its 95% CI so the plotted summary matches
  the text.
- `geom_violin(trim = FALSE)` lets the power density extend past 1 and below 0.
  Set `trim = TRUE` or clip to [0,1] for power and Type S: an unbounded
  probability display is Reviewer 2 §4.2's objection in visual form.
- The legend should state what each facet plots. The "Primary study" facet shows
  48 within-meta-analysis medians, not the 5,740 rows the Results' primary-study
  numbers come from — one label, two estimands.
- `S2_v2.R:2927, 3081, 3119, 3157` assign `MA_case` into the Type M and Type S
  objects **from the power object**. The construction order looks identical, so
  this is probably benign, but it is unverified and should be checked before
  resubmission.
- Reviewer 1 also found the lower panels hard to read at this size; worth
  considering alongside these changes.
