# Biology Open revision — final revision-analysis plan

> **Status (2026-08-10).** Approved by Ayumi as the working implementation plan.
> The audit corrections and the required sensitivity analyses are implemented in
> `R/` — see `R/README.md` for what has run and `outputs/` for results.
>
> **Explicitly left unresolved pending PI discussion:** decisions **D3** (whether
> `beta0_c3`'s one-directional gate is retained), **D4** (whether the extended
> dependence model is primary or a sensitivity), **D7** (whether a corrected
> optimistic scenario is built on `beta0_c2`), and **D8** (how the negative
> Type S CI lower limit is handled).
>
> Two items in this plan are **not** settled decisions and must not be treated as
> such: the proposed extended dependence model `(1 | study_ID) + (1 | MA_model)`
> in §4c, and the retention of meta-analysis-level power in §3a. The code
> computes the extended model only as a clearly-labelled comparison; the exact
> Yang et al. port `(1 | study_ID)` remains the reported specification.
>
> One item was added after approval: §0 now also records a genuine coding defect
> found while implementing the correction pipeline (`S2_v2.R:544,554`, the
> composite-moderator typo). See `docs/03_manuscript_corrections.md` §8.

## Context

`docs/01_audit_brief.md` asked for a forensic audit of `S2_v2.R` against Yang et al. (2023, BMC
Biology) and an evaluation of candidate revision analyses. The audit is complete (§1). Ayumi
then set seven decisions (2026-08-10) that this plan implements: investigate the sibling
`systematic_mapping_AnimCogn` repository for study identity; retain meta-analysis-level power as
a secondary analysis; drop the lognormal mean in favour of the model-based median + 95% CI; treat
the three PI sensitivity analyses as required; keep corrected leave-one-out open rather than
ruled out; reassess the dependence structure after harmonisation; and re-order everything by
dependency.

Nothing implemented. No analysis scripts written. `S2_v2.R`, the manuscript, and both
repositories unmodified — the sibling repo was read through the GitHub API and its two `.xlsx`
files inspected in the session scratchpad only.

---

## §0. Correction to the audit

My audit reported numeric/PMID-style `study_ID` values in `SMD/MA08.csv`. **That is wrong: they
are in `SMD/MA31.csv`** (212 unique values, 98% numeric). `MA08`'s identifiers are author-year
strings (`Abraham_Gruss_2010`, …). A fourth scheme I had not identified also exists: **author-only
with no year** in `SMD/des_stat/MA14_01.csv`, `MA14_02.csv`, and `zr/MA17.csv`
(`Piiroinen et al`, `Bensky_et_al`). This changes which datasets need key resolution, so it is
corrected here rather than folded silently into the plan.

---

## §1. Audit findings (established, carried forward)

Every uncorrected number in the Results reproduces from the CSVs. Full verification table and
component-by-component Yang comparison are in the audit section at the end of this file.

**Reproduced exactly:** MA-level power 82.21% (CI 71.72–94.23; lognormal mean 113.72%); MA-level
Type M 1.109 (CI 1.025–1.200; mean 1.298); MA-level Type S 0.058% (mean 0.425%); primary Type M
2.857 (mean 3.482); primary Type S 2.691% (mean 4.339%); primary power 17.15% (mean 23.10%, vs
17.20% reported — a 0.05 pp gap most parsimoniously explained by lme4/R version, testable by
rerunning in the declared environment).

**The port is faithful.** Yang's eight custom functions are character-identical; `beta0`,
bias detection, the four-scenario branch, the `beta0_c3` gate, both aggregation models, and both
back-transformations all match exactly. Three of Yang's defects came across with it: no
`set.seed` for Monte Carlo Type M; the unbounded lognormal mean; and corrected "means" computed
from the *uncorrected* dispersion (confirmed live: `0.449 × exp(0.5×0.649) = 0.6213` → the
reported 62.1%; `0.0906 × exp(0.5×0.595) = 0.1220` → the reported 12.2%).

**Ours alone:** the `zr/` vs `path = "Zr"` case mismatch (silently yields 37 models on a
case-sensitive filesystem, with nothing asserting the count); Yang's z-scaling computed but never
used, while `MMA_beta0_all_original` pools raw `beta0` across three scales; and Figure 2's "mean"
line being a different estimand from the text's.

**Documentation mismatches to fix in the manuscript:** "average sampling variance across all
effect sizes" (the code uses `se_beta0`, the SE of the pooled estimate — R2 §3.7); "linear
mixed-effects models" at Methods line 45 (the code uses weighted `lm`; line 46 is correct);
"5,740 primary studies" (they are effect-size estimates); "average" where the quantity is a
median; and the undocumented `beta0_c`/`beta0_c2`/`beta0_c3` selection and sign-and-significance
gate.

**The redundancy result (answers R1 directly).** With `t = |mu|/se`: power `= 2 − Φ(1.96−t) −
Φ(1.96+t)`; Type S is a function of `t` alone; and in Type M, `est = se·(t+Z)` makes `se` cancel.
All three metrics are deterministic re-encodings of one scalar. At the meta-analysis level `t` is
the pooled Wald statistic; within a meta-analysis `mu` is constant so the primary-study metrics
vary with `sei` alone.

---

## §2. Study-identity investigation (Ayumi decision 1)

Source: `Ayumi-495/systematic_mapping_AnimCogn`, `data/papers/primary_paper.xlsx`, sheet `data` —
1,927 rows across 48 `meta_id` values, with columns `meta_id`, `authors_id`, `primary_id`,
`author`, `title`, `journal`, `year`, `doi`, `find`, `confidence level`, `note`. Read-only.
Its `meta_id` scheme (`MA01`…) is the **same** scheme as the `MA__` filename prefixes in
`statistical_power_AnimCogn`, which is what makes the bridge possible.

**All 28 included papers are present in the mapping. 1,299 rows, 0 missing DOIs.**

### What can be established exactly, now

| Quantity | Value | Basis |
|---|---|---|
| Effect-size estimates | **5,740** | Our data, verified (5,845 raw − 105 dropped by the NA/zero-variance/year filter); identical to `sum(k)`, the aggregation weights |
| Meta-analysis models | **48** | Our data, verified |
| Source meta-analytical papers | **28** | Our data, verified; all 28 matched in the mapping |
| Unique primary studies **referenced by** the 28 papers | **1,187** | Mapping, DOI-normalized, 0 missing DOIs |
| Primary studies referenced by **>1 of the 28 papers** | **74** (one by up to 13) | Mapping, by DOI |
| Primary studies spanning **>1 effect-size metric** | **18** | Mapping, by DOI × metric |
| `(study_ID × dataset)` pairs in our data | **1,415** | Our data |
| …of which contribute **>1 effect-size row** | **903 (63.8%)**, max 115 rows | Our data |
| Unique `study_ID` labels summed per source paper | **1,260** | Our data |
| Source papers spanning >1 effect-size metric | **0** | Our data — each paper sits in exactly one metric |
| Mapping rows with `confidence level` 0 (primary list undeterminable in the source paper) | **4** (MA23, MA34) | Mapping |

The 74 / 18 duplication figures are **authoritative** — they come from DOIs, not from
surname/year inference — and they answer R2 §1.1 and the substance of R1's 1,824 → 5,740 question
without needing the row-level link.

### What cannot yet be established, and why

The mapping gives no key from our `study_ID` values to its DOIs. Comparing unique counts per
paper shows the mismatch is not uniform, and its **direction** is diagnostic:

- **9 papers have MORE `study_ID` labels than the mapping has primary studies** — MA41 +52,
  MA10 +11, MA08 +6, MA28 +3, MA05/MA12/MA45 +2, MA13/MA47 +1. A study count cannot exceed the
  number of studies referenced, so in these datasets `study_ID` is **finer than "primary study"**
  — plausibly an experiment, cohort, or treatment-arm identifier, or label variants of one paper.
  This also means `(1 | study_ID)` **under-pools** there: dependent rows from one paper receive
  separate random-effect levels.
- **17 papers have fewer**, which is expected (we retain only extractable effect sizes) — largest
  gaps MA09 −21, MA21 −21, MA39 −16, MA34 −11, MA43 −9.
- **2 papers match exactly** (MA31 212 = 212, and MA22/MA24/MA26/MA15 at 0 difference).

Four identifier schemes block a mechanical join:

| Dataset(s) | Scheme | Labels | Resolvability |
|---|---|---|---|
| `SMD/MA31.csv` | numeric, PMID-style (98%) | 212 | **High** — mapping has exactly 212 unique DOIs for MA31, so a PMID→DOI lookup should close it 1:1 |
| `lnRR/MA09.csv` | opaque codes `CD001…CD126` | 126 | **Blocked** — no `CD` string appears anywhere in the mapping's `primary_id`/`note`/`note_2` for MA09. Needs MA09's own appendix (the mapping's note names "APPENDIX S8. DATABASE OF STUDIES IDENTIFIED AS ELIGIBLE FOR INCLUSION – Table S6") |
| `SMD/des_stat/MA14_01.csv`, `MA14_02.csv`, `zr/MA17.csv` | author-only, **no year** | 20 + 17 + 22 | **Low** — surname alone cannot disambiguate; per Ayumi's constraint these must not be matched on surname |
| all others | author-year, with punctuation/format drift and mixed encodings (`Hämäläinen`, `Kverková`) | ~863 | **Medium** — within-`meta_id` matching on surname + year + title/DOI is tractable, but every match must be title- or DOI-confirmed, not surname/year-inferred |

### Reporting recommendation

Report **5,740 effect-size estimates from 48 models in 28 meta-analytical papers, which
collectively reference 1,187 unique primary studies (74 of them shared between papers, 18 across
effect-size metrics)**. Do not report a de-duplicated count of the studies *actually contributing
the 5,740 estimates* until the join is done; state plainly that identifier schemes differ across
source datasets and that a de-duplicated contributing-study count could not be established from
the available materials.

### Evidence still needed for an exact contributing-study count

1. MA09's Appendix S8 / Table S6, to decode `CD001…CD126`.
2. Publication years for MA14 and MA17 study labels (from those papers' own data files).
3. A PMID→DOI resolution for MA31's 212 numeric IDs.
4. For the 9 papers where `study_ID` is finer than "primary study", the source data's own
   definition of what a `study_ID` row is — most urgently MA41 (+52).
5. Title- or DOI-level confirmation for the author-year datasets. Surname+year alone is not
   acceptable evidence.

---

## §3. Reporting decisions (Ayumi decisions 2 and 3)

### 3a. Meta-analysis-level power — retained as secondary, narrower interpretation

**Retained. Not deleted.** Required framing:

- State explicitly in the Methods that the uncorrected meta-analysis-level quantity is a
  **monotone transformation of the pooled effect-to-SE ratio (the Wald statistic)** and is
  therefore **not an independent prospective design property**. Note that the same holds for
  MA-level Type M and Type S, since all three are functions of `t` alone.
- Present it as **retrospective and conditional**: conditional on the proxy effect and on the
  fitted precision, not a statement about how well a future meta-analysis would be powered.
- The informative content is the **corrected-vs-uncorrected contrast**, where `mu` moves to
  `beta0_c3` while `SE` is held fixed. Keep that contrast in the main text; move the uncorrected
  values to supplementary.
- **Preserve the Yang et al. comparison** (ours 82.2% → 44.9% vs ecology 55% → 36%, global change
  40–67%) as *context and precedent* — and state in the text that precedent is not the statistical
  justification, the conditional-retrospective framing is.
- Correct the abstract's "three indicators of statistical reliability", which overstates their
  independence.

### 3b. Reporting summaries compared (Ayumi decision A) — **"mean vs median" = how results are aggregated and reported**

This is issue **A** and is kept separate from the assumed-effect question in §3c. All figures below
are computed on the **uncorrected** side (`mu = beta0`), seeded, and reproduce the manuscript.

**Meta-analysis level (48 model values, current aggregation `lm(log(x) ~ 1, weights = k)`)**

| Metric | (1) Legacy lognormal mean | (2) Model-based median + 95% CI | (3a) Arithmetic mean, unweighted | (3b) Arithmetic mean, k-weighted |
|---|---|---|---|---|
| Power | **1.1372 — violates [0,1]** | **0.8221 [0.7172, 0.9423]** | 0.7100 | 0.8857 |
| Type M | 1.2957 | 1.1080 [1.0244, 1.1984] | 1.7892 | 1.1755 |
| Type S | 0.0043 | 0.0006 **[−0.0007, 0.0019]** — lower limit negative, currently floored to 0 | 0.0135 | 0.0016 |

**Primary-study level (5,740 rows, current aggregation `lmer(log(x) ~ 1 + (1|study_ID))`, unweighted)**

| Metric | (1) Legacy lognormal mean | (2) Model-based median | (3) Arithmetic mean, original scale |
|---|---|---|---|
| Power | 0.2310 | **0.1715** | 0.2266 |
| Type M | 3.4829 | 2.8579 | 3.8546 |
| Type S | 0.0438 | 0.0269 | 0.0499 |

**What each quantity is, and whether it qualifies:**

| | Estimates / describes | Respects bounds? | Same weighting as the model? | Inferential or descriptive? |
|---|---|---|---|---|
| (1) Legacy lognormal mean | the mean of a fitted lognormal *distribution* of across-unit values — not the mean of an estimate | **No.** MA-level power gives 1.1372 | **No** — k-weighted location combined with an unweighted `var(log x)` | Neither: it has no CI and falls outside the median's CI |
| (2) Model-based median | the back-transformed conditional mean of `log(x)`, i.e. a geometric-mean-type central value, with a CI from the fitted model | Power/Type S: point estimate yes; **the Type S CI's lower limit goes negative** because of the `+0.025` offset. Type M ≥ 1 respected | **Yes** — it *is* the model estimate | **Inferential** — the only candidate with a valid CI |
| (3a) Arithmetic mean, unweighted | the plain average of the 48 meta-analysis-model values, each counting once | Yes | **No** — the model weights by k | **Descriptive** |
| (3b) Arithmetic mean, k-weighted | the average effect-size-weighted value, matching the model's weighting | Yes | **Yes** | **Descriptive** |

**Numerical distance from the model-based median, and whether conclusions change:**

- MA-level power: the **weighting choice matters more than mean-vs-median**. Unweighted 0.710 vs
  k-weighted 0.886 vs median 0.822 — a 0.18 spread driven by small-`k` meta-analyses having lower
  power. Any arithmetic mean reported must state its weighting.
- MA-level Type M: unweighted 1.789 vs k-weighted 1.176 vs median 1.108. The unweighted mean is
  inflated by the heavy right tail (values > 20 in three models), which is exactly R1's "median
  2.03 vs mean 5.84" observation in the uncorrected data.
- Primary-study power: median 0.1715, arithmetic 0.2266, legacy 0.2310 — all comfortably below
  any conventional adequacy threshold. **The substantive conclusion (low power) is identical under
  all three summaries**, and the same holds for Type M (2.86 / 3.85) and Type S.
- Corrected-vs-uncorrected direction: unaffected by the choice of summary. Corrected values are
  lower for power and higher for Type M/Type S under every candidate. The **magnitude** of the
  corrected change does shift, because the legacy corrected "means" were built from the
  *uncorrected* dispersion (E3), so those specific numbers (62.1%, 12.2%, 9.49) are not comparable
  to anything and must be recomputed regardless.

**Confirmation: model-based median + 95% CI remains the most defensible principal summary.** It
is the only candidate that is an estimate from the fitted model with a valid interval, it uses the
model's own weighting by construction, and it is not distorted by the heavy right tails that
inflate the unweighted arithmetic mean. The arithmetic mean is not recommended as a replacement.

**Recommended for Supplementary Information (descriptive only):** the original-scale arithmetic
mean **alongside the minimum, quartiles, and maximum**, reported per metric and per level, with
the MA-level version given **both unweighted and k-weighted** and labelled as such. This is worth
doing: it is bounded, it needs no back-transformation assumption, it lets a reader see the heavy
tails that the log-scale median deliberately down-weights, and it answers R1's mean-vs-median
confusion directly by showing both rather than asserting one.

**Two further consequences.** The legacy lognormal mean is treated as a legacy quantity to
reproduce and diagnose, not a candidate — it is reproduced (`exp(−0.1959 + 0.5×0.649) = 1.1372`)
and diagnosed (§1, E2), then dropped. Separately, the **Type S CI lower limit is negative**
(−0.0007) and is currently floored to 0 in the manuscript ("0.06%, 95% CI: 0–0.19%"), with the
flooring done by a code comment rather than a stated method. Either state the flooring explicitly
or move Type S to a scale where the bound holds — the `+0.025` offset is what creates this.

### 3c. The >100% lognormal mean — removed

**The principal reported summary becomes the model-based median: the back-transformed log-scale
intercept with its 95% CI.** The `exp(intercept + 0.5·var(log x))` lognormal mean is dropped
everywhere. This also disposes of the E3 defect (corrected means built from uncorrected
dispersion), since the affected quantity no longer exists.

Consequential changes:

- **Results.** Every "(…; mean: X)" parenthetical is deleted — 12 of them across the two levels ×
  three metrics × corrected/uncorrected. "Average power was 82.2%" becomes "median power was
  82.2% (95% CI 71.7–94.2%)". Relabel "average" → "median" throughout, including the Abstract's
  "Average primary-study power declined from 17.2% to 9.1%" and "average meta-analysis-level
  power declined from 82.2% to 44.9%". R1's "1.11 what?" is answered by naming Type M a
  dimensionless exaggeration ratio at first use.
- **Table S1.** Regenerate in the declared environment (R 4.4.2, lme4 1.1.37, metafor 4.6.0) with
  a seeded Type M. Drop the `Mean` column; keep median + 95% CI. Retain the per-meta-analysis
  min/Q1/median/Q3/max summaries, which are bounded and legitimately descriptive. This also
  settles the 17.15% vs 17.20% question without hand-reconciliation.
- **Supplementary.** Add the descriptive original-scale table from §3b (arithmetic mean +
  min/quartiles/max; MA level both unweighted and k-weighted).
- **Figure 2.** Two changes. (i) Replace the `stat_summary(fun = mean)` crossbar in
  `scatter_plot.R` with the **model-based median and its 95% CI**, so the plotted summary is the
  same estimand as the text. (ii) Set `trim = TRUE` on `geom_violin` (or clip to `[0,1]` for
  power and Type S) so a probability display cannot visually extend outside valid bounds — with
  `trim = FALSE` the current power violin runs past 1 and below 0, which is exactly R2 §4.2's
  objection in visual form. Also state in the legend what each facet plots: the "Primary study"
  facet currently shows 48 within-meta-analysis medians, not the 5,740 rows the Results' text
  numbers come from — one label, two estimands. Verify the `MA_case` row alignment in panels (b)
  and (c), where `S2_v2.R:2927, 3081, 3119, 3157` copy the column from the *power* object.

---

### 3d. Assumed-effect sensitivity: pooled mean vs farther-from-zero 95% CI limit (Ayumi decision B)

This is issue **B** — sensitivity of the results to the assumed underlying effect `mu`. It is a
**different question** from §3b, which concerns how a given set of values is summarised. §3b varies
the *summary statistic*; §3d varies the *input* `mu` and keeps the summary fixed at the §3b
recommendation (model-based median + 95% CI).

**Scenario definition.** Baseline `mu = beta0`. Optimistic `mu =` the 95% confidence limit of
`beta0` **farther from zero** — the upper limit when `beta0 > 0`, the lower limit when `beta0 < 0`.
Not "the upper 95% CI". Verified: the sign is preserved in all 48 meta-analyses, so the scenario
never flips direction.

**Pre-implementation verification of the CI (as requested), all confirmed:**

- The interval is taken from `rma.mv`'s own `ci.lb`/`ci.ub`. Because the models are fitted with
  `test = "t"`, these are **t-based**: `beta0 ± qt(0.975, ddf) · se_beta0`. Reproduced to 9e-16.
  They are **not** `beta0 ± 1.96·se` — the largest discrepancy across the 48 is **0.825 on the
  effect-size scale**, because `ddf` ranges from **3 to 1,296** and `qt(0.975, 3) = 3.18`. Any
  reimplementation must read `ci.lb`/`ci.ub` from the fitted object, never reconstruct with 1.96.
- `beta0` and its CI are on the **same original effect-size scale**: the models are fitted on raw
  `es`, and the `es_zscore` variables are never used (A6). No back-transformation is involved.
- The interval is symmetric about `beta0` on that scale (verified).
- Scenario aggressiveness: the optimistic `mu` is a median **1.61×** `beta0`, ranging from 1.11× to
  **13.04×** — the extreme cases are meta-analyses with near-zero `beta0` and wide intervals.
  **17 of 48 CIs include zero**, so for those the scenario assumes an effect the meta-analysis
  itself cannot distinguish from nothing. Both facts must be reported with the result.

**Result (primary-study level, 5,740 rows, model-based median; Type M seeded):**

| Metric | Baseline (`mu = beta0`) | Optimistic (farther CI limit) | Absolute change | Proportional change |
|---|---|---|---|---|
| Power | 0.1715 | **0.2734** | **+0.1018** (+10.2 pp) | **+59.4%** |
| Type M | 2.8579 | **2.0341** | −0.8238 | −28.8% |
| Type S | 0.0269 | **0.0105** | −0.0164 | −60.9% |

**Distribution of study-level changes (power, percentage points):** min 0.00, Q1 +3.49,
median +7.50, Q3 +18.23, max +93.19. Effect sizes reaching power ≥ 0.8 rise from **269 to 670 of
5,740** (4.7% → 11.7%); those reaching ≥ 0.5 rise from **655 to 1,466** (11.4% → 25.5%).

**By metric (power, model-based median):** lnRR 0.1859 → 0.2831 (+52.3%, n = 1,602);
SMD 0.1697 → 0.2633 (+55.2%, n = 2,835); Zr 0.1785 → 0.3333 (+86.7%, n = 1,303). Zr is the most
sensitive, consistent with its meta-analyses having wider relative intervals.

**Substantive answer: the conclusions persist.** Under a deliberately optimistic assumed effect —
one that is on average 61% larger than the pooled estimate, and for 17 of 48 meta-analyses larger
than an effect they cannot distinguish from zero — median primary-study power is still only
**27%**, roughly **88% of effect sizes remain below 0.8**, and Type M is still **~2.0**, i.e.
significant results still exaggerate by about twofold. Low power and elevated Type M/Type S are
therefore **not artefacts of a conservative proxy effect.** This is the strongest single answer
available to R2 §3.6 and §4.3 and to the PI's question, and it should go in the main text.

Framing constraints: this is an **intentionally optimistic sensitivity scenario, never a better
estimate of the true effect**; it is a primary-study analysis; and at the meta-analysis level it is
near-vacuous by construction (setting `mu` to the CI limit forces the pooled Wald statistic up, so
power → 1), so report that there only as a stated arithmetic consequence.

**Corrected analogue — assessed, not assumed.** A coherent model-based CI **does** exist for
`beta0_c2` (the variance-moderator model's intercept: `rma.mv` returns `se`, `ci.lb`, `ci.ub`, and
`se_beta0_c2` is already extracted in `S2_v2.R`). It does **not** exist for `beta0_c3`, which is a
data-dependent selection between `beta0_c2` and `beta0` and therefore has no single fitted model
behind it. Two further empirical facts, from fitting the full `es ~ var + year_pub.l` model to all
48: **`SE(beta0_c2)` is a median 2.03× wider than `SE(beta0)`** (Q1–Q3 1.71–2.65, max 8.11; 25 of 46
exceed 2×), because the intercept is an extrapolation to zero sampling variance; and **2 of 48 did
not converge**. So an optimistic corrected scenario would be extremely aggressive and only
partially defined. Recommendation: **do not build a CI scenario for `beta0_c3`.** If a corrected
optimistic scenario is wanted, define it on `beta0_c2` restricted to the meta-analyses where the
gate actually selects `beta0_c2`, report it as conditional on that single model, and state the
SE-inflation and non-convergence. Decision deferred to D7.

---

## §4. Primary-study leave-one-out and dependence structure (Ayumi decisions 5 and 6)

### 4a. Exclusion unit

Per §2, the defensible unit is **the `(study_ID × dataset)` cluster** — all rows sharing a
`study_ID` within one CSV. It is well defined for all 48 datasets regardless of identifier scheme,
and 903 of the 1,415 clusters (63.8%) contain more than one row, so excluding a single row rather
than the cluster would leave most of the focal study's contribution in the proxy. Where §2's join
succeeds, additionally report a **paper-level** exclusion (all clusters sharing a resolved DOI
within that dataset) for the subset where the identifier is finer than "primary study" — that is
the honest unit for MA41 and the other 8 papers.

### 4b. Corrected leave-one-out — deferred, not ruled out

- **What it would estimate:** the primary-study power/Type M/Type S when the assumed effect is a
  bias-corrected pooled estimate that does not contain the focal study — i.e. it removes
  self-inclusion and publication-bias inflation simultaneously.
- **Additional computation:** the uncorrected version needs ~1,415 `rma.mv` intercept-only refits
  (hours). The corrected version needs, per exclusion, the full-moderator fit, the sign-and-
  significance test, the scenario assignment, the scenario-specific reduced fit, and the
  `beta0_c3` gate — roughly 3–4 fits per exclusion, so ~4,000–6,000 fits plus branch logic. Days,
  not hours, and it needs convergence handling per iteration (the Zr models already require
  `nlminb`).
- **How branch changes affect interpretation:** dropping one study can flip a `p < 0.05` slope
  test, moving a meta-analysis between scenarios, and can flip the `|beta0| > |beta0_c2|` gate so
  the "corrected" mean reverts to `beta0`. The resulting spread would then mix two sources —
  self-inclusion and scenario instability — and would **not** be interpretable as "the same
  analysis minus one study". If run, it must report the per-iteration scenario and gate
  assignments so branch flips can be separated from the self-inclusion effect; flip frequency is
  itself an informative robustness result about the correction procedure.
- **Whether to do it:** decide **after** the uncorrected LOO. If uncorrected LOO moves median
  primary power by little (say <1 pp on 17.2%), self-inclusion is demonstrably second-order — one
  study among `k` — and the corrected version can be declined in the response with that evidence.
  If it moves materially, the corrected version becomes necessary.

### 4c. Dependence structure — simplest defensible

Current: `lmer(log(metric) ~ 1 + (1 | study_ID_all))`, pooled across metrics, unweighted; faithful
to Yang. Three facts now bear on it: no source paper spans two metrics (so a metric random effect
is unnecessary); 18 primary studies span two metrics and are currently split; and `study_ID` is
finer than "primary study" in 9 papers, so the existing term under-pools.

**Recommendation: add one crossed term for the meta-analysis model —
`(1 | study_ID) + (1 | MA_model)` — and nothing more.** Rationale: rows from one meta-analysis
share a single `mu` by construction, which is the strongest unmodelled dependence in the design,
and 48 levels supports it comfortably. Do not add a source-paper level on top (28 levels nested
in 48, little extra information) and do not add a metric level (0 papers span metrics).

**This changes only the dependence modelling, not the estimand.** The fixed intercept still
estimates the same median log-metric across primary-study units; adding a grouping factor
reweights how units contribute and widens the CI. Report both the exact Yang port and the
extended model, and declare the extension as a deliberate deviation with its justification.
If the §2 join succeeds, re-fit with resolved study identity — that changes the *unit*, and
therefore is a genuine estimand change, so it must be reported separately rather than swapped in.

---

## §5. Implementation plan in dependency order

### A. Established corrections from the audit

| # | Item | Depends on |
|---|---|---|
| A1 | Manuscript text: "average sampling variance" → SE of the pooled estimate (R2 §3.7); delete "mixed-effects" at Methods line 45; "average" → "median" throughout; relabel 5,740 as effect-size estimates; document `beta0_c`/`beta0_c2`/`beta0_c3` and the sign-and-significance gate (R2 §3.4); explain the two effect-size routes — 43 datasets supply `es`/`var`, 5 use `escalc`, so R2 §3.3's apparent contradiction dissolves; explain why unitless metrics pool across effect-size metrics while effect sizes do not (R2 §3.2); name Type M's units at first use (R1) | nothing |
| A2 | New revision script (`S2_v2.R` untouched): fix `path = "Zr"` → `"zr"`; add `stopifnot` on 48 models and `sum(k) == 5740`; add `set.seed`; drop the lognormal mean in favour of median + 95% CI; run in R 4.4.2 / lme4 1.1.37 / metafor 4.6.0 | A1 decisions |
| A3 | Reproduce the §1 verified uncorrected numbers before trusting anything downstream | A2 |
| A4 | Check whether `MMA_beta0_all_original` (raw `beta0` pooled across three scales) reaches any figure or table; if so, restrict to within-metric or adopt Yang's z-scaling | A3 |
| A5 | Regenerate Table S1 and Figure 2 per §3b; verify panel (b)/(c) row alignment | A2, and last after all analyses |

### B. Required reviewer/PI sensitivity analyses

All use the aggregation from §4c and are reported per effect-size metric.

| Analysis | Request | Exact question | `mu` | `SE` | Unit varied/removed | Grouping / random effects | Held fixed | Recomputed | Expected output | Cost | Priority |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **B1. Bias-correction sensitivity** | R1 | Are conclusions driven by a few strongly corrected meta-analyses, or by the correction specification? | `beta0`; `beta0_c` (sei-model); `beta0_c2` (var-model); `beta0_c3` (gated) | `se_beta0` throughout (uncorrected, per C2) | nothing removed; the proxy is swapped | `lm(weights = k)` at MA level; §4c at primary level | data, all model fits, α | all six metrics under each of four proxies; plus per-meta-analysis shrinkage `\|beta0\|−\|beta0_c3\|`, which scenario each of the 48 took, and how many hit the gate | distribution of shrinkage across the 48; the four-proxy comparison; identification of the meta-analyses driving the corrected medians | Hours — all four proxies already exist in the pipeline | **Essential** — R1's central methodological question, and it also answers §2 of unresolved Q4 (gate frequency) |
| **B2. Leave-one-meta-analytical-paper-out** | PI | Are overall conclusions disproportionately driven by one source paper? | unchanged per retained model | unchanged | 1 of 28 papers, with all its models (1–4 each) | as retained | all 48 model fits — **no refitting** | the six aggregations, 28 times | 28 × 6 medians + CIs; jackknife influence per paper | Minutes | **Required** |
| **B3. Optimistic assumed-effect scenario** (§3d) | PI, R2 §3.6/§4.3 | Do low power and elevated Type M/S persist when `mu` is set near the larger plausible end of the meta-analytic CI? | the **t-based** 95% CI limit of `beta0` **farther from zero**, read from `rma.mv`'s `ci.lb`/`ci.ub` — never reconstructed as ±1.96·se | each row's own `sei`, held fixed | nothing removed; proxy replaced | §4c | `sei`, α, aggregation, data, all model fits | primary-study power, Type M, Type S under both scenarios | **already computed (§3d)**: power 0.1715 → 0.2734 (+59.4%); Type M 2.858 → 2.034; Type S 2.69% → 1.05%; per-row power change median +7.5 pp; ≥0.8 rises 269 → 670 of 5,740; per metric lnRR +52%, SMD +55%, Zr +87% | Minutes | **Required — primary-study level only.** Conclusions persist, which is the headline. At the MA level it is near-vacuous (power → 1 by construction): report as a stated arithmetic consequence only. **Never describe it as an estimate of the true effect.** Report that the optimistic `mu` is a median 1.61× `beta0` (max 13×) and that 17 of 48 CIs include zero. Corrected analogue: see D7 |
| **B4. Conventional / reference effect scenarios** | R2 §4.3, PI | How strongly do primary-study conclusions depend on the proxy? | **SMD:** d = 0.2, 0.5, 0.8. **Zr:** r = 0.1, 0.3, 0.5 converted as `atanh(r)` = 0.100, 0.310, 0.549. **lnRR:** reference response ratios — 10%, 25%, 50% → `ln(1.10)=0.095`, `ln(1.25)=0.223`, `ln(1.50)=0.405`, labelled **reference scenarios, not conventional benchmarks** | each row's own `sei` | nothing removed; proxy replaced | §4c | `sei`, α, aggregation, data | primary-study power, Type M, Type S per scenario | three scenarios × three metrics, reported separately | Minutes | **Required.** Never apply Cohen's `d` thresholds on a Fisher-`z` scale. **Do not pool benchmark results across SMD, Zr, and lnRR** — they share no common effect-size scale |
| **B5. Primary-study leave-one-out (uncorrected)** | R2 §3.6, rubric 3c | How much does self-inclusion inflate primary-study power? | `beta0` refit with the focal cluster excluded | the focal cluster's own `sei` (unchanged) | the `(study_ID × dataset)` cluster (§4a) | §4c | `sei`, α, aggregation, all other rows | ~1,415 intercept-only `rma.mv` refits, then all primary-study metrics | shift in median primary power/Type M/Type S vs the self-inclusive version | Hours | **Required** |
| **B6. Study-level descriptives** | R2 §4.1, §5.2–5.4, R1 | What does this evidence base actually cover? | — | — | — | — | — | mean/median primary-study sample size; taxa, cognitive tasks, outcomes, experimental contexts; a cited per-study table | main-text descriptives + supplementary study table; enables the construct-validity paragraph (R2 §5.3) to be interpretable | Days — largest non-statistical block; the mapping's `author`/`title`/`journal`/`year`/`doi` for all 1,299 rows is a ready starting point | **Required** — four reviewer points and three rubric items depend on it; start in parallel with A |

### C. Optional robustness analyses

| Analysis | Question | Priority |
|---|---|---|
| C1. Corrected primary-study LOO | as §4b | **Decide after B5** |
| C2. Paper-level (DOI-resolved) exclusion for the 9 papers where `study_ID` is finer than "primary study" | Does the exclusion unit change the LOO conclusion? | Optional, conditional on §2's join |
| C3. Bounded-scale aggregation (logit or beta) as a cross-check on the log-scale median | Is the median robust to the link function? | Optional — a reassurance check now that the unbounded mean is gone |
| C4. Extended vs exact-port dependence structure comparison (§4c) | Does adding `(1 \| MA_model)` change the medians or only the CIs? | Optional but cheap; falls out of A2 |

**Do not pursue:** leave-one-of-48-models-out (redundant with B2, and removing one model of a
four-model paper answers no reviewer question); any invented lnRR Cohen threshold; a pooled
small/medium/large summary across metrics; `model_est_all_corrected_original$MA.power*`
(`S2_v2.R:2449–2505`, row-mismatched `mu`/`SE` — traced as a dead end, no reported result affected,
but must not be reused); and describing anything as "robust" without naming what was varied.

### D. Decisions still requiring Ayumi / PI input

| # | Decision | Why it matters | What is needed |
|---|---|---|---|
| D1 | Whether to pursue the §2 evidence (MA09 Appendix S8; MA14/MA17 years; MA31 PMID→DOI; MA41's `study_ID` definition) or to report effect-size counts only | Determines whether a de-duplicated contributing-study count can appear at all | Effort call; items 1–3 are bounded, item 4 needs the source data |
| D2 | Where the 1,824 → 5,740 reconciliation goes — abstract or Methods only | R1 asked directly | Ayumi/PI preference |
| D3 | Whether `beta0_c3`'s one-directional gate is retained | It guarantees correction can only shrink magnitude, and collapses "no bias detected" with "correction moved away from zero" | B1 will report gate frequency; decide after |
| D4 | Whether to report the extended dependence structure as primary or as a sensitivity | Changes which CI is headline | After C4 |
| D5 | How far to soften the field-level framing (title, abstract, discussion) | R1's broadest concern; rubric 3a | Ayumi/PI — editorial, not derivable |
| D6 | DOI-archived data/code with a re-use licence (Zenodo/OSF) | R1's FAIR request and rubric 2c; A2's seed and count assertions are prerequisites | Ayumi to create the deposit once A2–A5 land |
| D7 | Whether to run a corrected optimistic scenario on `beta0_c2` (§3d) | `beta0_c3` has no valid CI; `beta0_c2` does, but its SE is a median 2.03× wider and 2 of 48 models do not converge | Ayumi/PI — decide after B1 reports gate frequency; declining it is defensible on the evidence in §3d |
| D8 | Whether the negative Type S CI lower limit is stated as floored, or Type S is moved to a bounded scale (§3b) | The manuscript's "95% CI: 0–0.19%" hides a −0.0007 limit floored by a code comment | Ayumi/PI |

### Execution order

**A1 → A2 → A3** (gate: verified numbers reproduce) **→ B1 → B2 → B5 → {B3, B4}** (one script,
swappable `mu` — they differ only in the assumed effect; do not write three pipelines) **→ C4 →
A4 → A5** (figures and Table S1 last). **B6 runs in parallel from the start.** C1 decided after
B5; D3 and D7 after B1; D4 after C4.

Note that §3b and §3d are already computed on the uncorrected side and reproduce the manuscript, so
A2/A3 amount to re-running them in the declared environment plus extending them to the corrected
side — not new derivations.

**Keep the two issues separate in the write-up.** §3b ("mean vs median") is a *reporting and
aggregation* decision: it changes how a fixed set of values is summarised, and it belongs in the
Methods' analysis-summary paragraph plus the Supplementary descriptive table. §3d ("pooled mean
effect vs farther-from-zero CI effect") is a *sensitivity analysis on the assumed underlying
effect* `mu`: it changes the input, keeps the summary fixed, and belongs in the Results as a
sensitivity subsection answering R2 §3.6/§4.3. Do not let one be presented as evidence about the
other.

---

## §6. Unresolved questions

1. **Exact contributing-study count** — bounded and localised, not open-ended: 1,187 unique
   primary studies are referenced by the 28 papers (authoritative, DOI-based), but the
   de-duplicated count of those actually contributing the 5,740 estimates needs the five evidence
   items in §2. Decision: D1.
2. **What a `study_ID` row means in the 9 papers where it is finer than "primary study"** — most
   urgently MA41 (106 labels vs 54 referenced studies). Affects the LOO exclusion unit and the
   under-pooling in `(1 | study_ID)`.
3. **Whether `beta0_c3`'s gate is retained** — B1 supplies the frequency; decision D3.
4. **The 17.15% vs 17.20% primary-power gap** — most parsimoniously an lme4/R version difference
   (my run used lme4 built under R 4.6.1; the manuscript declares v1.1.37 / R 4.4.2). Resolved by
   regenerating Table S1 in the declared environment (A5); not to be treated as a manuscript
   error before that check.
5. **Whether `MMA_beta0_all_original` reaches any figure or table** — A4 settles it.

---

## §7. Save-back status (open)

Ayumi authorized creating the AYUMI project hub and `worklog.md` for
`01_Projects/statistical_power_AnimCogn/` (currently an empty folder). **That write has not
happened** — the task prompt forbids modifying files and plan mode blocks writes — so
`AGENTS.md`'s "update the worklog before ending substantive work" is unmet and this plan file is
the only durable record. What the hub/worklog should carry: the §1 verified numbers; the A1
documentation mismatches; the three inherited Yang defects; the §0 correction (MA31 not MA08);
the §2 identity findings including the 1,187 / 74 / 18 figures and the four identifier schemes;
the §3 reporting decisions; and the §5 execution order. The vault's audit-only classification for
this repo (`Wiki/index.md:67`, `Wiki/log.md:192`) also needs updating, since that is what
conflicted with `AGENTS.md` in the first place.

---

# Appendix — full audit table

| # | Component | Yang et al. | Ours | Match? | Manuscript says | Statistical implication | Action |
|---|---|---|---|---|---|---|---|
| A1 | 8 custom functions (`power.ma_Shinichi`, `power.individual_Shinichi`, `error_S`, `error_M`, `error_M2`, `folded_es`, `folded_error`, `estimates.CI`) | `Rmd:46–113` | `S2_v2.R:15–84` | **Exact** (only `%>%`→`\|>`) | "as implemented by Yang et al." | none | none |
| A2 | Uncorrected `beta0` | `rma.mv(yi=es, V=var, random=~1\|study_ID/obs_ID, REML, test="t")` | `S2_v2.R:292–319`; `optim` for lnRR/SMD, `nlminb` for Zr | **Exact**; optimizer is a local adaptation | "the meta-analytic mean effect size" | none | none |
| A3 | Bias detection: `es ~ sei + year_pub.l`; effect present iff `p<0.05` **and** `beta0*beta1>0` / `beta0*beta2<0`; 4 scenarios | `Rmd:759–778, 932, 1016, 1053` | `S2_v2.R:429–443, 527, 591, 639` | **Exact** | full/reduced/uncorrected only | a wrong-signed slope counts as "no bias", never as evidence against bias | document (A1) |
| A4 | Effective-sample-size moderators | lnRR **and** SMD | **SMD `des_stat` only** (`:120–147, 253–256`); lnRR files lack `C_n`/`T_n` | **Intentional adaptation**, data-limited | not mentioned | 5 lnRR models keep the artefactual es–sei correlation Yang introduced `ess` to remove | state the limitation |
| A5 | `beta0_c` = sei-model intercept; `beta0_c2` = var-model intercept; `beta0_c3 = beta0_c2` if `\|beta0\|>\|beta0_c2\|` else `beta0` | `Rmd:1610–1613` | `S2_v2.R:1106–1111` | **Exact** | silent on all three | gate is one-directional; corrected distribution is downward-shifted by construction | document; report gate frequency (B1); D3 |
| A6 | z-scaling (Yang Eq. 7) + parallel scaled pipeline | present and fitted | computed at `:236–266`, **never used** | **Potentially important deviation** | code comment states the purpose | fine for unitless metrics; **not** for `MMA_beta0_all_original` (`:1131`), which pools raw `beta0` across three scales | A4 |
| B1 | Primary power: `mu` = the meta-analysis's `beta0`/`beta0_c3`, `se` = row `sei`, α=0.05 two-sided | same | `:1236, 1246`; positional alignment safe via `case = names(lnRR)` (`:397`) | **Exact** | matches | one `mu` per meta-analysis → within-meta-analysis variation is a monotone function of `sei` alone; focal study contributes to `mu` | B5 |
| B2 | Primary Type S: analytic | same | `:1278, 1287` | **Exact** | ✓ | none | none |
| B3 | Primary Type M: Monte Carlo, `N=10000`, **no `set.seed`** | `grep -c` = 0 | `grep -c` = 0 | **Exact — incl. the defect** | not mentioned | not reproducible; two runs differed by up to **2.625** on individual values, and the text highlights individual values ">20" | A2 seed |
| C1 | MA-level: `mu` = `beta0`/`beta0_c3`, **`SE = se_beta0`** (SE of the pooled estimate) | `Rmd:1760, 1769, 1783, 1803` | `:1167, 1178, 1611, 1620` | **Exact** | "**average sampling variance across all effect sizes**" | the manuscript describes a quantity neither implementation computes (R2 §3.7) | A1 |
| C2 | Corrected `mu` with **uncorrected** `se_beta0` | `Rmd:1769` + comment | `:1177–1178` + same comment | **Exact** | not mentioned | deliberate and defensible — holds precision fixed so only the proxy varies | state it |
| C3 | Combined-metric MA metrics | — | `:2449–2450` pairs `mu` and `SE` from differently-ordered frames | **Possible issue, non-propagating** | — | traced: only `beta0_c3` is taken from that object, always via `left_join(by="case")`; headline numbers come from `:2523` | do not reuse |
| D1 | MA aggregation: `lm(log(metric) ~ 1, weights = k)`, `k` = effect-size rows (`:1373`) | `Rmd:1991` | `:1397, 2529` | **Exact** | line 45 "mixed-effects models"; line 46 "weighted linear regression (lm)" | **manuscript contradicts itself**; line 46 is right. `lmer` alternative commented out at `:1384` under an unresolved `TODO` (`:1381`) | A1; resolve TODO |
| D2 | Primary aggregation: `lmer(log(metric) ~ 1 + (1\|study_ID))`, unweighted, no MA term | `Rmd:2094` | `:1513, 2631` | **Exact** | "study identity as a random effect" ✓ | per-study dependence modelled; per-meta-analysis and per-paper dependence not | §4c |
| D3 | What `study_ID` is | — | **four schemes**: author-year with format drift and mixed encodings; **numeric PMID-style (`SMD/MA31`, 212)**; **opaque codes (`lnRR/MA09`, `CD001…CD126`)**; **author-only, no year (`MA14_01/02`, `MA17`)** | **Verified from data** | "study identity" | fine within a dataset (`beta0` unaffected); across datasets the schemes block a join. 9 papers have more labels than referenced studies → `study_ID` is finer than "primary study" there | §2, §4a |
| D4 | Zr directory | — | `:98, 103` read `path = "Zr"`; disk has `zr/` | **Reproducibility issue** | — | works only on case-insensitive macOS; elsewhere silently yields 37 models, nothing asserts the count | A2 |
| E1 | Reported "median" | `exp(intercept)` | same | **Exact** | called "average" | it is a back-transformed conditional mean of `log(x)` ≈ median | A1 relabel |
| E2 | Reported "mean" | `exp(intercept + 0.5·var(log x))` | `:1406, 2539` | **Exact** | "mean: 114%" | reproduced: `exp(−0.1959 + 0.5×0.649) = 1.1372`. Unbounded; k-weighted location with unweighted dispersion; a distribution mean, not an estimate — hence outside the median's CI. Yang's power never hit the ceiling, so **precedent does not transfer** | §3b: removed |
| E3 | Corrected "mean" uses **uncorrected** dispersion | `Rmd:2015, 2119, 2170` | `:1425, 1540, 1567, 2557, 2683, 2716–2717` | **Exact — inherited defect** | reported as the corrected mean | confirmed live: `0.449×exp(0.5×0.649)=0.6213` → 62.1%; `0.0906×exp(0.5×0.595)=0.1220` → 12.2%; primary Type M ratios 1.2188 vs 1.2183 | dissolved by §3b |
| E4 | Type S offset `log(x+0.025)`, mean from `τ²+σ²` | same | same | **Exact** | not mentioned | two dispersion conventions for one purpose; code comment says "(25%)" but the offset is 2.5% | document |
| E5 | Figure 2 summary line | — | `scatter_plot.R:44–67, 134–150, 195–214`: `stat_summary(fun = mean)`, `geom_violin(trim = FALSE)` | **Different estimand from the text** | "the horizontal line denotes the overall mean" | figure mean = arithmetic mean of 48 bounded values; text mean = lognormal. Primary facet plots 48 within-meta-analysis medians, not 5,740 rows. `trim=FALSE` lets power run past 1 | §3b |
| F1 | Figure 2b/2c row labels | — | `:2927, 3081, 3119, 3157` copy `MA_case` from the power object | **Unconfirmed risk** | — | same construction order, so probably benign, but unverified | A5 |
| G1 | "28" | — | 28 distinct `MA__` prefixes; `old/` never read | **Verified** | ✓ | — | none |
| G2 | "48" | — | 48 CSVs = 48 fits; 1–4 models per paper | **Verified** | ✓ | models nested in papers; MA-level `lm` ignores the nesting | D1 |
| G3 | "5,740" | — | **effect-size rows** (5,845 − 105); equals `sum(k)` | **Manuscript mismatch** | "5,740 primary studies" | R1's and R2 §1.1's first question | A1, §2 |
