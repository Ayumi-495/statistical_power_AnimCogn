# Verification panel — findings

Panel run against `05_verification_panel_briefs.md`. **Lovelace** (code and
reproducibility), **Turing** (statistics and formal reasoning) and **Wald**
(adversarial re-derivation) have all reported. **Nightingale**'s readiness gate
has NOT run, so this is not yet a completed protocol run.

> ## The panel's most consequential finding
>
> **"Power is low" survives everything. "And bias correction makes it much worse"
> does not.** Wald established that the corrected-side headlines are
> substantially artifacts of a corrected mean that is statistically
> indistinguishable from zero, and which in 20 of 48 cases has the **opposite
> sign** to the uncorrected mean. See §2m. This is a larger issue than anything
> in the earlier audit, and it bears directly on four specific Discussion
> sentences.
>
> Separately, **both** the original B1 claim and my correction to it were wrong.
> See §2e.

Full reports: Turing's is durable at
`/private/tmp/claude-501/audit_turing2/turing_report.md` with scripts
`01_ids.R`–`13_fixed.R` in the same directory. Lovelace worked in
`/private/tmp/claude-501/audit_lovelace2/`.

**Nothing has been applied.** No fix, no model change, no regenerated figure or
table, no commit. D3, D4, D7, D8 remain open.

---

## 1. Confirmed coding facts

**The `~` defect — proven, and every falsification attempt failed.** Lovelace
established structurally that `mods` is argument 4 of `rma.mv` while `...` is
argument 28, so `mod =` partial-matches `mods`; the "silently ignored" hypothesis
is refuted. Fitting `S2_v2.R:552-556` verbatim on `lnRR/MA09.csv` gives
`nrow(beta) == 2` with `rownames == c("intrcpt","mods")` and
`beta0_c2 = 0.1059904365`, against **0.0681099550** for the intended
`mods = ~ var + year_pub.l` and 0.2194470305 intercept-only. `R/01_estimates.R`'s
`fit_composite()` reproduces the defect **bit-exactly** (0.000e+00 difference), so
the legacy path is a faithful emulation, not an approximation. Scope: lnRR
scenario counts s1=1, s2=0, s3=3, s4=1; `composite_moderator_used` is TRUE for
MA09 only. `scenario`, `beta0`, `beta1`, `beta2` are bit-identical between legacy
and corrected runs, so the defect cannot propagate through scenario assignment.
Exhaustive search: two bare `mod =` in `rma.mv` in `S2_v2.R` (`:544`, `:554`);
zero in Yang's Rmd, whose `mod =` hits are all `orchard_plot`. **Ours, not
inherited.**

**Optimizer table is exactly right, including `uses_ess`.** Lovelace inventoried
all 29 `rma.mv` calls in `S2_v2.R`: null lnRR/SMD `optim`, null Zr `nlminb`
(`:314`), detect lnRR `optim` (`:341`), detect SMD `sei` **`nlminb`** (`:352`),
detect SMD `ess.sei` `optim` (`:362`), detect Zr `optim` (`:375`), and all 20
correction models `optim`. `R/00_setup.R:161-165` encodes precisely this. All 20
correction models also carry `test = "t"` and `sparse = TRUE`, matching
`fit_rma()` argument for argument.

**No object-alignment error.** Canonical order lnRR → SMD direct → SMD computed →
Zr matches `S2_v2.R:172, 1117`. `S2_v2.R`'s scenario subsetting is name-based and
`beta0_c3` returns via `left_join(by = "case")`, so the positional uses are safe.
`k = length(obs_ID)` equals `nrow` for all 48. Row counts recomputed
independently: lnRR 1602 / SMD 2835 / Zr 1303 = 5,740.

**B4's `dat_m` is correct** (`n_units` 2835/1303/1602), `mods_for()` is correct
for all 16 scenario × moderator combinations, the tibble data-masking fix is
complete, all three joins are 1:1 on verified keys, and `S2_v2.R:2927, 3081, 3119,
3157` are benign (identical `MA_case` vectors, same block order; the RHS is
merely re-levelled).

**`MMA_beta0_all_original` reaches no reported figure or table.** Its only
consumers are two orchard plots that are never saved or composed; `figure/`
contains only the three firepower panels.

**The lnRR `ess` deviation is genuine and correctly documented.** None of the five
lnRR datasets carries `C_n`/`T_n`, so `(4·C_n·T_n)/(C_n+T_n)` is not computable.
Yang does apply `ess` to lnRR (`Rmd:153, 168, 343`).

**The 105 filtered rows are not a precision-based bias — refuted in direction.**
All 105 come from four datasets (MA06_01 ×1 missing `year`; MA41_02 ×36,
MA41_03 ×11, MA41_04 ×57 missing `es`), **every one has a usable `var`**, and
median `sei` of dropped ÷ kept is 0.742 / 0.950 / 0.968 / 1.000 — dropped rows
are if anything *more* precise.

**B5 coverage and distribution reproduce exactly** (1,415 clusters, 5,740 rows,
0 skipped; −30.057 to +29.496 pp; mean |Δ| 0.8735 pp; 21.882% >1 pp, 3.240%
>5 pp, 0.732% >10 pp). **Both sign flips confirmed.**

**Non-negative shrinkage is an analytic identity, not an empirical accident.**
`shrinkage ≡ max(0, |beta0| − |beta0_c2|)`, so the `na.rm` is inert and 14 models
have shrinkage exactly 0. State it as an identity.

**Claim 2's two figures both reproduce** — 16.9702% under the coded ordering
(from MA39_1, MA26, MA01_02, MA13_03, MA06_02) and **58.8608%** ranked by
`shrinkage × k` (MA09, MA08, MA26, MA34, MA28_02). Fixing the `~` defect
**strengthens** concentration: MA09's share 17.52% → 22.07%, top five 56.455% →
58.861%. The defect does not create the concentration.

**Claim 4's arithmetic reproduces to every digit**, and **no fit is singular** —
`isSingular()` FALSE across 30 refits (5 specifications × 3 metrics × 2
scenarios). Turing found **no singular-fit warning anywhere** in the artifacts;
the brief's premise about "many singular-fit warnings" is unsupported. **No
reported CI is invalidated by singularity.**

**The `case` variance component is genuine structure, not an identifier
artifact.** `v(case)` = 0.7455–0.8477 and always dominant across four subsets,
including one with all cross-`case` label collisions removed and one restricted to
datasets whose labels all carry a 4-digit year. Crossed and nested are
numerically indistinguishable (0.22231 vs 0.22264). Both falsification attempts
fail.

**Monte Carlo error does not contaminate anything structural.** Median within-row
SD of `log(type_M)` = 0.00482 over 30 re-runs = **0.02%** of the fitted residual
variance. "Type M exceeded 20 in three models" is stable in 30/30 seeds under
both specifications. The `round(·, 3)` grain matters for only 103 of 5,740 rows.

**The 13 gate targets are transcribed correctly** — Lovelace extracted them from
`word/document.xml` independently; all 13 match, and all 12 CIs reproduce except
the primary-power-uncorrected upper limit (17.913 vs a reported 18.0).

---

## 2. Remaining discrepancies

### 2a. D4 is TWO separable decisions, and the audit bundled them (Turing)

This is the most consequential finding of the panel. The 5 pp point shift and the
interval widening have **different causes**:

| target population | point | interval | se(log) |
|---|---|---|---|
| these ~1,250 study clusters as the sample | 17.154% | 16.43–17.91% (A Wald) | 0.02207 |
| study clusters, generalising over meta-analyses | 17.15–17.56% | 14.34–22.87% (case bootstrap) | — |
| **these 48 meta-analyses weighted equally (`case` FIXED)** | **22.371%** | **21.33–23.46%** | **0.02421** |
| a superpopulation of meta-analyses (`case` random) | 22.231% | 17.28–28.60% (B Wald) | 0.12849 |

**Equal weighting of meta-analyses accounts for the entire point shift and none of
the widening** (se 0.0221 → 0.0242). The whole 5.8-fold widening comes from the
separate step of treating the 48 as a random sample. So D4 decides two things:
which unit is weighted equally, and whether the 48 are a sample or the
population. The paper's own framing ("28 published meta-analytical papers,
yielding 48 meta-analysis models") is finite-population language.

**"7.6-fold widening" is wrong on the model scale — it is 5.82× on the SE.** The
7.6 figure conflates SE inflation with the point shift.

Mechanism: because `mu` is constant within `case` **by construction**, the
between-`case` spread is a *fixed* design feature (48 different assumed effects),
not exchangeable variation. `(1|case)` does not "account for dependence" — it
changes which unit gets equal weight. Model A's implied cluster weights span only
1.4652–1.8681 (ratio 1.275), so A ≈ the unweighted mean over ~1,250 study
clusters; in B, `v(case)` dominates, so B ≈ the unweighted mean over 48
meta-analyses. **The shift is a change of estimand, not a correction.**

Extended to all metrics and both scenarios: `case` takes 62–95% of the
between-cluster variance, SE inflates 5.2–6.9×, and the point moves **up** for
power but **down** for Type M and Type S — as monotonicity in `t` requires.
Substantive conclusion unchanged under every specification (highest power bound
anywhere 31.19%; Type M ≥ 2.1 uncorrected, ≥ 4.4 corrected). But corrected
primary power's **upper** limit moves 9.07% → 13.22%, which bears on the
Discussion's comparison to "the lowest bias-corrected estimate in medicine (9%)".

### 2b. A defect in the *reported* specification, absent from the audit (Turing)

130 `study_ID` labels appear in ≥2 `case`, covering **1,098 of 5,740 rows
(19.1%)**. Of these, 118 labels (827 rows) collide only *within* one source paper
— MA22_01/02/03 share the same 30 studies, MA41_01–04 overlap — and only 12
labels (271 rows) collide across two papers. Collapsing to `case:study_ID` moves
the reported headline **+0.20 pp** to 0.17354 [0.16637, 0.18103]. Separable from
D4. Not obviously a bug: when three meta-analyses in one paper cover the same
studies, pooling under one `study_ID` is what a primary-study random effect
*should* do.

### 2c. "Model-based median" is materially wrong, not a simplification (Turing)

It is a weighted geometric mean, and log-symmetry fails badly:

| level | metric | reported | raw median | skew(log x) |
|---|---|---|---|---|
| MA | power unc | 0.8221 | **0.9129** | −1.368 |
| MA | Type S unc | 0.000583 | **7.22e−08** | — |
| primary | power unc | 0.17154 | **0.13578** | +0.706 |
| primary | Type S cor | 0.10431 | **0.13580** | — |

Off by 0.09, four orders of magnitude, 26%, and in **both directions** — not even
a consistent bias. The intervals are valid CIs for the weighted geometric mean by
monotone transformation, but not for the median, the mean, or (once the Type S
offset is subtracted) any conventional summary.

### 2d. The `k`-weighted `lm` selects for high power by construction (Turing)

`cor(log k, log se_beta0) = −0.766`: larger meta-analyses are more precise → larger
`t` → higher power mechanically. The top five by `k` (MA09, MA31, MA34, MA08,
MA26) carry **62.3%** of the weight, and their MA power values are
1.000, 1.000, 1.000, 1.000, 0.973 — **62% of the weight sits on units pinned at
the ceiling**, and 16 of 48 are ≥0.999 where the log model cannot distinguish
them. Unweighted residual SD by `k` tercile is 1.097 / 0.612 / 0.555 — variance
falls with `k` but nothing like the `1/k` that WLS assumes. Shapiro–Wilk on
`sqrt(k)·resid`: **p = 0.00041**.

Consequence: meta-analysis-level power ranges **57%–89%** on choice of weighting
and transformation alone (unweighted geometric 0.569, unweighted arithmetic
0.710, k-weighted geometric 0.822, k-weighted arithmetic 0.886), with the
reported 82.2% near the top. "Average meta-analysis-level power was 82.2%" is a
statement about the largest meta-analyses.

### 2e. B1: the audit's proposed fix is still the wrong statistic (Turing)

The exact leverage decomposition is `pull_i = (k_i/Σk)·(log m_i^cor − log m_i^unc)`,
which sums to **−0.74472 = log(0.3903771/0.8220737)** to all digits. Ranked by
|pull|: **MA26 + MA09 alone = 57.7%**, top five = **73.95%**. That ranking shares
only 3 of 5 members with the `shrinkage × k` ranking.

The diagnostic reviewers actually asked for — drop those five and re-aggregate:

| dropped | MA power unc | MA power cor | cor/unc ratio |
|---|---|---|---|
| none | 0.8221 | 0.3904 | **0.475** |
| top 5 by \|pull\| | 0.7511 | **0.5254** | **0.700** |

Dropping 5 of 48 turns a 52.5% correction-induced reduction into **30.0%**.
**Reviewer 1's concern is substantively correct at the MA level, more strongly
than 58.86% conveys.** At the primary level it is modest (+9.8%) because that
aggregation is unweighted — and that level-dependence is the answer: the
correction looks concentrated *because the MA-level summary is k-weighted*.

### 2e-bis. Wald: BOTH the original B1 claim and my correction are wrong

Re-verified by me directly after Wald raised it, from
`outputs/estimates_per_meta_analysis.csv`:

| | value |
|---|---|
| Top five hold | **2,926 of 5,740 rows = 50.98% of total `k`** |
| Their share of k-weighted shrinkage | 58.86% |
| **Excess concentration over the `k` floor** | **+7.9 pp** |
| k-weighted mean shrinkage inside vs outside the five | 0.1789 vs 0.1300 = **factor 1.38** |
| Meta-analyses that shrank | **34 of 48, covering 5,025 rows = 87.5%** |

A `k`-weighted share is dominated by `k`: five units holding 51% of the data with
exactly average shrinkage would score ~51%. So 58.86% is **not** evidence of
concentration. The correct summary is **pervasive in coverage, mildly concentrated
in magnitude**. My self-correction ("the earlier conclusion was wrong; it is
concentrated, not pervasive") overshoots in the opposite direction and is also
wrong. **The 58.86% figure must never appear without the 50.98% `k` floor beside
it.**

**And the underlying quantity is wrong for 20 of 48 units** (Wald W1).
`shrinkage = abs(beta0) - abs(beta0_c3)` understates the change in the assumed
effect wherever the corrected mean crosses zero. MA13_02: `shrinkage` = 0.00762
against an actual **|Δµ| = 0.7506**; it ranks 33rd by `shrinkage·k` and 15th by
`|Δµ|·k`. Ranking by `|Δµ|·k` changes the membership of the top five to
**MA08, MA09, MA26, MA28_01, MA28_02** (MA34 out, MA28_01 in), share 55.12%
— verified independently by me.

So there are now **three** candidate diagnostics and they disagree on membership:
`shrinkage·k` (58.86%, five members), `|Δµ|·k` (55.12%, different five), and
Turing's exact leverage `(k_i/Σk)·Δlog(metric)` (73.95%, MA26+MA09 = 57.7%). Only
the last decomposes the reported change exactly. **Unresolved; needs a decision on
which question is being asked.**

### 2m. Wald: the corrected-side conclusion does not survive

**A. The gate accepts a zero *crossing* as "shrinkage."** `R/01_estimates.R:157`
compares magnitudes only. The flagship is not MA39_1 but **MA08**: k = 531,
`beta0` = +0.3359 → `beta0_c3` = **−0.1214**, `t_c3` = **2.76** — a *significant
reversed* corrected effect, numerically perfectly stable, scored by the gate as
64% shrinkage. The abstract's own wording is "bias correction to reduce potential
overestimation"; a sign flip is not a reduction of overestimation, it is a claim
that the effect runs the other way. Of the 20 flips, **3 have |t_c3| ≥ 1.96
covering 725 rows** (MA08, MA28_01, MA41_01) and cannot be dismissed as noise.

**B. All three metrics diverge as µ → 0** (Type M → ∞, Type S → 0.5, power → α),
and the corrected analysis drives a quarter of the data into that limit. Under the
corrected effect **1,320 of 5,740 rows (23.0%) have power < 0.055** and **552
(9.6%) have Type S > 0.45**, against 177 and 8 uncorrected. Minimum corrected
power is exactly 0.05000000; maximum corrected Type S 0.4998962; maximum
corrected Type M **26,326**. **37 of 48 corrected means have a CI including
zero — 3,158 of 5,740 rows, 55%.**

**The manuscript's most quotable corrected result is entirely mechanism B.**
Recomputed closed-form for all 48, exactly three MA-level Type M exceed 20, and
all three are sign flips whose corrected CI spans more than the uncorrected
effect:

| MA | k | `beta0` | `beta0_c3` | `t_c3` | Type M | CI of `beta0_c2` |
|---|---|---|---|---|---|---|
| MA39_1 | 37 | −0.11816 | +0.0000649 | 0.00055 | **4,259.6** | [−0.4605, +0.4606] |
| MA26 | 459 | +0.22590 | −0.001757 | 0.030 | **77.5** | [−0.2621, +0.2586] |
| MA13_03 | 14 | +1.09765 | −0.029896 | 0.095 | **24.5** | [−1.1051, +1.0453] |

Type M = 4,260 is derived from an estimate consistent with any true effect in
±0.46. It measures 1/noise, not exaggeration. **Figure 2b's "values > 20 shown as
20+" means no reader ever sees 4,260.**

**How much of the corrected headline is this?**

| Variant | primary power | primary Type M | primary Type S | MA n>20 |
|---|---|---|---|---|
| uncorrected | 0.17154 | 2.857 | 0.02691 | 0 |
| **as reported** | **0.08774** | **8.123** | **0.10431** | **3** |
| floor \|µ\| at 0.5·se (no correction withdrawn) | 0.08832 | **6.017 (−26%)** | 0.10011 | **0** |
| drop 9 MAs with \|t_c3\| < 0.5 | 0.09516 | **5.059 (−38%)** | 0.08298 | — |
| sign-preserving gate (**confounded**) | 0.13463 | **3.711 (−54%)** | 0.04265 | **0** |

Bounding the denominator alone — changing no sign, withdrawing no correction —
costs a quarter of the corrected Type M and eliminates all three ">20" models;
flooring at just 0.25·se already gives zero (max 9.42). The sign-preserving gate
moves most but is **confounded**: it withdraws the correction for 20 of 48 MAs
(1,932 of the 5,025 shrunk rows), so part of its effect is simply "less correction
applied" — that confound must be stated wherever the number appears. Three
independent variants bracket corrected Type M at **5.06–6.02** against the
reported 7.79/8.12.

**Specific Discussion claims that break:**
- ms 56, "Type M exceeded 20 in three models after bias correction" → **zero**
  under the floor variant and under the sign gate. Entirely flip-driven.
- ms 52/57, primary "some values exceeding 20" → 796 rows (13.87%), of which
  **520 (65%) lie in the 20 flipped MAs**. 276 remain; a separate claim that
  survives partially.
- ms 58, primary Type S "9.85%… exceeding typical estimates in ecology and
  evolution (5–8%)" → at 4.27% the **comparative claim reverses**; at 7.26%
  (drop-five) or 10.01% (floor) it sits inside or at the edge of the comparator.
  **The cross-field comparison is not robust.**
- ms 57, primary Type M "7.79 vs eco-evo 2.5 to 4.0" → at 3.71 it lands *inside*
  the comparator range; at 5.06–6.02 it is above but far less dramatic.

Wald's recommended rule, which is data-driven rather than a chosen threshold: do
not report Type M for any meta-analysis whose corrected mean is not
distinguishable from zero, using `ci_lb_c2`/`ci_ub_c2` already in the outputs.
That removes **37 of 48** corrected Type M values — which is itself the finding.

### 2n. Wald: three source CSVs are latin-1, not UTF-8

`SMD/MA24.csv`, `SMD/MA26.csv`, `SMD/des_stat/MA18.csv`; and
`outputs/B5_leave_one_study_out_rows.csv` and `outputs/B3_per_effect_size.csv` are
not valid UTF-8. `readr::read_csv` assumes UTF-8, so `study_ID` values with
diacritics are mangled, and because `summarise_primary` pools `study_ID` across
datasets a mangled label **splits** a level a clean label would join. Impact is
small (it only ever splits, never merges) but it makes the artifacts non-portable
and it is what produced my own false-positive "32 missing clusters". Resolving
check: `read_csv(..., locale = locale(encoding = "latin1"))` for those three files
and re-tabulate `n_distinct(study_ID)`.

### 2f. Missing uncertainty larger than D4, discussed nowhere (Turing)

`mu = beta0` is treated as **known** at the primary-study level; `se_beta0` is
discarded. So every primary-level interval reflects only the spread of `sei`
across rows. For the **17 of 48** meta-analyses whose CI includes zero, `mu` is
not distinguishable from 0 and power is not distinguishable from 0.05.

### 2g. Type S: the point estimate is not identified without the offset (Turing)

`exp(β₀) = 0.025583` against an offset of **0.025** — the reported 0.00058 is the
difference of two numbers agreeing to three significant figures. The log-scale CI
is ±5.3% of `exp(β₀)`; after subtraction it becomes **[−222%, +234%]** of the
point estimate. Across defensible offsets (1e−6 to 0.1) the MA point estimate
spans a factor of **~490** and the primary one a factor of ~7; the lower limit
turns negative at any offset ≥ 0.001 at MA level. Cause: **40 of 48** MA-level
Type S values are below 10% of the offset (raw median 7.22e−08), and at primary
level 56.2% are below the offset — which was introduced to handle **one** exact
zero out of 5,740 and none out of 48.

Yang's own comment attributes the negative limit to "the negative variance"
(`Rmd:2029` and eleven other sites). **That diagnosis is wrong**; the cause is
subtracting an offset larger than the estimand. Do not repeat it.

### 2h. The manuscript's twelve "mean" values use three different formulas

Both reviewers independently. `S2_v2.R` uses `var(log(uncorrected))` for both
power means and both primary Type M means, `var(log(corrected))` for MA Type M
and MA Type S, and `σ²_study + σ²_resid` of the **uncorrected** fit for primary
Type S. Turing: **8 of 12 reproduce within 2%; 4 do not** — MA power corrected
(+29.8%), primary power corrected (−10.0%), primary Type M corrected (+194%),
primary Type S corrected (+24.7%) — and all four reproduce exactly under Yang's
uncorrected-dispersion rule. **The 13-check gate tests exactly 1 of the 12.** So
"the legacy path reproduces the manuscript on all 13 checks" is true, while "the
legacy path reproduces the manuscript" is false.

*Reviewer disagreement, unresolved:* Lovelace reports the 13.4% (primary Type S
corrected) reproduces under **no** variant of nine tried (required variance
0.5047); Turing reproduces it as 0.1378 vs 0.134 under the uncorrected
`varcor + sigma²` rule. Turing's is within ~3%, Lovelace's tolerance was tighter.
**Wald to arbitrate.**

### 2i. Gate tolerances do not discriminate (Lovelace)

Measured Monte Carlo SD of the Type M medians is **0.00181** (MA, 30 reps) and
**0.00080** (primary, 8 reps), so tolerances of 0.20 and 0.40 are **110× and
500×** the real error. Of the six checks that could distinguish the legacy
specification from the corrected one, **four do not** — MA Type M corrected
(2.1799 passes, 77 SD off), primary Type M corrected (8.1238 passes, 417 SD off),
MA Type S corrected, and primary power corrected (passing by 5%). Only MA power
corrected and primary Type S corrected actually bind. Separately,
`MA type S, uncorrected` has a tolerance of 0.0005 on a reported value of 0.0006 —
83% of the value.

### 2j. The seed does not make Type M reproducible across scripts (Lovelace)

| quantity | `03_baseline.R` | `04_sensitivity.R` | diff |
|---|---|---|---|
| MA Type M, uncorrected | 1.10870628 | 1.10738912 | −1.3e-03 |
| MA Type M, corrected | 2.17987249 | 2.17612608 | −3.7e-03 |
| primary Type M, uncorrected | 2.85755132 | 2.85718159 | −3.7e-04 |

`primary Type M, corrected` agrees to 8 dp only by arithmetic coincidence — both
scripts reach that call at RNG offset **17,412**. `set.seed()` in `00_setup.R`
pins Type M only when the number of preceding draws happens to match, and
`04_sensitivity.R:109` computes a *third* independent realisation of the baseline
Type M. **`R/README.md`'s claim that the seed fixes reproducibility is wrong.**

### 2k. 20 of 48 corrected means flip sign (Lovelace)

Covering **1,932 of 5,740 rows (34%)**. The gate compares magnitudes only — no
sign constraint at `R/01_estimates.R:157-158`, `S2_v2.R:1108-1111`, or Yang
`Rmd:1610-1613`. Examples: MA08 0.3359 → −0.1214; MA26 0.2259 → −0.00176;
MA28_02 −0.7055 → +0.1548; MA41_01 0.2086 → −0.1310. Not near-zero jitter: 17 of
20 have |`beta0_c3`| > 0.05 and 13 have > 0.10. All metrics use `|mu|`, so nothing
errors. Related: **MA39_1**'s `beta0_c3` = 6.49e−05 (t = 0.00055) yields a Type M
of **4,164–4,304** depending on seed — and this is one of the three "values over
20" the manuscript reports. **D3 evidence.** Wald is investigating whether the
manuscript's Type M claims are driven by near-zero corrected denominators.

### 2l. Smaller items

- **`R/00_setup.R:9-11` and `R/README.md` state the wrong failure mode** for the
  `zr` path: with `zr/` unreadable, `length(Zr) == 0` makes `1:length(Zr)` equal
  `c(1, 0)` and `Zr[[1]]` raises "subscript out of bounds" at **`S2_v2.R:189`**.
  The failure is **loud, not silent**; 37 models is never reached. The fix remains
  correct; only the stated rationale is wrong — and it is the text most likely to
  be reused in a reviewer response.
- **`outputs/legacy_mean_diagnosis.csv` is an orphan** — no script in `R/` writes
  it; it is stale from a superseded version of `03_baseline.R`.
- **`stopifnot(identical(names(L), idx$case))`** at `R/01_estimates.R:64` and
  `R/05_loo_study.R:37` cannot fail: both sides derive from the same expression.
  `R/02_metrics.R:54` *is* load-bearing.
- **`confint` default differs.** `S2_v2.R:1525, 2641` and Yang `Rmd:2105` call
  `confint(fit)` with no method, and `lme4` defaults to `"profile"`; `R/` uses
  `"Wald"`. Numerically identical here (16.4276–17.9134 both ways), but
  undocumented — and `S2_v2.R:1381` carries the author's own TODO about it.
- **`nlminb` is a deviation from Yang**, who uses `optim` in all 69 `rma.mv`
  calls. The Zr use is commented at `S2_v2.R:311`; the SMD-detection use at `:357`
  is undocumented in `S2_v2.R`.
- **Latent, currently inert:** `S2_v2.R` has no lnRR scenario-2 and no Zr
  scenario-4 block. Both are empty under the present data, but would silently drop
  datasets if either became non-empty.
- **`fit_rma` swallows errors into `NULL`** (`R/00_setup.R:175`), caught by
  `stopifnot` for null/detection stages but only by a `warning` at the correction
  stage.
- **The 17.20% vs 17.154% discrepancy is explained, and it is transcription.**
  `var(log power) = 0.59565`, `exp(0.5v) = 1.34692`, and
  0.17154 × 1.34692 = 0.23105 → the reported 23.1%, whereas
  0.1720 × 1.34692 = 0.23167 → 23.2%, which the manuscript does not report. The
  manuscript's own pair is internally inconsistent and 17.154 is the consistent
  value, so the lme4-version hypothesis is unnecessary. **The 18.0% upper limit
  remains unexplained by anyone.**

---

## 3. Methodological choices, not coding facts

- **The corrected MA-level metric is a hybrid.** Pairing `beta0_c3` with its own
  SE gives MA power **0.3096**, not the reported 0.3904 (medians 0.1286 vs
  0.2614) — because `se_beta0_c2` is a median 1.84× and up to 7.04× the
  uncorrected SE. The reported quantity is "the precision we had, applied to a
  smaller effect": the power of no fitted model. Defensible as a counterfactual
  isolating `mu`, and it is Yang's choice, but the natural reading — "power after
  accounting for publication bias" — would need the corrected SE and is ~21%
  lower.
- **All three metrics are one scalar.** Verified numerically: holding `t` fixed and
  varying `se` over {0.01, 0.1, 1, 10} leaves power, Type S and Type M unchanged
  to 10 dp. Power increases in `t`; Type S and Type M decrease. So row rank
  orderings, all threshold counts, and any claim of the form "low power **and**
  high Type M **and** elevated Type S" are **one finding stated three times**. The
  three reported *summaries* are not redundant, because a weighted mean of logs is
  not equivariant under the monotone maps.
- **MA-level power is the pooled Wald statistic relabelled.** Spearman correlation
  with `|t|` is **1.000000** exactly; with the normal-approximation p-value
  −0.999729. Median `|t| = 3.3198` → median p ≈ 0.00091. It answers "how
  decisively did this meta-analysis reject the null?", not "was this design
  adequate?". 31 of 48 have p < 0.05; 29 of 48 have MA power > 0.8.
- **Yang implies model A — as an omission, not an argument.** Yang fits
  `lmer(log(power) ~ 1 + (1|study_ID))` per effect-size type (`Rmd:2095, 2524,
  2955`) with no meta-analysis random effect anywhere, and their `study_ID` is
  likewise unlisted across their own meta-analyses. Inheriting A is faithful; it is
  not evidence A is right.
- **The B5 lower-bound claim is sound in direction, and should not be
  quantified.** It holds in expectation, not per cluster. The focal cluster's share
  of its dataset's rows has median 1.54%, mean 3.39%, and 81.8% of clusters remove
  <5%. A 2× coarser identifier would roughly double that; scaling mean |Δlog| from
  4.1% to ~8% is the order of magnitude, but the map from removed fraction to `Δmu`
  is non-linear. **Do not put a single number on it.**
- **The ±30 pp B5 tails are real, not a near-zero-`mu` artifact.** In the 22
  clusters with >10 pp change, median `|mu_loo|` = 0.5860 and median `|mu_self|` =
  0.5791, both far from zero, and in 14 of 22 the LOO mean is *farther* from zero.
  Only 1 of 1,415 clusters has `|mu_loo|` < 0.02, and its maximum change is
  0.11 pp. Genuine small-dataset sensitivity.
- **The two sign flips do not invalidate those clusters' numbers**, only the
  *interpretation* of Type S, which presupposes a known sign. Exclude from any
  highlighted-value list; do not call them invalid.
- **B5's aggregate is a cancellation, not smallness.** 61.60% of rows move up,
  38.24% down; mean Δlog +0.00900 against mean |Δlog| 0.04125. Three distinct
  shift measures: reported `lmer` summary **+0.407%**, raw row-level geometric mean
  **+0.904%**, typical row **+4.21%** — a factor of ~10. The reason corrected LOO
  should be reconsidered is that the uncorrected aggregate is *uninformative*
  about the corrected case (which re-triggers scenario assignment and the gate),
  not that the effect is small.

---

## 4. Decisions for Ayumi and Shinichi

### D4 — now two questions, not one

1. **Which unit is weighted equally?** Study clusters (17.15%) or meta-analyses
   (22.37%)? This accounts for the whole point shift and none of the widening.
2. **Are the 48 meta-analyses a sample or the population?** Fixed `case` gives
   [21.33, 23.46]; random `case` gives [17.28, 28.60].

Four coherent packages, each with a consistent interval, and each answering a
different question. Whichever is chosen, the wording should name the unit — "the
typical effect size's study" versus "the typical meta-analysis" — rather than
"average primary-study power", which does not distinguish them. **Settled by the
panel and no longer open:** 48 levels do support the term (no singular fit in 30
refits), the variance reallocation is genuine (not an identifier artifact), and
crossed versus nested is immaterial. What remains open is only the population
target. **Not decided.**

### D3 — the one-directional gate

New evidence: shrinkage non-negativity is an **analytic identity**; 34 of 48 gate
to `beta0_c2` and 14 have shrinkage exactly 0; **20 of 48 corrected means flip
sign**, covering 34% of rows, with no sign constraint anywhere; and MA39_1's
near-zero corrected mean produces a Type M of ~4,200 that is reported among the
manuscript's three "values over 20". **Not decided.**

### D7 — corrected optimistic scenario on `beta0_c2`

Unchanged, with the relevant magnitude now quantified: `se_beta0_c2` is a median
1.84× and up to 7.04× the uncorrected SE, and pairing `beta0_c3` with its own SE
would move corrected MA power from 0.3904 to 0.3096. **Not decided.**

### D8 — Type S bound and offset

Now much larger than a negative lower limit: the point estimate is **not
identified** without the offset convention (factor ~490 across defensible
offsets at MA level). Options laid out by Turing, none recommended: report offset
sensitivity as the headline diagnostic; drop the log model at MA level for raw
median + IQR; use a logit scale (`qlogis(x)` gives 0.00148 [0.00065, 0.00335]);
report threshold counts (1,601 of 5,740 rows exceed 0.05, 898 exceed 0.1, 232
exceed 0.25, maximum 0.4912); or keep the offset and state the flooring as a
method. **Not decided.**

---

## 4-bis. Additional D3 evidence from Wald (still not decided)

A sign-preserving gate (`|beta0| > |beta0_c2|` **and** same sign) would select
`beta0_c2` for **14 of 48** and withdraw the correction from 1,932 of the 5,025
shrunk rows. With that confound stated: MA power 0.39038 → 0.60839; MA Type M
2.1771 → 1.3093; MA Type S 0.01223 → 0.00185; primary power 0.08774 → 0.13463;
primary Type M 8.1231 → 3.7114; primary Type S 0.10431 → 0.04265; MA n>20 3 → 0;
primary rows >20 796 → 276.

The mechanism-isolating alternative withdraws **no** correction: floor |µ| at
0.5·se → primary Type M 6.0174, power 0.08832, Type S 0.10011, MA n>20 = 0.

D8 note: under the sign-preserving gate the MA-level Type S lower limit is
**−0.000883**, slightly *worse* than the current −0.000710 — so the flooring
problem is not an artifact of the current specification.

D7 note: 37 of 48 `beta0_c2` intervals include zero, so a "farther-from-zero
confidence limit" built on `beta0_c2` would usually be a limit on an estimate
consistent with no effect, with an arbitrary sign for the 20 flipped cases.

## 5. Still outstanding

- **Nightingale**'s readiness gate has not run, so the protocol run is incomplete.
- **The 13.4% is unreproduced by all three reviewers.** Wald's resolving check is
  not another formula variant but `git log --follow` on the four data directories
  and `S2_v2.R` against the .docx modification date, to test whether the twelve
  parentheticals were generated from different data states. If so, no single
  formula will ever reconcile them.
- **Which B1 diagnostic answers Reviewer 1's question** — three candidates
  disagree on membership (§2e-bis).
- **Wald's five could-not-check items**, notably: whether the 104 missing `es`
  values in MA41 are missing at random *with respect to effect magnitude* (all
  three reviewers flag this; the `sei` test is null but insufficient, and MA41_01,
  MA41_02 and MA41_03 are all in the 20 sign-flipped set); whether Figure 2b's
  "overall mean" line uses the arithmetic mean (93.86 corrected) or the lognormal
  mean (5.86); and a cheap corrected-LOO probe on the 12 datasets with
  `n_study_id < 10` (~120 clusters × 4 fits) to measure how often the gate or
  scenario flips when one study is removed.
- **Whether the 104 missing `es` values in MA41_02/03/04 are missing at random**
  with respect to the source meta-analyses' reporting. Both reviewers flag it as a
  first-order threat that no code check can settle; Lovelace showed it is *not* a
  precision effect, which is as far as these files go.
- **The manuscript's 18.0% CI upper limit** is unexplained by anyone: not the
  identifier collapse (0.18103), not profile-versus-Wald (0.17913). Needs the
  declared environment or the original session objects.
- **Table S1 was not available** to the panel, and it is where the 17.20% question
  would be settled definitively.
