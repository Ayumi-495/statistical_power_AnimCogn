# Verification audit — PROVISIONAL record

**Status: provisional. Verified by the analyst using independent methods, NOT by
independent reviewers.** Ayumi requested a multi-agent verification panel
(Lovelace, Turing, Wald, Nightingale from the AYUMI roster). All three wave-1
agents terminated early on an account spend limit before reporting, so their
output is unusable and is not recorded here. The panel must be rerun from
scratch; the falsification briefs are in `05_verification_panel_briefs.md`.

What follows was checked by the same person who wrote the code, using
deliberately different *methods*: a different toolchain (Python, hand-computed
weighted means of logs, no `lm`/`lmer` at the meta-analysis level), closed-form
Type M instead of Monte Carlo, and reconciliation from saved artifacts rather
than from the counters inside the pipeline. That is weaker than independent
review. It did catch two errors in the analyst's own work, both recorded below.

**Nothing has been changed on the basis of this audit.** `S2_v2.R`, the
manuscript, `R/`, `outputs/`, and the reporting specification are untouched. The
B1 fix is NOT applied. The main model is unchanged. No figures or tables
regenerated. Nothing committed.

---

## 1. Confirmed coding facts

Each row was re-derived by the method named, not by re-running the pipeline.

| Claim | Pipeline | Independent re-derivation | Method |
|---|---|---|---|
| Hierarchy | 48 models / 28 papers / 5,740 rows = `sum(k)` | 48 / 28 / 5,740 | Python CSV parse with own filter |
| MA-level power median, uncorrected | 0.82207 | **0.82207** | hand-computed weighted mean of logs, no `lm()` |
| `var(log(MA.power))` | 0.649 | **0.6490** | Python |
| Models at power ≥ 0.999 | 16 of 48 | **16** | Python |
| Legacy lognormal mean | 1.13719 | **1.13719**, confirmed > 1 | Python |
| Arithmetic mean, unweighted | 0.7100 | **0.71004** | Python |
| Arithmetic mean, k-weighted | 0.8857 | **0.88571** | Python |
| MA-level Type S median | 0.0005831 | **0.0005831** | Python, offset applied by hand |
| MA-level Type M median | 1.10871 (Monte Carlo) | **1.10856** | closed-form `E[abs(X) given abs(X) > 1.96·se] / abs(mu)` |
| Primary-level Type M median | 2.8573 (Monte Carlo) | **2.8574** | closed form |
| B3 primary power | 17.154% → 27.335% | **exact match** | recomputed from `B3_per_effect_size.csv` |
| B3 threshold counts | 269→670 (≥0.8), 655→1466 (≥0.5) | **exact match** | recomputed |
| B5 medians | 17.154% → 17.224% | **17.1538% → 17.2236%** | recomputed from artifact |
| MA09 weight | 1,297 of 5,740 rows | **1,297** | artifact |

### 1a. Monte Carlo error is safe at the aggregate level, not at the row level

Closed form vs the pipeline's simulated Type M: the difference on the medians is
**0.00015** at both levels — negligible. Per row, however,
`median |MC − closed form| = 0.0082` and **`max = 1.115`**. So the medians are
stable, but the manuscript's individually highlighted Type M values "exceeding
20" are not stable to about ±1, and Figure 2b colours individual cells. A closed-
form Type M would remove this entirely and is cheaper than the simulation.

### 1b. The `~` defect is OURS, not inherited from Yang et al.

`S2_v2.R` contains exactly two bare `mod =` arguments inside `rma.mv` calls,
at **`S2_v2.R:544`** (`mod = sei + year_pub.l`) and **`S2_v2.R:554`**
(`mod = var + year_pub.l`). The only other `mod =` hits in the file are
`orchard_plot(mod = "1")` and `orchard_plot(mod = "es_type")`, where `mod` is a
legitimate string argument to a different function, not a formula.

Yang et al.'s `EcoEvo_PB_script.Rmd` contains **65** `mods = ~` and **zero** bare
`mod =` inside `rma.mv`; its only `mod =` hits are also `orchard_plot` calls.

**Conclusion: the defect was introduced in our adaptation. It is not a defect
inherited from Yang et al.** This matters for the response to reviewers, and it
corrects nothing in the earlier audit — the earlier audit already called it
"ours alone" — but it is now independently established rather than asserted.

Mechanism and effect, both re-confirmed: `metafor` partial-matches `mod` to
`mods`, but the value is an expression rather than a formula, so a single
composite moderator equal to the arithmetic sum of the two variables is fitted.
Affected: the lnRR scenario-1 block only, whose sole member is `MA09.csv`
(1,297 of 5,740 rows, 126 `study_ID` values). `beta0_c2` = 0.1060 as written,
**0.0681** under the intended two-moderator model. Every downstream direction is
as previously reported: corrected power down, corrected Type M up, corrected
Type S up.

### 1c. B5 coverage reconciles — nothing silently dropped

The concern was that `R/05_loo_study.R` has two `next` branches (datasets with
fewer than two distinct `study_ID`, and null fits) and collects into a list that
`purrr::list_rbind` drops NULLs from without counting. Reconciled independently
against the raw data rather than against the pipeline's own counter:

- `B5_leave_one_study_out_rows.csv` holds **5,740 rows** across exactly **1,415**
  distinct `(case, study_ID)` clusters.
- 5,740 is the total number of effect-size rows, so **every** row's focal cluster
  is present. A skipped cluster would have removed its rows from the artifact.
- **No dataset has fewer than two distinct `study_ID`**, so the first `next`
  branch could not fire, and no fit was null.

The reported "1,415 of 1,415, 0 skipped" is therefore correct.

*(A first Python reconciliation appeared to show 32 missing and 16 unexpected
clusters. Both were artifacts of that comparison, not of the pipeline: literal
`"NA"` strings that `readr::read_csv` parses as missing but a naive Python check
retains, and mojibake in non-ASCII `study_ID` values such as `Hämäläinen` under a
lossy decode. Recorded because it is exactly the kind of false positive a
verification pass must not propagate.)*

---

## 2. Remaining discrepancies — both are errors in the analyst's own work

### 2a. B1 concentration: wrong ranking, and it reverses the conclusion

**`R/04_sensitivity.R:54`** sorts `b1_units` by `dplyr::desc(shrinkage_pct)`, and
the statistic then sums the first five rows' `shrinkage × k`:

```r
top5_share_of_kweighted_shrinkage =
  sum((b1_units$shrinkage * b1_units$k)[seq_len(5)], na.rm = TRUE) /
  sum(b1_units$shrinkage * b1_units$k, na.rm = TRUE)
```

That takes the top five by **percentage** shrinkage — which are predominantly
small-`k` meta-analyses — and reports their share of a `k`-weighted total. It does
not answer the question it was built to answer.

| | Value |
|---|---|
| Reported (top five by percentage shrinkage) | 16.97% |
| **Correct (top five by `k`-weighted contribution)** | **58.86%** |
| The five contributors | **MA09, MA08, MA26, MA34, MA28_02** |

**The earlier conclusion — "the correction is pervasive rather than concentrated
in a few meta-analyses" — is wrong.** Roughly three fifths of the `k`-weighted
shrinkage comes from five of 48 meta-analyses. Reviewer 1's hypothesis that a few
heavily corrected estimates drive the result is closer to correct than the
earlier report stated. The largest single contributor is **MA09**, which is also
the dataset affected by the `~` defect — so the two findings interact and must be
reported together.

The statistic is at least well defined as a share: **0 of 48** meta-analyses have
negative shrinkage, because the one-directional gate guarantees
`|beta0_c3| <= |beta0|`. So there is no sign cancellation in the denominator.

**Fix NOT applied**, per instruction. The minimal correction would be to rank by
`shrinkage * k` before taking the top five, and to report both the concentration
share and the identity of the contributing meta-analyses.

### 2b. B5: "negligible" is right for the aggregate and wrong as stated

The aggregate holds: median primary-study power moves 17.154% → 17.224%
(+0.07 pp, +0.41%), Type M 2.8571 → 2.8569, Type S 2.691% → 2.694%. Verified.

But the aggregate is a **near-cancellation of offsetting shifts**, not a set of
uniformly tiny shifts. Paired per-row change in power:

| min | Q1 | median | Q3 | max | mean abs. change |
|---|---|---|---|---|---|
| **−30.06 pp** | −0.11 pp | +0.06 pp | +0.42 pp | **+29.50 pp** | 0.873 pp |

**21.88%** of the 5,740 effect sizes shift by more than 1 pp, **3.24%** by more
than 5 pp, and **0.73%** by more than 10 pp.

The defensible claim is therefore: *self-inclusion is negligible for the pooled
estimate, but not for individual studies.* The earlier phrasing overstated it.

The two sign flips also deserve reporting rather than a footnote:

| Dataset | study_ID | mu self-inclusive | mu leave-one-out | Type S self | Type S LOO |
|---|---|---|---|---|---|
| MA15_04.csv | Wu_etal_2019 | −0.1213 | **+0.5122** | 0.3161 | 0.0367 |
| MA39_2.csv | Diehl_2014 | +0.0520 | **−0.1365** | 0.4444 | 0.3576 |

Because power depends on `|mu|`, a sign flip is invisible in power but moves
Type S sharply. In MA15_04 the leave-one-out mean is four times the magnitude of
the self-inclusive one and of opposite sign — a small dataset where one cluster
dominates.

**Consequence for the earlier recommendation.** The claim that corrected
leave-one-out is unnecessary rested on the aggregate being small. That inference
is weaker than stated: the aggregate is small, but individual-level movement is
not, and the corrected pipeline additionally re-triggers scenario assignment and
the gate. This should be reconsidered rather than treated as settled.

---

## 3. Methodological choices, not coding facts

Distinguished here so they are not mistaken for defects, and none of them is
decided by this audit.

- **`se_beta0` as the meta-analysis-level precision term.** A choice matching
  Yang et al. Separate from the documentation error (the manuscript says "average
  sampling variance across all effect sizes", which is a different quantity).
- **Pairing the corrected `mu` with the uncorrected `se`.** Deliberate: precision
  held fixed so only the assumed effect varies.
- **Log-scale aggregation of a bounded probability.** Respects the lower bound,
  not the upper. This is what produced the impossible 114%.
- **"Model-based median."** `exp(intercept)` is a geometric-mean-type central
  value; it equals the median only under log-symmetry. The label is a
  simplification.
- **Exclusion unit for B5.** The `(study_ID x dataset)` cluster, forced by the
  identifier problem: `study_ID` runs four incompatible schemes and is finer than
  "primary study" in 9 of 28 papers. This makes B5 a **lower bound** on the true
  self-inclusion effect.
- **The optimistic scenario.** A t-based confidence limit, with `ddf` as low as 3,
  used as an assumed effect for 17 meta-analyses whose interval includes zero. A
  bounding exercise, never an estimate.

---

## 4. Decisions for Ayumi and Shinichi — all OPEN

### D4 — extended dependence model. New evidence; no longer cosmetic.

Primary-study power, uncorrected assumed effect:

| Model | Estimate | 95% CI | CI width | Variance components |
|---|---|---|---|---|
| `(1\|study_ID)` — current, exact Yang port | 17.15% | 16.43–17.91% | 1.48 pp | study 0.534, residual 0.149 |
| `(1\|study_ID) + (1\|case)` — extended | **22.23%** | **17.28–28.60%** | **11.32 pp** | study 0.160, **case 0.766**, residual 0.123 |

Meta-analysis identity absorbs most of the variance (0.766 against 0.160 for
study). The current interval is therefore materially overconfident: the CI widens
**7.6-fold** and the point estimate rises about 5 pp. The substantive conclusion
(low power) survives either specification; the reported number and its precision
do not.

Ayumi has asked that the extended model be checked very carefully before this is
decided. Open questions for the panel: whether `(1|study_ID) + (1|case)` is the
right crossed/nested expression given that `study_ID` labels recur across
datasets and are non-harmonised; whether 48 levels support the term; and whether
the variance reallocation is genuine structure or an artifact of the identifier
problems. **Not decided.**

### D3 — the one-directional `beta0_c3` gate

Evidence now available: 0 of 48 meta-analyses have negative shrinkage, i.e. the
gate guarantees shrinkage is non-negative by construction; the gate selects
`beta0_c2` for 34 of 48 and reverts to `beta0` for 14; and the corrected
concentration figure is 58.86%, not 16.97%. **Not decided.**

### D7 — corrected optimistic scenario on `beta0_c2`

Unchanged: a coherent model-based CI exists for `beta0_c2` but not for
`beta0_c3` (a data-dependent selection). `se_beta0_c2` is a median **1.84×** the
uncorrected SE, maximum **7.04×**. **Not decided.**

### D8 — negative lower CI limit for Type S

Unchanged: the uncorrected meta-analysis-level Type S interval has a lower limit
of **−0.0007**, floored to 0 in the manuscript by a code comment rather than a
stated method. Options are to state the flooring or move to a bounded scale.
**Not decided.**

---

## 5. Not verified by anyone yet

The panel died before reporting, so these were checked by the analyst only, or
not at all:

- line-by-line object alignment of `R/` against `S2_v2.R`'s positional matching;
- completeness of the optimizer table;
- joins, grouping variables, and factor levels, particularly the per-metric
  `dat_m` construction in `R/04_sensitivity.R`'s B4 block;
- whether the 13 reproduction-gate tolerances are honest rather than tuned to
  pass;
- whether any **undocumented** deviation from Yang et al. remains;
- whether the seed makes Type M reproducible across separate script invocations
  (each script sources `00_setup.R` afresh, and the number and order of RNG draws
  differs between scripts);
- the missingness mechanism of the 105 filtered rows;
- Nightingale's readiness gate.

These are the panel's brief. See `05_verification_panel_briefs.md`.
