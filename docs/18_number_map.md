# What every number means, and which file holds it

Verified 2026-08-15. Read this with `results/revision/README.md`, which explains the
*methods*; this file explains the *numbers*.

---

## 0. The verification, in one place

| check | what it compares | result |
|---|---|---|
| **Reproducibility** | full `run_all.R` from scratch, against the previous outputs | **all 14 CSVs and Table S1 byte-identical** |
| **Manuscript text** | every number in `docs/17`, against the result files (`21_audit_manuscript_claims.py`) | **90 of 90 claims agree** at the precision they are written to |
| **Second derivations** | 11 checks in `20_verify_reported_numbers.R` | **11 of 11 pass** → `verification_audit.csv` |
| **Scenarios** | all 123 rows re-derived from scratch in Python (`17_verify_scenarios.py`) | **123 of 123**, max 1e-6 relative |
| **Yang-2024 CR2** | the authors' published worked example (`06_...`) | **9 of 9 published values** |
| **Frozen files** | `S2_v2.R`, `outputs/`, `figure/` | md5 `a910b4cc…` unchanged, 32 files unchanged |

**One error was found and fixed.** `docs/17` said equal weighting changed corrected
meta-analysis-level power to **25.1%**. The stored value is 0.250474, which is **25.0%**.
The 25.1% came from rounding Table S1's already-rounded 0.2505 a second time. This is
exactly the class of error the manuscript-claim audit exists to catch, and it is the only
one it found.

**Nothing else changed.** The pipeline was already correct; this pass confirmed it and
closed four gaps in *coverage* — the bias-correction gate, the metric definitions, the two
model-level estimands, and the meta-analysis-level Part A rows — none of which had a
second derivation before today.

---

## 1. The headline numbers

### Primary-study level — this is the paper's main result

5,740 effect-size observations, each judged against **its own** sampling standard error.
Summary = back-transformed intercept of `lmer(log(metric) ~ 1 + (1|cluster))`.

| | uncorrected | Yang-2023 corrected | Yang-2024 bias-robust |
|---|---|---|---|
| **power** | 17.4% | **9.0%** | 13.4% |
| **Type M** | 2.89 | **7.88** | 3.86 |
| **Type S**, model-based | 2.76% | **10.2%** | 4.84% |
| **Type S**, raw median | 1.8% | **13.6%** | 5.2% |

**What it means.** A typical primary study in this corpus had about a 1-in-6 chance of
detecting the effect its own meta-analysis estimates, and about a 1-in-11 chance once
publication bias is corrected for. When such a study *did* reach significance, its effect
was on average 2.9 times too large uncorrected, 7.9 times too large corrected.

**File:** `results/revision/primary_level_sensitivity.csv` (36 rows = 4 assumed effects ×
3 metrics × 3 estimands). Also Table S1 part A.

### Meta-analysis level — a descriptive summary, not a principal result

48 models, each judged against the standard error of its own pooled estimate. Summary =
`k`-weighted geometric mean.

| | uncorrected | Yang-2023 corrected | Yang-2024 bias-robust |
|---|---|---|---|
| **power** | 82.2% | **39.0%** | 47.9% |
| **Type M** | 1.11 | **2.18** | 1.53 |
| **Type S**, raw median | 7×10⁻⁸ | **0.22%** | 8×10⁻⁷ |
| **Type S**, model-based | 0.06% | **1.22%** | 0.56% |

**What it means, and its main limitation.** The assumed effect and the standard error it
is judged against are outputs of the same fitted model, so these metrics restate that
model's own effect-to-uncertainty ratio rather than adding information. The estimand is a
descriptive summary of the 48 models included.

**File:** `results/revision/meta_analysis_level_sensitivity.csv` (39 rows; 24 reportable,
15 marked `diagnostic_*` and not for quoting).

---

## 2. Why the corrected numbers are what they are

Three facts do all the work, and each is worth stating in the paper.

**The correction gate compares magnitudes only, and is blind to sign.** It replaces the
uncorrected mean when the corrected one is smaller in absolute value. It selected the
corrected mean in **34 of 48** models, and **20 of those 34 reversed sign**. Those 20
models carry **1,932 of the 5,740** observations. Under the bias-robust estimator the
reversal count falls to **6 of 48**.
→ `reversal_counts.csv`, and `per_meta_analysis_estimates.csv` model by model.

**Type M diverges as the assumed effect approaches zero.** Type M → **2.3378/|t|** as
t = |μ|/se → 0, verified to 1.3e-13. So a corrected mean near zero produces an arbitrarily
large Type M — this is a property of the definition, not an anomaly, and it is why three
models exceed Type M = 20 after correction (none do before) and why the primary-study
maximum reaches 26,326. It also inverts intuition: the *smallest* reversals give the
*largest* Type M.
→ `model_level_metrics.csv` (192 rows = 48 models × 4 specifications).

**Type S cannot be summarised by the log-scale model at the meta-analysis level.** The
model needs an offset of 0.025 to handle exact zeros, and at that level almost every
observed value is far below it (median 7×10⁻⁸). The fitted summary then reflects the
offset more than the data. Every affected row carries
`summary_dominated_by_offset = TRUE`, and both the raw median and the model value are
reported.

---

## 3. The sensitivity analyses, and what question each answers

| analysis | question it answers | headline | file |
|---|---|---|---|
| **Optimistic assumed effect** | if we assume the most favourable effect the data support, is power still low? | 28.1%, Type M 2.02 — yes | `assumed_effect_scenarios.csv` |
| **External assumed effects** | against a conventional effect, chosen without looking at our data? | SMD 7.7 / 20.5 / 41.4%; Zr 8.7 / 33.1 / 65.3%; lnRR 12.5 / 27.3 / 46.8% | same |
| **External, at the meta-analysis level** | is the low corrected value imprecision, or corrected means near zero? | the same 48 models reach **89.2 / 97.4 / 98.6%** against a medium effect — so it is corrected means near zero | same, `aggregation = meta_analysis_level` |
| **Bias-robust estimator** (Yang 2024) | does the regression-based correction drive the result? | primary power 13.4%, Type M 3.86, reversals 20 → 6 | `per_meta_analysis_estimates.csv` |
| **Leave-one-cluster-out** | each observation helps estimate the mean used to judge it — does that inflate the result? | 17.35% → **17.45%**; Type M 2.891 → 2.883. No. | `leave_one_cluster_out.csv` |
| **Leave-one-model-out** | does any one of the 48 models carry the summary? | Yang-2023: MA26 +19.5%. Bias-robust: **MA09 +44.7%**, the only value above 20% anywhere | `loo_influence.csv` (576 rows) |
| **Leave-one-paper-out** | does any one of the 28 *papers* carry it? | `k`-weighted MA09 +44.7%; **equal-weighted the largest is a different paper at +8.4%**, median 1.4% — so it is the weighting, not the paper | `leave_one_paper_out.csv` (672 rows) |
| **Alternative estimands** | typical *observation* or typical *meta-analysis*? | 17.4% vs 22.4% (equal) vs 22.3% (random effect) | `primary_level_sensitivity.csv` |
| **rho sensitivity** | does the assumed within-study correlation matter? | ρ ∈ {0, 0.25, 0.5, 0.75} | `rho_sensitivity.csv` |

---

## 4. Every file, and what it is for

### Results you would quote

| file | rows | what it holds |
|---|---|---|
| `supplementary/TableS1_reported_metrics.csv` | 153 | **The supplementary table.** Part A reported results, B primary-level sensitivity, C meta-analysis-level external effects. |
| `supplementary/TableS2_evidence_base.csv` | 28 | **The supplementary table.** One row per source paper. |
| `primary_level_sensitivity.csv` | 36 | primary-study-level summaries, all three estimands |
| `meta_analysis_level_sensitivity.csv` | 39 | meta-analysis-level summaries, both weightings, plus 15 diagnostic rows |
| `assumed_effect_scenarios.csv` | 123 | optimistic and external assumed effects, both levels |
| `per_meta_analysis_estimates.csv` | 48 | one row per model: every estimate, standard error, interval and flag |
| `model_level_metrics.csv` | 192 | power / Type M / Type S per model per specification |
| `reversal_counts.csv` | 5 | sign reversals under each correction approach |
| `evidence_base_characteristics.csv` | 28 | the full-width source behind Table S2 |

### Influence and robustness

| file | rows | what it holds |
|---|---|---|
| `leave_one_cluster_out.csv` | 12 | the self-reference check |
| `loo_influence.csv` | 576 | leave-one-**model**-out |
| `leave_one_paper_out.csv` | 672 | leave-one-**paper**-out, both weightings |
| `rho_sensitivity.csv` | — | the assumed within-study correlation |
| `ma_level_uncertainty.csv` | 60 | three ways of estimating the meta-analysis-level interval |

### Verification

| file | what it holds |
|---|---|
| `verification_audit.csv` | the 11 second-derivation checks, with the tolerance each was held to |
| `yang2024_reference_validation.csv` | our CR2 implementation against the authors' published example |
| `supplementary/metadata_columns.csv` | what every column of Tables S1 and S2 means |
| `supplementary/metadata_files.csv` | the file index for the archived deposit |
| `supplementary/captions.md` | the two supplementary captions |

### Computed but **not** used

| file | why it exists, and why it is not quoted |
|---|---|
| `ma_level_paired_contrasts.csv` | Paired paper-level bootstrap contrasts. Correct and gated, but **deliberately not used as a sensitivity analysis for now** (your decision, 2026-08-15). Kept so it need not be recomputed. Do not quote it without deciding to. |
| the 15 `diagnostic_*` rows in `meta_analysis_level_sensitivity.csv` | alternative critical value, CR1 variant, and the hybrid pairing. Retained so the choices are visible rather than assumed. **Never quote a `diagnostic_*` row.** |

---

## 5. Three numbers that look wrong and are not

**Type M of 26,326.** The primary-study-level maximum under the Yang-2023 correction. Real
and expected: the assumed effect is the denominator, and one corrected mean is very close
to zero. It is a reason to interpret Type M cautiously near zero, not a computational
error.

**Type S of 7 × 10⁻⁸.** The meta-analysis-level median before correction. Also real: a
pooled estimate with a large effect-to-standard-error ratio has an essentially zero chance
of reaching significance in the wrong direction.

**Confidence limits outside a metric's range.** The intervals are back-transformed from a
log-scale model and are not bounded, so a power limit can exceed 1 and a Type S limit can
fall below 0. In Table S1 these are constrained to the bound and marked with an asterisk,
as in the submitted supplement. **No point estimate ever required constraining.**

---

## 6. Two numbers that WERE wrong, now corrected

Recorded so they are not reintroduced.

**"mean: 114%"** in the current Results paragraph. A statistical power cannot exceed 100%.
The figure is `exp(μ̂ + σ̂²/2)` (`S2_v2.R:1391`) — the arithmetic mean of a *lognormal*
variable, which is unbounded while a probability is not. `docs/17` deletes every such
parenthetical.

**The stale clustering figures.** Before the clustering decision of 2026-08-15, the
primary-study-level numbers used the raw `study_ID`, which merges observations from
different meta-analyses that happen to share a label. Those values (0.17154 power,
optimistic 0.27335, and the corresponding external ladder) appeared in
`results/revision/README.md` §4a and `docs/16_handoff.md` §3. **Both are now corrected in
place.** The current values are 0.17354 and 0.28057.

---

## 7. How to check any of this yourself

```bash
Rscript R/revision/run_all.R                        # regenerates everything
Rscript R/revision/20_verify_reported_numbers.R     # 11 second-derivation checks
python3  R/revision/17_verify_scenarios.py          # 123 scenario rows, from scratch
python3  R/revision/21_audit_manuscript_claims.py   # 90 manuscript claims
```

The first is deterministic — a rerun reproduces every file byte for byte. The last is the
one to rerun after **any** edit to `docs/17`, because it is the only check that sees the
manuscript text.
