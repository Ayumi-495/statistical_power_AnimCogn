# Revision handoff — current state (2026-08-15)

**Supersedes `docs/10_revision_handoff.md` entirely.** That file describes the state
before the Yang-2024 work, the clustering decision and the estimand decisions, and its
Section 2 numbers are stale. Read this one.

> **Durability warning.** `.gitignore` contains `/docs`, so `docs/09`–`docs/16` are
> untracked and exist only in this worktree. `docs/02`–`docs/08` were force-added and are
> committed. Force-adding the rest needs Ayumi's approval.

---

## 1. Repository state

| | |
|---|---|
| branch | `claude/ayumi-methodological-audit-3061e3` |
| HEAD | `a503967` "Add leave-one-source-paper-out to the revision workflow" |
| pushed | **no** — commits from `0ec8f5b` onwards are local only |
| submitted analysis | **frozen and verified unchanged.** `S2_v2.R` md5 `a910b4cc1ac134b9792f2da5d0558ef9`; `outputs/` and `figure/` all 32 files unchanged, checked after every run |

**Workflow.** `R/revision/00`–`15` plus `run_all.R`, which runs 01–05, 07, 08, 09, 11,
12, 14, 15, 13 in that order. Scripts 06 (external validation, needs network) and 10
(evidence-base table, needs the sibling systematic map) are run separately.

Results in `results/revision/`: 14 CSVs, `supplementary/` (Tables S1 and S2, captions,
column metadata, file index) and `figures/` (`main_metrics`, `model_level_metrics`).

---

## 2. Decisions taken, and by whom

All are implemented. Do not reopen without the PI.

| # | decision | date |
|---|---|---|
| 1 | Yang et al. (2023) remains the primary analysis; Yang et al. (2024) FE + VCV is the reported sensitivity analysis; UWLS supplementary | 12 Aug |
| 2 | At the meta-analysis level the bias-robust estimate is paired with **its own** CR2 standard error, not with `se_beta0` | 12 Aug |
| 3 | CRVE is **CR2 with Satterthwaite df**, verified from Yang's tutorial, their 448-model script and metafor's source; validated against their published worked example (9 of 9 values) | 12 Aug |
| 4 | Near-zero assumed effects handled by interpretation, not by any floor, exclusion or sign rule | 13 Aug |
| 5 | Critical value is **z = 1.96** everywhere, matching Yang's own code; stated explicitly in Methods | 13 Aug |
| 6 | The `k`-weighted meta-analysis-level aggregate is a **secondary descriptive summary**; the paper centres on the primary-study level; MA09 is not excluded, its influence is reported | 13 Aug |
| 7 | **Type S** reported as raw median with quartiles, model-based estimate retained for comparability, difference explained | 15 Aug |
| 8 | **Weighting**: effect-size-weighted estimate is the main result with its estimand stated; equal-weight and random-effect model-level summaries are substantive sensitivity analyses, not competing estimates | 15 Aug |
| 9 | **Clustering** prefixed by meta-analysis (`MA_model × study_ID`) as the primary definition; duplicate primary publications across meta-analyses noted as a limitation, not pursued | 15 Aug |
| 10 | Out-of-range confidence limits constrained to the bound and marked with an asterisk, as in the submitted supplement | 15 Aug |

---

## 3. The numbers, as they now stand

Source: `results/revision/supplementary/TableS1_reported_metrics.csv` (153 rows).
Everything includes the `~` correction, closed-form Type M, and the adopted clustering.

**Primary-study level**, each effect-size estimate weighted equally (the reported estimand):

| assumed effect | power | Type M | Type S (raw median [IQR]) |
|---|---|---|---|
| uncorrected | 0.1735 [0.1664, 0.1810] | 2.891 | 0.0183 [0.0020, 0.0573] |
| **Yang 2023 corrected** | **0.0899 [0.0869, 0.0930]** | **7.879** | 0.1358 [0.0402, 0.2580] |
| Yang 2024 bias-robust | 0.1339 [0.1284, 0.1396] | 3.860 | 0.0524 [0.0071, 0.1292] |

Other estimands, uncorrected / corrected power: equal per meta-analysis 0.2237 / 0.1054;
meta-analysis as random effect 0.2226 / 0.1052.

**Meta-analysis level**, `k`-weighted:

| assumed effect | power | Type M | Type S model / raw median |
|---|---|---|---|
| uncorrected | 0.8221 | 1.109 | 0.00058 / 7.2e-08 |
| **Yang 2023 corrected** | **0.3904** | **2.177** | 0.0122 / **0.0022** |
| Yang 2024 bias-robust | 0.4789 | 1.526 | 0.0056 / 7.6e-07 |

**Sensitivity to the assumed effect**, primary-study level:

- optimistic (confidence limit farther from zero): power **0.2806**, Type M 2.0166
- external, SMD d = 0.2/0.5/0.8: 0.0767 / 0.2053 / 0.4144
- external, Zr r = 0.1/0.3/0.5: 0.0865 / 0.3312 / 0.6531
- external, lnRR 10/25/50%: 0.1254 / 0.2726 / 0.4677

> **Corrected 2026-08-15 (second pass).** This block previously carried 0.2734,
> 0.0759/0.2017/0.4082, 0.0860/0.3287/0.6501 and 0.1151/0.2514/0.4437. Those were
> computed with the *raw* study identifier and predate the clustering decision, so they
> were inconsistent with Part A of the same section. `09_assumed_effect_scenarios.R`
> prints the pair on every run: `power, optimistic  raw 0.27335 | prefixed 0.28057`.
> Table S1 was always right; only this summary was stale. **Do not quote the old set.**

**External assumed effects at the meta-analysis level** (Table S1 part C), the answer to
the self-reference objection: against a conventionally medium effect the same 48 models
have power **0.892** (SMD), **0.974** (Zr), **0.986** (lnRR), against a reported corrected
summary of 0.390. So the low corrected value reflects corrected means near zero, not
imprecise meta-analyses.

**Leave-one-cluster-out** (adopted clustering): power 0.1735 → **0.1745**; equal per
meta-analysis 0.2237 → 0.2279; Type M 2.891 → 2.883. **The older 17.15% → 17.22% pair
must not be quoted**; it used the raw identifier for the second-stage summary.

**Sign reversals**: 33/48 (Eq. 2), 26/48 (Eq. 3), **20/48 (as reported)**, 6/48 (FE + VCV),
5/48 (UWLS).

**Verification standard.** Every headline has two derivations: `lme4` against `nlme`
(agreement ≤ 1e-8 across the nine primary-level estimand × assumed-effect cells), and an
independent re-implementation in Python for the meta-analysis-level summaries and the
Table S1 transcription (max difference 4.9e-05, which is display rounding).

---

## 4. OPEN — must be settled before the manuscript text is final

> **Update, 2026-08-15 second pass.** Ayumi revised the Methods in the `.docx` after this
> section was written. **4.2 is now done** (¶0058 discloses the log scale). **4.1 is done
> in the Methods and still open in the Results and Discussion.** **4.3 is untouched.**
> **4.5 (scenario rows were single-derivation) is closed** — `16_export_scenario_inputs.R`
> plus `17_verify_scenarios.py` re-derive all 123 rows in Python, agreement < 1e-6.
> Replacement text for every remaining item is in **`docs/17_results_and_figure_text.md`**,
> which supersedes the RESULTS / ABSTRACT / DISCUSSION / METHODS sections of `docs/15`.

These emerged on 15 August and are **not** resolved.

### 4.1 Type S in the manuscript text contradicts the figure and Table S1 (blocking)

Decision 7 is implemented in `11_main_figure.R` and Table S1 but **not in the manuscript**.
`rev.txt:63` still reports meta-analysis-level Type S as 0.06% and 1.21%, which are the
offset-dominated model values; the raw medians are 7.2e-08 and 0.22%. Figure and text
would report the same quantity differently, by a factor of about 8,000 for the
uncorrected value. This changes numbers, not wording, and should be fixed first.

### 4.2 The Methods do not disclose the log transformation

`rev.txt:59` says only "weighted linear regression (`lm` function in base R), with the
number of effect sizes per meta-analysis as weights". It never says the values are
log-transformed. A reader implementing the Methods as written gets a different number.
Yang et al. (2023) have the identical omission and it was inherited.

### 4.3 The summary statistic is mislabelled

`exp(intercept)` of a weighted log-scale regression is a **k-weighted geometric mean**. It
estimates a median only if the log-scale values are symmetric, and they are not
(k-weighted skewness of `log(power)` = −2.69 at the meta-analysis level, +0.71 at the
primary level — **the sign flips between levels, so no single label is correct for both**).

Verified magnitudes for uncorrected meta-analysis-level power: reported k-weighted
geometric mean **0.8221**; like-for-like k-weighted median **0.9999**; unweighted median
0.9129; unweighted geometric mean 0.5693. The reported value sits at the **17.1st
k-weighted percentile** of the 48 values it summarises, because 19 models have power
≥ 0.99 and carry 66.5% of the weight.

The manuscript body says **"average"** (`rev.txt:63-64`), which is wrong under every
reading: the arithmetic mean is 0.710 unweighted and 0.886 k-weighted. Only the Figure 3
caption says "Median power".

Proposed fix, not yet applied: keep the estimator (decision 8 requires it for
comparability), correct the label level-specifically — "weighted geometric mean" at the
meta-analysis level, "back-transformed model intercept" at the primary level — and
disclose the log transformation.

### 4.4 The confidence interval rests on an unverified weight assumption

`confint()` on the weighted fit treats `k` as an inverse-variance weight. That fails:
weighted residual variance by `k` tertile is 20.2 / 15.8 / 43.4. A nonparametric
bootstrap over the 48 models (B = 20,000) gives [0.652, 0.914] against the reported
[0.717, 0.942], about 28% narrower on the log scale. **Decision needed**: report as is
with a caveat, or report the bootstrap interval, or note it as a limitation. Single
derivation so far.

### 4.5 The leave-one-MODEL-out analysis is `k`-weighted only

`07_influence_loo.R` computes only the `k`-weighted summary. The equal-weight summary is
now a reported sensitivity analysis, so 07 should be extended to both weightings, as
`15_leave_one_paper_out.R` already is. Small job, not blocking.

---

## 5. Corrections that must not be resurrected

Beyond the 16 recorded in `docs/10` section 3, all still valid, this session added:

1. **"R2C13 (conventional effect sizes) was never done."** False. It was done in the
   audit phase (`outputs/B4_*`) and had merely been left out of the revision workflow.
2. **"The optimistic scenario is near-vacuous at the meta-analysis level because power
   → 1."** False. It is an exact additive shift of the Wald statistic by
   `qt(0.975, ddf)`, so the arithmetic floor is 0.50–0.89 and the observed minimum is
   **0.652**. It is uninformative because it is *deterministic*, not because it saturates.
3. **"Scenario assignment uses p < 0.05."** The *correction model selection*
   (`S2_v2.R:529, 581, 587, 642`) uses **sign only**, identically to Yang's own script.
   A separate descriptive detection (`:429, :436`) uses p < 0.05 and sign; it is printed,
   not used to select the correction. Two different objects; do not conflate them.
4. **"11 tiles were censored in the submitted figure."** Recomputation gives **10**, and
   the exact number is not reproducible anyway because Type M was an unseeded Monte Carlo
   estimate. Quote no count.
5. **"The manuscript calls the summary a median."** It says "average". Only the Figure 3
   caption says median.
6. **"The reply to R1C13 fails to name the geometric mean."** It does name it: "equivalent
   to a weighted geometric mean on the original scale". The defect is narrower — it
   asserts "model-based median", which the data do not support, and gives no raw median.
7. **The 22.37% equal-weight figure** is `lmer(log(x) ~ 0 + case + (1|cluster))` with the
   unweighted mean of the 48 case effects on the log scale. A plain mean of the 48
   per-model means gives **22.0%**, a different number. The procedure must be stated or
   it cannot be reproduced.

---

## 6. Leave-one-out: what exists, verified 2026-08-15

**Leave-one-MODEL-out** — `R/revision/07_influence_loo.R`, object `loo`, file
`results/revision/loo_influence.csv`, **576 rows**.

- all **48** models, each dropped in turn — completeness asserted, one row per
  (model × specification × metric)
- all **3** metrics
- **4** specifications: uncorrected, Yang-2023 corrected, FE + VCV, UWLS
- **`k_effect_sizes` weighting only.** `wgeo(v, o$k, off)` at lines 42–50. The
  equal-weight summary is **not** covered
- meta-analysis level only

Largest influences on power: Yang-2023 summary — MA26 +19.5%, MA31 −12.9%, MA08 −6.9%,
median |change| 0.6%. FE + VCV summary — MA09 +44.7%, MA31 −10.2%, MA08 −7.2%. MA09 is the
only model anywhere above 20%.

**The unit is the model, not the paper.** `dropped_MA_model` holds values like
`MA01_01.csv`, `MA01_02.csv`. There are 48 models from 28 papers, and **12 papers
contribute more than one model**, so dropping a model does not remove a paper.

**Leave-one-PAPER-out** — `R/revision/15_leave_one_paper_out.R`, object `lopo`, file
`results/revision/leave_one_paper_out.csv`, **672 rows** = 28 papers × 4 specifications
× 3 metrics × 2 weightings. Built 15 August to the current pipeline: current estimates
and model definitions, closed-form Type M, Type S through `aggregate_ma()` so the offset
flag and raw quartiles travel with every row, z = 1.96, and the aggregate recomputed from
only the models remaining after the paper is dropped. Gates: every combination appears
exactly once; dropped plus remaining always 48 models and 5,740 effect sizes; the 24
baselines match the canonical table to 2.2e-16. 44–47 models remain after a drop.

**A finding worth carrying into the text.** Under `k`-weighting the largest single-paper
influence on the bias-robust summary is MA09 at **+44.7%**; under equal weighting the
largest is a different paper (MA39) at **+8.4%**, with a median absolute change of 1.4%.
So MA09's leverage is a property of weighting by effect-size count, not of that paper
being unusual. This strengthens decision 6.

**The audit-phase `outputs/B2_leave_one_paper_out.csv` is superseded and must not be
quoted.** It predates the adopted clustering, used a Monte Carlo Type M and covered only
two specifications.

**Still `k`-weighted only: `07_influence_loo.R`** (leave-one-model-out). The equal-weight
summary is now a reported sensitivity analysis, so extending 07 to both weightings, as 15
already does, is a small outstanding job.

---

## 7. Deliverables ready to use

- `docs/15_manuscript_edits.md` — every Methods, Results, Abstract and Discussion change,
  as "current text → replace with". Needs updating for sections 4.1–4.3 above.
- `docs/13_crossreference_fixes.md` — eight `Figure 2a/b/c` callouts that should be
  `Figure 3`, `Figure S1` that should be `Figure 1`, the Table 1 collision and the missing
  Table 2 caption.
- `docs/14_methods_and_reply_drafts.md` — the sensitivity-analyses Methods subsection and
  the replies to R1C11 and R1C15.
- `docs/12_pi_question_draft.md` — the three questions, now answered.
- Tables S1 and S2, captions, column metadata, file index — in
  `results/revision/supplementary/`. Table S1 also lives in a Google Sheet tab that Ayumi
  maintains; the CSV is the source of truth.

---

## 8. Not done

- **Response to reviewers**: R1C12, R2C8, R2C9, R2C13, R2C15 and the entire Biology Open
  rubric section are unwritten. All the numbers they need now exist.
- **Zenodo deposit, licence, and the GitHub URL mismatch** (`Ayumi-495` in the manuscript
  against `Ayumi495` in the reviewer's comment) — Ayumi's side.
- **The manuscript itself has not been edited.** Everything is paste-ready text in
  `docs/15`; the `.docx` is untouched.
- **Nothing has been pushed**, and nothing has been posted to issue #1 since Ayumi's own
  message of 15 August.

---

## 9. The AYUMI vault is current as of 2026-08-15

- `01_Projects/statistical_power_AnimCogn/worklog.md` — new entry covering 12–15 August,
  including a section on three of my own defects.
- `01_Projects/statistical_power_AnimCogn/project hub.md` — synthesis, bottleneck,
  durable decisions and next actions all rewritten. **Nothing is waiting on the PI.**
- `03_Resources/Wiki/Methods/repeated mistakes and safeguards.md` — two safeguards added
  under Active: "Apply a decision to every consumer of a shared definition, or to none",
  and "Make display mappings fail closed, and gate artefact-against-artefact consistency".

---

## 10. How to resume after compaction

1. **Read this file first.** It supersedes `docs/10_revision_handoff.md`.
2. **Do not re-audit.** Four reviewers and Ayumi have been through this pipeline. The
   outstanding items are decisions and prose, not unknowns.
3. **Check nothing drifted**: `git log --oneline -5`, `git status`, and
   `md5 -q S2_v2.R` (must be `a910b4cc1ac134b9792f2da5d0558ef9`).
4. **The three things to settle before manuscript text is final** are in section 4:
   Type S in the text, the undisclosed log transformation, and the summary-statistic
   label. Sections 4.4 and 4.5 are smaller and not blocking.
5. **Then write the six outstanding replies** (section 8) using the numbers in section 3.
6. Standing rules that still apply: `S2_v2.R`, `outputs/` and `figure/` are frozen; every
   headline needs two derivations by different routes; nothing is posted or pushed
   without Ayumi saying so.
