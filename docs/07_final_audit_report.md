# Final consolidated audit report

> ## ADDENDUM (2026-08-11): Table S1 located — three open items resolved, one arbitration reversed
>
> **Correction to this report and to the panel's record: Table S1 exists.** It is
> in `docs/Supplemental1.pdf`, pages 3–5. Nightingale's claim that it "exists
> nowhere in the repository or either .docx" was wrong, and every downstream
> statement that depended on its absence is withdrawn. Table S1 reports median,
> 95% CI, mean, `k` and `N` for power, Type M and Type S at both levels, for all
> three effect-size metrics, under both `β_0[overall]` and
> `β_0[bias-corrected]`.
>
> ### 1. The 13.4% is resolved, and Turing was right — Nightingale's arbitration is reversed
>
> **Table S1 reports 0.1378, not 0.134.** The pipeline reproduces it exactly:
> under Yang's rule (the *uncorrected* model's variance components) primary
> Type S corrected mean = **0.13775**. Table S1 = **0.1378**. Median 0.09854
> against Table S1's 0.0985; uncorrected row 0.02691 / 0.04339 against
> 0.0269 / 0.0433.
>
> So the number was never unreproducible. **Turing's 0.1378 was exactly correct**,
> and Lovelace's "reproduces under no variant" was true only because it targeted
> the manuscript *text*'s 13.4% rather than Table S1's 0.1378. Nightingale's
> arbitration ("Lovelace right, Turing wrong") is reversed, as is its conclusion
> that the required variance 0.5047 "matches nothing" — the correct comparison was
> never to a fifth variance.
>
> **The defect is a single transcription error in the manuscript text:
> 0.1378 → "13.4%" instead of "13.8%".** Not an analysis problem.
>
> ### 2. All twelve reported means reproduce against Table S1
>
> The panel's "8 of 12 reproduce, 4 do not" was an artifact of comparing against
> the manuscript text. Against Table S1, all twelve reproduce, and they confirm the
> three-formula pattern exactly: MA power corrected 0.62 and primary power
> corrected 0.12 and primary Type M corrected 9.49 use the **uncorrected**
> dispersion; MA Type M corrected 5.84 and MA Type S corrected 0.0370 use their
> **own**; primary Type S corrected 0.1378 uses the **uncorrected** variance
> components. `outputs/reproduction_check.csv`'s 13 checks therefore understate the
> pipeline's fidelity rather than overstating it.
>
> ### 3. The 17.20% and the 18.0% upper limit are resolved — no discrepancy remains
>
> Table S1 reports primary uncorrected power to **two decimals**: median **0.17**,
> CI **[0.16, 0.18]**, mean **0.23**. The pipeline gives 0.17154 [0.16428,
> 0.17912], mean 0.23105 — which rounds to Table S1's values exactly, including
> **0.17912 → 0.18**.
>
> The manuscript text then re-expressed these with **inconsistent precision**:
> "17.20%" is 17.154 rounded to three significant figures (17.2) written with a
> spurious trailing zero; "16.4%" is 16.428 at three significant figures, i.e.
> taken from the underlying value; but "18.0%" is Table S1's already-rounded 0.18
> re-expressed as a percentage, when the underlying value is 17.9%.
>
> So the lme4-version hypothesis is unnecessary (as already established), *and* the
> transcription account is now exact rather than inferred. **The 18.0% is not an
> unexplained discrepancy — it is a 2-dp value quoted to 3 significant figures.**
> `R/03_baseline.R:19-24`'s stale text should be replaced with this, and the
> reproduction gate's `primary power` tolerance can be tightened accordingly.
>
> ### 4. New from Table S1: the 0/1 constraint IS documented, and the CIs do breach the bound
>
> Table S1's caption states: *"For statistical power and type S error, values below
> 0 or above 1 were constrained to 0 and 1 (denoted as 0\* or 1\*)."*
>
> This **partly retracts the D8 criticism**: the flooring is disclosed as a method
> in the supplementary caption, not only in a code comment. The remaining valid
> criticism is narrower — it is not stated in the Methods, and Yang's "negative
> variance" explanation is still the wrong diagnosis.
>
> But it also reveals something the panel could not see, and it is worse than what
> was reported: **the constraint is actually exercised on the CI upper limits for
> power.** Table S1 shows `1*` for MA-level power in four rows (lnRR uncorrected
> and corrected, SMD uncorrected, Zr uncorrected). Recomputed: the MA-level power
> CI upper limit is **1.635 for lnRR** and **1.041 for SMD** — genuinely above 1,
> then constrained. So the bound violation is not confined to the discarded
> lognormal mean; it reaches the **retained** summary's intervals. This strengthens
> the case for a bounded-scale alternative and belongs in the D8 evidence.
>
> ### 5. New from Table S1: `N` is labelled "number of primary studies" and is the cluster count
>
> Table S1 reports `N = 1415` for the All rows, and 224 / 945 / 246 for lnRR / SMD
> / Zr, summing to 1,415 — with the caption defining `N` as "number of primary
> studies". That is the **(study_ID × dataset) cluster count**, which
> double-counts a study appearing in several meta-analyses of the same paper (130
> labels do so, and MA22_01/02/03 share the same 30 studies). Against the sibling
> repository's authoritative 1,187 unique referenced primary studies, Table S1's
> `N` overstates by ~19%. So the "5,740 primary studies" error has a companion in
> the supplement, and fixing the abstract alone will not make the two consistent.
>
> ### 6. Corrections to §5 of this report
>
> - Item 12 ("Table S1 cited five times and exists nowhere") is **withdrawn**.
> - Item 1 ("mean: 13.4%" reproduces under no rule) is **replaced** by: the text
>   mis-transcribes Table S1's 0.1378 as 13.4%; it should read 13.8%.
> - A new item: **Table S1's `N` = 1,415 "number of primary studies"** is a
>   cluster count, not a study count.
> - §7's "Could not check" items on Table S1, the 18.0% limit and the fifth
>   variance are all **withdrawn as resolved**.
>
> ### 7. The scenario-selection rule: coding error, documentation error, or Yang's method?
>
> **Verdict: not a coding error. A documentation error in our manuscript, on top of
> a deliberate and self-acknowledged methodological feature of Yang's procedure.**
>
> **(i) What Yang intends.** Yang's paper defines the rule by **sign, explicitly
> and deliberately**: *"a slope (β1[small−study] or β2[time−lag]) with opposing
> direction (unexpected sign) indicates no detectable publication bias and
> subsequently does not require correction for such a bias"*, and the Fig. 3 legend
> defines all four scenarios purely by "expected sign" / "unexpected sign". They
> keep a *separate*, significance-based notion for reporting prevalence ("a
> small-study effect is statistically detected if…"). Decisively, they acknowledge
> the consequence in their own Results: *"β1[small−study] from 54 (62%)
> meta-analyses were in the expected direction, indicating that these
> meta-analyses exhibited a **(statistically non-significant)** tendency for a
> small-study effect (note that the likelihood of a meta-analysis to show this
> tendency is 50% if there is no real effect)."* And their "Deviations and
> additions" section records that a significance-based two-step procedure **was**
> the Stage 1 registered plan and was abandoned at Stage 2. So sign-only is
> intentional, documented, and its 50%-under-null property was known to Yang.
>
> **(ii) What our code implements.** Exactly Yang's rule.
> `S2_v2.R:529-530, 581-582, 587-588, 642-643` subset on
> `beta0Tbeta1`/`beta0Tbeta2` signs alone, matching Yang `Rmd:932, 1016, 1053,
> 1132` and their SMD and Zr analogues; and `S2_v2.R:429-443` computes the
> `pval < 0.05` flags separately for reporting, matching Yang `Rmd:765-779`.
> `R/01_estimates.R:112-117` reproduces this. **No coding error here.**
>
> **(iii) What our manuscript says.** *"When evidence of publication bias was
> detected, we applied the full bias-correction model; if only one type of bias …
> was detected, we used a reduced version of the model including only the
> **significant** bias term. If no bias was detected, the uncorrected mean was
> used."* (`_rev.docx`, para 41 — unchanged from the original.)
>
> "detected", and especially "the significant bias term", describe a
> **significance-based** rule. The code applies a **sign-based** one. This is a
> documentation error, and it is the reason the census reads as alarming: 25 of 48
> meta-analyses have both slopes non-significant, and 15 of the 34 that gate to
> `beta0_c2` had no *significant* bias in either test — all fully consistent with
> Yang's method and flatly inconsistent with our description of it.
>
> **What this changes.** The finding is downgraded from "the code does not do what
> the manuscript says, and that may be a bug" to "the manuscript misdescribes a
> faithfully implemented method". It does **not** dissolve the substantive issue:
> under a sign-only rule roughly half of bias-free meta-analyses will be
> "corrected" by construction — Yang say so themselves — and that is what drives
> the 20 sign reversals and the near-zero corrected means. So D3 remains open and
> remains the most consequential decision, but it is now a question about adopting
> Yang's method rather than about a defect in our implementation of it. The
> honest options are to correct the Methods wording to describe the sign rule and
> discuss its 50%-under-null property, or to depart from Yang deliberately and say
> so.

Full panel complete: **Lovelace** (code and reproducibility), **Turing**
(statistics), **Wald** (adversarial re-derivation), **Nightingale** (gate).
Nightingale's status: **`revise before Pancho`** — nothing found invalidates the
pipeline, but the record required the corrections applied below.

**Nothing has been applied.** Verified by Nightingale independently, and by mtime
rather than git (`R/` and `outputs/` are untracked, so `git status` cannot detect
an overwrite): `S2_v2.R` md5 `a910b4cc1ac134b9792f2da5d0558ef9` identical in MAIN
and worktree with an empty `git diff HEAD`; **all 29 `outputs/` files have mtimes
15:24–15:39 on 2026-08-10**, i.e. before `docs/04` (18:27), `docs/05` (18:28) and
`docs/06` (21:06) — so the artifacts the reviewers cited are the artifacts that
were gated; `R/04_sensitivity.R:54` still `arrange(desc(shrinkage_pct))`;
`B1_shrinkage_concentration.csv` still `0.16970206383797365`; no commit since
1f14ab1.

Prior reports: Turing `/private/tmp/claude-501/audit_turing2/turing_report.md`;
Wald `/private/tmp/claude-501/audit_wald3/f1.md`; Nightingale
`/private/tmp/claude-501/audit_nightingale/gate.md`. Panel detail in
`06_panel_findings.md`; corrections to that record are listed in §6 below.

---

## 1. Findings confirmed by the full panel

**The `~` defect.** `S2_v2.R:544, :554` omit the `~`, so `metafor` partial-matches
`mod` → `mods` and fits ONE composite moderator (the arithmetic sum of a sampling
variance and a centred publication year) instead of two. Proved structurally
(`mods` is argument 4, `...` argument 28) and empirically
(`rownames(beta)[2] == "mods"`). On `lnRR/MA09.csv`: as-written **0.1059904365**,
intended **0.0681099550**, intercept-only 0.2194470305 — so the moderator was
fitted, not ignored. MA09 is the sole affected dataset (lnRR scenario counts
s1=1, s2=0, s3=3, s4=1); `beta0`, `beta1`, `beta2` and `scenario` are bit-identical
between paths, so it cannot propagate through scenario assignment. Two occurrences
in `S2_v2.R`, **zero in Yang's Rmd** — the defect is ours. All six downstream
deltas reproduce independently.

**The pipeline's arithmetic is correct.** Every reported median reproduces to the
last printed digit under closed-form Type M and independently written fits. The
optimizer table is exactly right including `uses_ess` (all 29 `rma.mv` calls
inventoried). No object-alignment error. No fit silently lost — all `NA` counts
zero, `fit_rma`'s `tryCatch` never fired. **No singular or boundary `lmer` fit
anywhere** in ~30 refits. B5 coverage is *provably* complete: clusters partition
each dataset and per-case counts equal `k` for all 48. Wald's verdict: *"The
problem is the specification, not the code."*

**Sign reversals.** **20 of 48** corrected means oppose `beta0`, covering
**1,932 of 5,740 rows (33.66%)**; all 20 pass the gate; all 20 have a `beta0_c2`
CI including zero; **three have |t_c3| ≥ 1.96** (MA08 2.758, MA28_01 2.215,
MA41_01 2.062; 725 rows). **37 of 48** corrected means have a CI including zero
(3,158 rows, 55.0%). MA08: +0.3359 → −0.1214, `t_c3` = 2.76, scored by the gate as
64% "shrinkage".

**The corrected side is driven by near-zero denominators.** All three MA-level
Type M > 20 are sign flips: MA39_1 **4,259.6**, MA26 77.5, MA13_03 24.5. Under the
corrected effect 23.0% of rows sit at the power = α asymptote and 9.6% at
Type S ≈ 0.5. Flooring |µ| at 0.5·se — withdrawing no correction — cuts primary
Type M 8.123 → 6.017 and eliminates all three >20; three variants bracket it at
**5.06–6.02** against the reported 7.79/8.12.

**B1: pervasive in coverage, mildly concentrated in magnitude.** The top five hold
**2,926 of 5,740 rows = 50.98%** of `k`, so their 58.86% share is only **+7.9 pp**
over that floor; inside/outside shrinkage 0.1789/0.1300 = **1.38×**; and **34 of 48
shrank, covering 87.5% of rows**. Both the original claim ("pervasive") and its
correction ("concentrated") were wrong as stated.

**D4 arithmetic.** `(1|study_ID)` 17.1538% [16.428, 17.912], se(log) 0.02207,
v(study) 0.534 / v(resid) 0.1485. `+(1|case)` 22.2310% [17.282, 28.598],
se(log) 0.12849, v 0.1599 / **0.7659** / 0.1233. SE inflation **5.822×** (so
"7.6-fold" is wrong). Nested ≈ crossed (22.2644 vs 22.2310). `v(case)`
0.7455–0.8477 across four subsets including a collision-free one, so the
reallocation is **genuine structure, not an identifier artifact**. Nightingale
independently reproduced the discriminating fixed-`case` fit: **22.3708%,
se(log) 0.02421, [21.33, 23.46]**.

**B5.** 1,415 clusters, 5,740 rows, 0 skipped. 0.17153780 → 0.17223559
(+0.0698 pp). Per-row −30.057 to +29.496 pp; mean |Δ| 0.8735 pp; >1 pp 21.88%,
>5 pp 3.24%, >10 pp 0.73%. Movement is systematic in dataset size: mean |Δpp|
**9.61 / 4.82 / 1.32 / 0.51** for <5 / 5–10 / 10–25 / ≥25 studies. Tails are
genuine self-inclusion, not numerical instability.

**Type S is not identified without the offset.** `exp(β₀)` = 0.025583 against an
offset of 0.025. The MA point estimate spans **497×** across offsets 1e−6…0.1
(1.885e−06 → 9.369e−04); the lower limit is negative at **every** offset ≥ 0.001;
40 of 48 MA values lie below 10% of the offset and 56.2% of rows below it — and
the offset was introduced for **exactly one exact zero in 5,740**.

**The `k`-weighted `lm` selects for high power mechanically.**
`cor(log k, log se_beta0) = −0.7662`; the top five by `k` carry **62.30%** of the
weight with MA power 0.99993 / 1.0000 / 1.0000 / 1.0000 / 0.9725; 16 of 48 are
≥0.999; Shapiro–Wilk W = 0.89432, **p = 0.000414**. MA-level power spans
**0.569 / 0.710 / 0.822 / 0.886** on summary choice alone.

**All three metrics are one scalar.** `se` cancels in Type M; Spearman(MA power,
|t|) = **1.000000** exactly. So "low power **and** high Type M **and** elevated
Type S" is one finding stated three times.

**Missing uncertainty.** `se_beta0` is discarded at the primary level; **17 of 48**
`beta0` CIs include zero, so for those meta-analyses power is not distinguishable
from α.

**Other confirmed:** the 105 filtered rows are 104 missing `es` (MA41_02/03/04) +
1 missing `year` (MA06_01), all with usable `var`, and dropped rows are marginally
*more* precise — a precision-based upward bias is refuted in direction. Latin-1
in exactly `SMD/MA24.csv`, `SMD/MA26.csv`, `SMD/des_stat/MA18.csv` of all 48. The
seed does not reproduce Type M across scripts. 17.20% vs 17.154% is transcription
(0.17154 × 1.34692 = 0.23105 = the reported 23.1%; 0.1720 gives 0.23167).

### New from Nightingale — the first audit of the main figure

**Figure 2's six "overall mean" crossbars disagree with every value they
illustrate.** The plotted data are **192 values** (4 per meta-analysis × 48), and
the "Sampling"/"cSampling" rows are *within-meta-analysis medians*
(`S2_v2.R:2935, 2941`), **not** the 5,740 rows. `stat_summary(fun = mean)` acts on
the mapped `y`, so for Type M — where `aes(y = log(power))` — the crossbar is an
unweighted **geometric** mean; for power and Type S an **arithmetic** mean. Legacy
spec:

| metric | facet | figure crossbar (unc / cor) | manuscript text | statistic |
|---|---|---|---|---|
| power | primary | 0.3261 / 0.1611 | 0.1720 / 0.0906 | arithmetic |
| power | MA | 0.7100 / 0.4136 | 0.8220 / 0.4490 | arithmetic |
| Type M | primary | 2.609 / 7.007 | 2.86 / 7.79 | geometric |
| Type M | MA | 1.388 / 2.921 | 1.11 / 2.03 | geometric |
| Type S | primary | 0.0534 / 0.1552 | 0.0269 / 0.0985 | arithmetic |
| **Type S** | **MA** | **0.01346 / 0.07007** | **0.00058 / 0.0121** | arithmetic |

The MA Type S line is **23× the reported uncorrected value**; primary power is off
by 1.9×. A defect rather than a design choice, because the legend says "the
horizontal line denotes the overall mean" and the text labels its own
parentheticals "mean" — one word, two quantities.

**The "20+" truncation does not exist.** `S2_v2.R:3013-3021` is
`scale_fill_gradient(limits = c(0, 20))` with **no `oob` argument**, so ggplot2's
default `censor` sets out-of-range to `NA` and `geom_tile` draws
`na.value = "grey50"`. **10 of 192 tiles render grey** (14890.4, 4260.9, 615.6,
77.7, 66.2, 42.4, 42.0, 24.4, 22.7, 22.5) and the colourbar carries no "20+"
label — in a white→teal ramp they read as *missing data*. Wald was right about the
heatmap and wrong about the violin: `limits = c(-1, 10)` drops nothing, so 4,259.6
is plotted at y = 8.36, unlabelled.

**Scenario assignment is by sign with no significance test** (`R/01_estimates.R:112-117`
= `S2_v2.R`). Census: `beta1` non-significant in **27 of 48** (1,433 rows, 25.0%);
`beta2` in **42 of 48** (3,527 rows, 61.4%); **both non-significant — scenario set
by noise — in 25 of 48 (1,293 rows, 22.5%)**; at least one in 44 of 48 (63.9%).
And **of the 34 that gate to `beta0_c2`, 15 (688 rows) had no detected bias in
either test.** *The manuscript describes detect-then-correct; the code corrects
regardless of detection.* Demonstrated fragility: MA41_04's `beta1` is **+0.0099**,
and imputing its 57 missing `es` flips the scenario 3 → 4 under three of four
imputations and to 1 under the fourth.

**MA41's missing data, settled as far as the files allow.** No magnitude-bearing
column survives (`obs_ID, study_ID, year, es, var, study_ID_org, meta_ID` only).
Missingness is **study-level, not effect-level**: 18 studies touched, **13 wholly
missing**, 5 partial — consistent with unrecoverable primary reports rather than
effect-dependent censoring. `year` null (p = 0.2475); within-study paired implied-n
null (p = 0.50). **The two MA41 sign flips are robust to the missing rows.** So
the threat the panel feared is smaller than feared, while the threat to the
*scenario mechanism* is larger than anyone had stated.

**A revised manuscript exists.** `docs/AnimalCogn_statistical power_BiologyOpen_rev.docx`
— revision 16, `lastModifiedBy` Ayumi Mizuno, modified 2026-08-11T00:22Z,
TotalTime 92 min. **Hand-edited by Ayumi, not an agent, so not a constraint
breach.** Diff: 4 changed paragraph blocks, **zero numeric changes**; the abstract
already reads "5,740 effect-size observations from primary studies", so one audit
correction is partly applied. All twelve means, all thirteen gate targets and the
Figure 2 legend are identical in both files. **`docs/03_manuscript_corrections.md`
targets the original — re-diff before applying anything.**

---

## 2. Findings where reviewers disagreed, and how Nightingale arbitrated

**(1) The reported 13.4% — Lovelace right, Turing wrong; and the escape route is
closed.** No single dispersion rule reproduces all twelve manuscript means.
Switching to the uncorrected `var(log)` repairs three of the four failures to
≤0.06% but breaks MA Type M corrected (−65.6%) and MA Type S corrected (−52.8%).
13.4% misses by **2.80%** at best — 30–100× the reproduction error of its own
comparators and 7× the ±0.4% rounding tolerance of a 3-s.f. figure. Lovelace's
required variance **0.5047** is exactly right; the four available variances are
0.5513, 0.5631, 0.8835, 0.9024. **Wald's git test was run and it refutes the
"different data states" hypothesis**: the data directories were committed once
(2025-07-18) and never modified, `S2_v2.R`'s last commit is 2025-08-21 with only
cosmetic changes since 2025-08-07, and the .docx was created 2026-07-27 — one
static script, one static dataset. Remaining explanations: a typing error, or a
fifth ad hoc variance in an interactive session. **Report as unreproduced;
recompute. Do not present 0.1378 as agreement.**

**(2) B1 — report the counterfactual, not a share.** Every share invites the
k-floor objection: `shrinkage·k` 58.86% against a 50.98% floor (+7.9 pp);
`|pull|` 67.97% against a **45.71%** floor (+22.3 pp — better, not floor-free).
Use `|pull|` to *select* the five (its decomposition is exact:
`sum(pull) = −0.605717 = log(0.4486/0.8221)`) and report the re-aggregation:
MA level **0.8221 → 0.4486 (ratio 0.546) becomes 0.7511 → 0.5254 (ratio 0.700)**.
Two label corrections: Turing's 73.95% / 57.7% / −0.74472 are **corrected-spec**
quantities; on the manuscript's legacy spec the same ranking gives **67.97% /
48.05%**. And Turing's primary-level "+9.8%" **does not reproduce** (+1.8% under
one reading, +6.3% under another).

**(3) D4 — compatible, not competing; stays open.** Turing's decomposition is
correct and Nightingale reproduced the disputed fixed-`case` row exactly. Wald's
"the width is the real question" is a priority judgement, not a rival finding. One
nuance: the two decisions are *less* separable than stated, because the
equal-weight estimand's interval is narrow only if you also declare the 48 to be
the population — 22.371% [21.33, 23.46] against a sample-based 22.00%
[16.92, 28.61].

**(4) B5 — "negligible" is not defensible.** Both reviewers right about different
halves. The defensible sentence: *"Excluding the focal study from the assumed
effect changes the pooled primary-study power estimate by +0.07 pp, but this
aggregate is a near-cancellation (61.6% of rows move up, 38.2% down; mean |change|
0.87 pp, range −30.1 to +29.5 pp). The movement is systematic in the number of
studies: mean |change| is 9.61, 4.82, 1.32 and 0.51 pp for meta-analyses with <5,
5–10, 10–25 and ≥25 studies. 21.9% of rows move more than 1 pp and 3.2% more than
5 pp."* Do not attach a number to the coarser-identifier bound.

**(5) Sign flips — two phenomena, on the remedy test; stays open.** A 0.5·se floor
fixes the 17 low-|t| cases and removes all three MA-level Type M > 20 while
withdrawing no correction, but does nothing for MA08, MA28_01 and MA41_01, whose
corrected estimate is a *significant reversed* effect — a substantive claim, not a
denominator artefact. Different remedies ⇒ two phenomena. New third remedy:
**9 of the 20 flips (416 rows) had no detected bias in either test.**

---

## 3. Coding bugs that must be fixed (ranked by Nightingale)

1. **`S2_v2.R:3013-3021`** — `scale_fill_gradient(limits = c(0,20))` with no
   `oob`: 10 tiles censor to grey and read as missing data. Published-figure
   defect. Either `oob = scales::squish` with a "20+" break, or change the legend.
2. **Figure 2's crossbars use a different statistic from the text** (six lines,
   one wrong by 23×). Recompute to the reported statistic, or name the figure's
   statistic in the legend.
3. **`S2_v2.R:3010`** — a stray literal `054042` on its own line, introduced in
   ac1e5cd. Harmless at runtime, but this is the supplementary code reviewers read.
4. **`S2_v2.R:544, :554`** — the missing `~`.
5. **`R/03_baseline.R:19-24`** — still asserts the lme4-version hypothesis that
   §2l retired as transcription. Stale text in a file likely to be quoted.
6. **`SMD/MA26.csv`** byte 0x88 is not printable latin-1, so
   `locale(encoding = "latin1")` would yield a control character; likely
   windows-1252 or mac-roman. Determine the encoding per file.
7. **`outputs/legacy_mean_diagnosis.csv`** is an orphan —
   `03_baseline.R:120` writes `summary_statistic_comparison_power.csv` from the
   same block.

Also outstanding from the earlier panel, not re-ranked: the B1 sort key
(`R/04_sensitivity.R:54`); the seed's failure across scripts; the two tautological
`stopifnot` guards; the `confint` default (`profile` in `S2_v2.R` and Yang, `Wald`
in `R/`); the wrong `zr` failure-mode rationale in `R/00_setup.R:9-11` and
`R/README.md`.

---

## 4. Methodological choices for Ayumi and Shinichi — all four still OPEN

**D3 — the one-directional gate. Now the most consequential of the four.** The
gate is sign-blind, and that choice drives the corrected headlines and four
Discussion sentences. Evidence: 20 of 48 flips / 33.66% of rows, three of them
significant reversals; 37 of 48 corrected CIs spanning zero; and Nightingale's new
census — scenario assignment uses **sign only, no significance test**, so 25 of 48
meta-analyses (22.5% of rows) have their scenario set by noise, and 15 of the 34
that gate to `beta0_c2` had **no detected bias in either test**. Three costed
remedies: floor |µ| at 0.5·se (Type M 8.12 → 6.02, no correction withdrawn);
drop the 9 with |t_c3| < 0.5 (→ 5.06); sign-preserving gate (→ 3.71, but
confounded — withdraws correction from 20 of 48). A fourth, uncosted: do not
correct where no bias was detected.

**D4 — two questions.** Which unit is weighted equally (study clusters 17.15% vs
meta-analyses 22.37%), and whether the 48 are a sample or the population
([21.33, 23.46] vs [17.28, 28.60]). Settled by the panel and no longer open: 48
levels support the term, no singular fits, the variance reallocation is genuine,
crossed ≈ nested.

**D7** — 37 of 48 `beta0_c2` CIs include zero, so a farther-from-zero limit built
on `beta0_c2` would usually bound an estimate consistent with no effect, with an
arbitrary sign for the 20 flipped cases. `se_beta0_c2/se_beta0` median 1.844,
max 7.043.

**D8** — the estimate is not identified without the offset (497× span), and the
negative lower limit is *slightly worse* under the sign-preserving gate, so it is
not an artifact of the current specification. Options: report offset sensitivity;
raw median + IQR at MA level; logit scale (`qlogis` gives 0.00148 [0.00065,
0.00335]); threshold counts (1,601 of 5,740 rows exceed 0.05); or keep the offset
and state the flooring as a method. Do not repeat Yang's "negative variance"
explanation — it is wrong.

**Also for decision, not previously framed as such:** whether the reported
summary should remain a `k`-weighted geometric mean, given that it spans 57–89% on
weighting/transformation choice and puts 62% of the weight on ceiling-pinned
units; and whether `se_beta0` should propagate into primary-level uncertainty.

---

## 5. Manuscript claims that no longer appear defensible

1. **"mean: 13.4%"** (primary Type S corrected). Reproduces under no rule.
2. **"Type M error exceeded 20 in three of the included meta-analytic models after
   bias correction… may substantially overstate the assumed underlying effect
   sizes."** All three are sign flips with corrected CIs spanning more than the
   uncorrected effect; the values measure 1/noise. Zero under the floor variant.
3. **Primary-level "some values exceeding 20."** 796 rows, of which **520 (65%)
   lie in the 20 flipped meta-analyses**. 276 remain — survives partially.
4. **Type S "9.85%… exceeding typical estimates in ecology and evolution
   (5–8%)."** The comparison **reverses** at 4.27% under the sign gate and sits
   inside or at the edge of the band under other variants.
5. **Type M "7.79 vs eco-evo 2.5 to 4.0."** Lands *inside* the comparator range at
   3.71; above but far less dramatic at 5.06–6.02.
6. **Figure 2's legend** — "values > 20 shown as '20+'" describes behaviour the
   code does not implement, and "the horizontal line denotes the overall mean"
   describes a statistic that is not the one reported.
7. **"5,740 primary studies"** — effect-size estimates. Already partly fixed by
   hand in `_rev.docx`.
8. **"average sampling variance across all effect sizes"** — the code uses
   `se_beta0`.
9. **"average" throughout** — the quantity is a weighted geometric mean, off from
   the raw median by 26% (primary power) to four orders of magnitude (MA Type S),
   in both directions.
10. **Methods line 45's "linear mixed-effects models"** at the MA level — the code
    uses weighted `lm`; line 46 is correct.
11. **The detect-then-correct description** — the code assigns scenarios by sign
    alone and corrects 15 meta-analyses (688 rows) where no bias was detected.
12. **"Table S1" is cited five times and exists nowhere** in the repository or
    either .docx. Until it does, five citations point at nothing and the
    17.20%/18.0% question cannot be closed.

---

## 6. Corrections applied to `06_panel_findings.md`

Per Nightingale's `revise before Pancho`:

- **§2h contradicts itself** — "all four reproduce exactly" under the uncorrected
  rule is **false**; that rule repairs three and breaks two. Corrected in §2 above.
- **§2e's 73.95% / 57.7% / −0.74472 are corrected-spec**, presented in a
  manuscript-reproduction context; the legacy-spec values are **67.97% / 48.05%**.
- **§2e's primary-level "+9.8%" is unreproduced** (+1.8% / +6.3%).
- **The 0.1697 label "uncorrected" is wrong** — it is the coded-`desc(shrinkage_pct)`-sort
  value, nothing to do with publication-bias correction.
- **§2j survives with a caveat**: `04_sensitivity.R`'s 1.10738912 is a genuine
  second realisation, but it coincides with `03`'s *legacy* value to 13 dp only
  because `error_M` rounds to 3 dp and 44 of 48 Type M values sit near 1.00 — which
  is itself evidence of how coarse the grain is.
- **Figure 2 is built from 192 values, not 5,740** — nowhere stated in the record.
- **Wald's Figure 2b claim must be split**: the heatmap conceals (grey NA), the
  violin plots 4,259.6 unlabelled.

---

## 7. Could not check

- **Table S1** — does not exist in the repository or either .docx; `docs/` is
  untracked in MAIN so there is no git history. The 17.20% and the **18.0% upper
  limit** cannot be closed without it or the original session objects.
- **Whether MA41's missing `es` are missing at random with respect to effect
  magnitude** — unanswerable from these files (no magnitude-bearing column
  survives). Consequences are now bounded instead; settling it needs the MA41
  source paper.
- **Which fifth variance produced 13.4%** — would need the original R session.
- **The rendered PDFs** — the grey-NA behaviour was established from ggplot2
  semantics and the scale call, not from pixels. Worth one visual check.
- **Corrected leave-one-out** — not implemented. Wald's cheap probe: the 12
  datasets with `n_study_id < 10` (~120 clusters × 4 fits), reporting how often the
  gate or scenario flips when one study is removed.
