DRAFT comment for https://github.com/Ayumi-495/statistical_power_AnimCogn/issues/1
Not posted. Review and edit before sending.

---

@itchyshin - the audit is now finished, so here is the promised update. The short
version: the implementation is a faithful port, the uncorrected results are solid,
and almost everything I raised earlier turned out to be documentation rather than
analysis. But the **bias-corrected** side has a real problem, and it is the one
place where I need your judgement rather than mine — partly because you are an
author on Yang et al., so you can tell me what the procedure was actually intended
to do.

## Resolved since my last comment — no input needed

Table S1 was in the supplementary PDF all along (I had missed it), and it settles
three of the loose ends:

* the "mean: 13.4%" for corrected primary-study Type S is a **transcription
  error**. Table S1 reports 0.1378, and the code reproduces it exactly (0.13775).
  The text should read 13.8%;
* the "17.20% (95% CI 16.4%–18.0%)" is consistent with the code once Table S1's
  2-decimal rounding is taken into account (code: 0.17154 [0.16428, 0.17912];
  Table S1: 0.17 [0.16, 0.18]). The text just mixes precisions — "18.0%" is
  Table S1's rounded 0.18 quoted to three significant figures when the underlying
  value is 17.9%;
* **all twelve** reported "mean" values reproduce against Table S1. My earlier
  worry that several were unreproducible was an artifact of checking them against
  the main text rather than the supplement.

Two things I will simply fix, no decision needed: the missing `~` in the MA09
moderator specification (uncorrected results unchanged; corrected
meta-analysis-level power moves 0.449 → 0.390), and the reproducibility items
(seed, counts, the "average sampling variance" wording, the effect-size vs
primary-study count — note Table S1's `N` = 1,415 is also labelled "number of
primary studies" and has the same problem as the abstract).

I also found two defects in Figure 2 that I will fix: the colour scale has no
`oob` argument, so ten cells above 20 are being censored to `NA` and rendered grey
rather than as "20+"; and the "overall mean" crossbars are computed on a different
statistic from the values quoted in the text (the meta-analysis-level Type S line
is about 23× the reported figure).

## 1. The bias-corrected analysis — the one I really need your view on

This is where the audit changed my mind, and it bears directly on Reviewer 1's
question.

The scenario-selection rule uses the **sign** of the small-study and decline
slopes, with no significance test. Our code does this correctly — it is what the
Yang et al. code does, and your paper states it explicitly ("a slope with opposing
direction (unexpected sign) indicates no detectable publication bias and
subsequently does not require correction"). Yang et al. also note the consequence:
62% of meta-analyses had β1 in the expected direction, described there as a
"(statistically non-significant) tendency", with the observation that the
probability is 50% if there is no real effect. So this is a deliberate feature of
the method, not a bug in our port, and I was wrong to describe it as an
implementation problem in my earlier comment. **Our Methods, however, describe a
significance-based rule ("including only the significant bias term"), which is not
what either code does — so at minimum the Methods wording has to change.**

What concerns me is how that plays out in this corpus:

* **20 of 48** corrected means have the **opposite sign** to the uncorrected mean,
  covering 1,932 of 5,740 effect sizes (34%). All 20 pass the magnitude gate,
  because the gate compares `|beta0|` with `|beta0_c2|` and is blind to sign;
* **37 of 48** corrected means have a confidence interval that includes zero (55%
  of effect sizes);
* because Type M diverges as the assumed effect approaches zero, this drives the
  headline corrected numbers. The three meta-analyses we report as "Type M > 20"
  are **all** sign flips, with values of 4,260, 78 and 25 — the first from a
  corrected mean of 0.00006 whose own CI spans ±0.46. That is measuring 1/noise,
  not effect-size exaggeration;
* 25 of 48 meta-analyses have **both** slopes non-significant, and 15 of the 34
  that the gate sends to the corrected estimate had no significant bias in either
  test.

I looked at what it would take to make the corrected results interpretable
without abandoning the method:

| | primary Type M | meta-analysis Type M > 20 |
|---|---|---|
| as currently reported | 8.12 | 3 |
| floor \|mu\| at 0.5·SE (no correction withdrawn) | 6.02 | 0 |
| exclude the 9 meta-analyses with \|t\| < 0.5 | 5.06 | – |
| sign-preserving gate (confounded: withdraws correction from 20 of 48) | 3.71 | 0 |

So the corrected Type M is somewhere around 5–6 rather than 7.8–8.1 once the
near-zero denominators are handled, and the "Type M > 20" claim does not survive
any of them. The claim that corrected Type S (9.85%) exceeds ecology and evolution
(5–8%) also reverses under the sign-preserving variant.

**What I would like your view on:** was the magnitude-only gate intended to admit
sign reversals, or was the implicit assumption that the corrected estimate stays on
the same side of zero? And would you be comfortable with the minimal fix — not
reporting Type M and Type S for meta-analyses whose corrected mean is not
distinguishable from zero, and reporting the floored variant as the headline with
the unfloored numbers in the supplement? That keeps Yang et al.'s estimator intact
and only changes which values we are willing to interpret. My worry about doing
nothing is that the single most quotable number in the Discussion is the one least
supported by the data.

## 2. What are we summarising — a typical effect size, or a typical meta-analysis?

This is the point 3 from my original post, and it is now clearly an estimand
question rather than a dependence-modelling detail. For primary-study power:

| model | estimate | 95% CI |
|---|---|---|
| `(1 \| study_ID)` — current, as in Yang et al. | 17.15% | 16.43–17.91% |
| meta-analyses weighted equally, treated as fixed | 22.37% | 21.33–23.46% |
| `(1 \| study_ID) + (1 \| MA_model)` — 48 as a sample | 22.23% | 17.28–28.60% |

The point estimate moves because equal weighting of meta-analyses replaces
effective weighting by size; the interval widens ~6× only if we also treat the 48
as a sample from a broader population. Meta-analysis identity takes most of the
between-cluster variance (0.77 vs 0.16), and that survives every subset I tried, so
it is real structure rather than an artifact of the messy `study_ID` labels. No fit
is singular, and nested and crossed specifications agree.

The conclusion is unchanged either way — power is far below 80% under all three —
but the number and its precision are not. My inclination is still to keep the Yang
model as primary for comparability and report the others as sensitivity analyses,
**but** I would rather hear which quantity you think we are actually claiming,
since the paper's framing ("28 papers, 48 models") reads like a finite population.

## 3. Meta-analysis-level results — keep in the main text?

Following our discussion I had planned to retain these with a much more explicit
interpretation. Two things from the audit make me want to check that again:

* the aggregation weights by number of effect sizes, and `k` is strongly correlated
  with precision (r = −0.77 on the log scale), so the five largest meta-analyses
  carry 62% of the weight and four of them sit at power ≈ 1.00. The reported 82.2%
  spans **57%–89%** depending only on weighting and transformation;
* meta-analysis-level power is a monotone transformation of the pooled Wald
  statistic (Spearman = 1.000), as we discussed.

I still think there is a case for keeping it for cross-field comparability, but I
am no longer sure it belongs in the abstract. Would you prefer main text with the
caveat, or supplement with the comparison to Yang et al. retained in the
Discussion?

## 4. Type S at the meta-analysis level

The log-scale model needs an offset because Type S can be zero, and at the
meta-analysis level the offset (0.025) is larger than the estimate itself
(`exp(intercept)` = 0.0256, giving 0.00058 after subtraction). Across defensible
offsets the point estimate moves by a factor of ~500, and the lower CI limit is
negative for any offset above 0.001 — which is why Table S1 shows `0*`. The same
constraint is doing real work on the power CIs too: the upper limits for lnRR and
SMD are genuinely 1.64 and 1.04 before being constrained to 1.

Options are a logit-scale model (no offset needed), reporting the raw median and
IQR, or reporting threshold counts (28% of effect sizes exceed 5% Type S). Do you
have a preference? Meta-analysis-level Type S is effectively zero on any of them,
so this may be an argument for reporting it only as "negligible" and not as a
number.

## What I am not asking about

For completeness, these came out of the audit but I do not think they need your
input: the leave-one-out analysis is done and self-inclusion turns out to be
negligible for the pooled estimate (17.15% → 17.22%), though individual clusters in
small meta-analyses move by up to 30 percentage points, so I will report it with
that caveat rather than as a flat "negligible". On Reviewer 1's
correction-dependence question, the honest answer is that the five most influential
meta-analyses also hold 51% of the effect sizes, so the correction is broad in
coverage and only mildly concentrated in magnitude; I will answer him with the
counterfactual (dropping those five moves corrected meta-analysis-level power from
0.45 to 0.53) rather than with a share statistic.

If it is easier, I am happy to go through 1 and 2 in a call rather than here — they
are the only two that change the results rather than the wording.
