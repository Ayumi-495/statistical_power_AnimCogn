# Results, Discussion, Abstract and Figure text — paste-ready

**Reads against** `docs/AnimalCogn_statistical power_BiologyOpen_rev.docx` as of
2026-08-15 15:48 (248 paragraphs, no tracked changes, two live comments).
**Supersedes** the RESULTS / ABSTRACT / DISCUSSION sections of `docs/15_manuscript_edits.md`.
The METHODS section of `docs/15` is superseded too — see §0, most of it is already done.

Every number below is from `results/revision/supplementary/TableS1_reported_metrics.csv`
and the CSVs behind it. **All of them now carry two independent derivations**: the
assumed-effect scenario rows were the last single-derivation table and were re-derived
from scratch in Python today (`R/revision/17_verify_scenarios.py`, 123 of 123 rows agree
to better than 1e-6 relative).

> **`docs/16_handoff.md` §3 has a stale block.** Its "Sensitivity to the assumed effect"
> figures (optimistic 0.2734, SMD 0.0759/0.2017/0.4082, Zr 0.0860/0.3287/0.6501, lnRR
> 0.1151/0.2514/0.4437) were computed with the *raw* study identifier and predate the
> clustering decision. Under the adopted clustering they are the values used below —
> optimistic **0.2806**, and so on. Confirmed by running `09` today: it prints
> `power, optimistic  raw 0.27335 | prefixed 0.28057`. §3 is corrected in place.

---

## 0. What the current draft already fixes

Of the three items §4 of the handoff listed as blocking:

| | item | status in the current `.docx` |
|---|---|---|
| 4.2 | Methods never disclose the log transformation | **done** — ¶0058: "We analysed power and Type M error on the log scale, then back-transformed the fitted intercepts and confidence limits" |
| 4.1 | Type S contradiction | **half done** — the Methods (¶0058) now state the offset problem and that the raw median and quartiles are the main descriptive summary. The Results (¶0070) and Discussion (¶0077) still report the offset-dominated model values as though they were the headline. §2 and §4 below finish it. |
| 4.3 | summary statistic mislabelled | **open** — the Abstract and Results still say "average", and the Methods still never name the estimator. §1 and §3 below. |

The Methods also already contain, and need no further change: the α = 0.05 / z = 1.96
statement (¶0057), the closed-form Type M statement (¶0057), the deletion of the
"average sampling variance" sentence, the whole Sensitivity analyses subsection
(¶0060–0066), the clustering definition and its limitation (¶0058), and the three
estimands (¶0065). That is the great majority of `docs/15`'s Methods section.

**Three Methods items remain** — §1. Two are small; the third (M-c) corrects what the
primary-study-level summary is a summary of, found by an external review on 15 August.

---

# 1. METHODS — three remaining edits

## M-a. Name the estimator, and put the log scale before the model that uses it

The current text introduces `lm` and `lmer` first and mentions the log scale four
sentences later, so a reader implementing it in order fits the wrong model. It also
never names the resulting summary, which is what makes "average" possible downstream.

**Current** (¶0058, first two sentences and the later log sentences):

> At the meta-analysis level, we summarised estimates of statistical power, Type M error,
> and Type S error using weighted linear regression (lm function in base R), with the
> number of effect-size observations in each meta-analytic model used as the weight. At
> the primary-study level, after calculating the metrics for each effect-size observation,
> we summarised them using linear mixed-effects models (lmer function, lme4 v1.1.37 [47])
> including study-cluster identity as a random intercept. […] We analysed power and Type M
> error on the log scale, then back-transformed the fitted intercepts and confidence limits
> to report them on the original scale.

**Replace with:**

> We summarised all three metrics on the log scale and back-transformed the fitted
> intercept and its confidence limits to the original scale, following Yang et al.
> [13,14]. At the meta-analysis level we fitted an intercept-only weighted linear
> regression (lm function in base R) to the log-transformed values, weighting each
> meta-analytic model by the number of effect-size observations it contributed. The
> back-transformed intercept of this model is a weighted geometric mean of the 48 model
> values, and we refer to it as such throughout. At the primary-study level we fitted an
> intercept-only linear mixed-effects model (lmer function, lme4 v1.1.37 [47]) to the
> log-transformed values, with study-cluster identity as a random intercept. Its
> back-transformed intercept is not a simple geometric mean, because the random intercept
> gives observations in large clusters less individual weight; we therefore refer to it
> as the back-transformed model intercept. Both summaries are smaller than the
> corresponding arithmetic mean of the same values.

Then keep the existing Type S sentences unchanged ("Type S error was also modelled on the
log scale after adding the 0.025 offset…"), and delete the now-duplicated sentence "We
analysed power and Type M error on the log scale, then back-transformed…".

*Why the label has to be level-specific:* the two estimators are genuinely different, and
the skewness of `log(power)` even changes sign between the levels (−2.69 at the
meta-analysis level, +0.71 at the primary-study level). No single word is correct for both.

## M-b. Give the sign-reversal count

**Current** (¶0050): "Since this rule is based on absolute value, the corrected estimate
does not have to keep the same sign as the uncorrected pooled mean."

**Add immediately after:**

> This occurred in 20 of the 48 models, which together contribute 1,932 of the 5,740
> effect-size observations.

This is the single most useful sentence in the revision: Reviewer 1's question about
whether the correction drives the results turns on exactly this fact, and stating it in
the Methods is much stronger than conceding it in the response letter. (Source:
`results/revision/reversal_counts.csv`; under the bias-robust estimator the count falls
to 6 of 48, which §3 reports.)

## M-c. Correct what the primary-study-level summary is a summary *of* — NEW, 2026-08-15

An external adversarial review found that the reported primary-study-level estimand was
described wrongly, in the Methods and in every table. **The number does not change** —
17.4% stands, and the model fitted is the right one, Yang et al.'s — but what we say it
means was not true.

**Current** (¶0058):

> The primary-study-level summary was fitted at the effect-size-observation level, so it
> reflects the experience of a typical effect-size observation. Meta-analytic models that
> contribute more effect-size observations have a greater influence on this summary.

**Replace with:**

> Because this model includes a study-cluster random intercept, its intercept does not
> weight each effect-size observation equally: observations in large clusters receive
> less individual weight, and the summary is close to one in which every study cluster
> counts equally. In our data the implied cluster weights span a factor of 1.2 while the
> number of effect-size observations per cluster spans a factor of 115, so this summary
> describes a typical study cluster rather than a typical effect-size observation, and a
> meta-analysis has greater influence through the number of study clusters it contributes
> than through the number of observations.

*Why this matters and how it was found.* `lmer(log(x) ~ 1 + (1|cluster))` weights cluster
*i* by n/(1 + n·λ) with λ = τ²/σ². Here λ = 4.66 and cluster sizes run 1 to 115, so the
weights run only 0.177 to 0.214. Measured against the two candidate estimands for
uncorrected power: fitted intercept **0.17354**, equal per study cluster **0.17601**
(1.4% away), equal per effect-size observation **0.16153** (7.4% away). The gap reaches
28% in the worst cell. The old description named the estimand the model is *furthest*
from.

`20_verify_reported_numbers.R` now measures this rather than asserting it, and refuses to
run if the retired label reappears.

## M-c2. Housekeeping already flagged in the file

`clubSandwich [REF]` in ¶0066 still needs its citation, and the comment on the software
paragraph notes package versions will be refreshed from the final session info. Neither
is a decision.

---

# 2. RESULTS — replacement paragraphs

## R1. Meta-analysis-level paragraph (replaces ¶0070)

Five things change: the corrected numbers predate the missing-`~` correction; every
"mean:" value is deleted (see the box below); Type S is reported the way the figure and
Table S1 report it; the figure callouts move from Figure 2 to Figure 3; and "some error
values were over 20" becomes exact.

**Replace ¶0070 with:**

> At the meta-analysis level, the weighted geometric mean of statistical power to detect
> the uncorrected mean effect was 82.2% (95% CI: 71.7–94.2%), falling to 39.0% (95% CI:
> 29.6–51.4%) after bias correction (Figure 3a, Table S1). Type M error rose from 1.11
> (95% CI: 1.02–1.20) to 2.18 (95% CI: 1.48–3.21), and exceeded 20 in three of the 48
> models after correction, where the corrected mean lay very close to zero (Figure 3b,
> Table S1). Type S error was very low at this level throughout: the median across models
> was 7 × 10⁻⁸ before correction and 0.22% (interquartile range 0.001–5.5%) after it
> (Figure 3c, Table S1). The corresponding model-based summaries, 0.06% (95% CI: 0–0.19%)
> and 1.22% (95% CI: 0.35–2.36%), are substantially larger because the log-scale model
> requires an offset of 0.025 that exceeds almost every observed value at this level; the
> two summaries differ for that reason alone. Under either, the median probability of a
> sign error in a pooled estimate was well below 1%, although the upper quartile after
> correction reached 5.5%.

## R2. Primary-study-level paragraph (replaces ¶0071)

**Replace ¶0071 with:**

> At the primary-study level, power was low before bias correction, at 17.4% (95% CI:
> 16.6–18.1%), and fell to 9.0% (95% CI: 8.7–9.3%) when the bias-corrected mean was used
> as the assumed effect (Figure 3a, Table S1). Type M error rose from 2.89 (95% CI:
> 2.79–2.99) to 7.88 (95% CI: 7.25–8.56), exceeding 20 for 70 of the 5,740 effect-size
> observations (1.2%) before correction and for 796 (13.9%) after it (Figure 3b,
> Table S1). Type S error rose from 2.76% (95% CI: 2.57–2.97%) to 10.2% (95% CI:
> 9.6–10.8%), with medians of 1.8% (interquartile range 0.2–5.7%) and 13.6% (4.0–25.8%)
> respectively (Figure 3c, Table S1). These summaries fit study cluster as a random
> effect and therefore describe a typical study cluster. Weighting each meta-analysis
> equally instead gave 22.4% (95% CI: 21.5–23.3%) before correction and 10.5% (95% CI:
> 10.3–10.8%) after it, and treating the 48 models as a sample from a broader population
> gave 22.3% (95% CI: 17.3–28.6%) and 10.5% (95% CI: 8.4–13.2%). Power was low under all
> three summaries; the two model-level summaries are higher because they give equal weight
> to meta-analyses contributing few observations.

> **Note on the Type S clause.** Both summaries are given at this level, model-based first,
> because unlike the meta-analysis level the two are of the same order (2.76% against a
> median of 1.8%) and because the cross-field comparison in the Discussion has to be made
> on the model-based value — that is what the comparison papers report. Figure 3c draws
> the model-based value at the primary-study level and the raw median at the
> meta-analysis level, exactly as you specified, and the caption in §5 says so explicitly.

## R3. New paragraph — sensitivity analyses

Place immediately after R2. This is currently missing entirely and is the substance of
the replies to R1C12, R2C8, R2C9 and R2C13.

> The low-power result did not depend on how the assumed underlying effect was obtained.
> Under a deliberately optimistic assumed effect, set to the confidence limit farther from
> zero, primary-study power rose only to 28.1% (95% CI: 26.9–29.3%) and Type M error
> remained at 2.02. Under externally specified effects, primary-study power was 7.7%,
> 20.5% and 41.4% for standardised mean differences of 0.2, 0.5 and 0.8; 8.7%, 33.1% and
> 65.3% for correlations of 0.1, 0.3 and 0.5; and 12.5%, 27.3% and 46.8% for log response
> ratios corresponding to 10%, 25% and 50% change (Table S1). Replacing the
> regression-based correction with the bias-robust estimator of Yang et al. [61] gave a
> primary-study power of 13.4% (95% CI: 12.8–14.0%) and a Type M error of 3.86, and
> reduced the number of models whose corrected mean opposed the uncorrected mean from 20
> to 6 of 48. Primary-study power therefore remained below 30% under every assumed effect
> except the largest conventional values, whereas Type M error ranged from 2.0 to 7.9
> depending on that choice.
>
> Excluding each study cluster from the pooled mean used to evaluate it, so that the
> assumed effect no longer contains the observations being assessed, changed primary-study
> power only from 17.35% (95% CI: 16.6–18.1%) to 17.45% (95% CI: 16.7–18.2%) and Type M
> error from 2.89 to 2.88. Dependence between an observation and the mean used to evaluate
> it therefore does not account for the low power we report.
>
> *(Give this pair to two decimal places — 17.35% and 17.45% — rather than 17.4% and
> 17.5%. The change is 0.096 percentage points and straddles a rounding boundary, so the
> one-decimal version invites the reader to compute a difference of 0.1 in a sentence
> whose point is that nothing moved. "by less than 0.1 percentage points" also works.)*
>
> The meta-analysis-level summary was less informative than its value alone suggests.
> Because the assumed effect and its standard error come from the same fitted model, this
> quantity is a monotone function of the ratio between them. When the assumed effect was
> instead specified from outside the corpus, the same 48 models had power of 89.2%, 97.4%
> and 98.6% against conventionally medium effects for standardised mean differences,
> correlations and log response ratios respectively (Table S1). The low corrected value
> therefore reflects corrected means lying close to zero rather than imprecise
> meta-analyses. This summary was also sensitive to weighting and to individual
> contributions: giving each meta-analysis equal weight changed corrected power from 39.0%
> to 25.0%, and omitting the single meta-analysis contributing the most observations
> (1,297 of 5,740) changed the bias-robust summary from 47.9% to 69.3%. That leverage is a
> property of weighting by observation count rather than of that paper being unusual —
> with each meta-analysis weighted equally, the largest single-paper influence was a
> different paper at 8.4%, and the median absolute change across all 28 papers was 1.4%.
> We therefore report the meta-analysis-level results as a descriptive summary rather than
> as a principal result.

## R4. One-word change in the composition paragraph (¶0069)

The sentence "The 28 included meta-analytical papers were concentrated in particular
parts of animal-cognition research (Table 1)" points at the recommendations table. See
§6 — this is a numbering decision, not a wording one.

---

> ### Why every "mean:" value is deleted
>
> The manuscript reports, for example, "82.2% (95% CI: 71.7%–94.2%; **mean: 114%**)". A
> statistical power of 114% is impossible, and a reviewer will see it immediately.
>
> That figure is `exp(μ̂ + σ̂²/2)` — the arithmetic mean of a lognormal variable with the
> fitted intercept and residual variance (`S2_v2.R:1391`, and the same construction at
> :1582, :2003, :2420, :2700). It is a legitimate estimator of the arithmetic mean of a
> lognormal quantity, but power, Type S error and any probability are bounded above,
> lognormal variables are not, and the estimator is free to exceed the bound. It does.
>
> If an arithmetic mean is wanted, the plain arithmetic mean is in
> `results/revision/assumed_effect_scenarios.csv` (`arithmetic_mean_unweighted`): 22.7%,
> not 23.1%, for uncorrected primary-study power. The replies to R1C13 and R2C12 already
> undertake to report the model-based summary, so the simplest course is to drop the
> parenthetical entirely, which is what R1 and R2 above do.

---

# 3. ABSTRACT

**Current:** "Average primary-study power declined from 17.2% to 9.1%, while average
meta-analysis-level power declined from 82.2% to 44.9%."

**Replace with:** "Primary-study power declined from 17.4% to 9.0%, and
meta-analysis-level power from 82.2% to 39.0%."

The word "average" goes because the quantities are a back-transformed model intercept and
a weighted geometric mean, and because the arithmetic means are different numbers. Check
the Abstract for any other occurrence.

**Current:** "Type S errors were low at the meta-analysis level but increased among
primary studies after correction." — accurate under the new numbers; keep.

**Also current:** "when evidence of publication bias was detected, applied bias correction"

This does not describe the procedure: a correction model is fitted for every
meta-analysis, which moderators enter depends on the direction of their slopes rather
than on any significance test, and the corrected estimate is then used only when it is
smaller in absolute magnitude.

**Replace with:** "we estimated a bias-corrected mean for each meta-analysis and used it
in place of the uncorrected mean when it was smaller in magnitude".

---

# 4. DISCUSSION

Four paragraphs carry Results numbers. Each replacement below is minimal — only the
figures and the sentences that depend on them change.

## D1. ¶0075, comparison with other fields

> …estimated power among the models included here was 82.2% using uncorrected mean
> effects and declined to **39.0%** using bias-corrected means (Table S1).

and later in the same paragraph:

> Type M error exceeded 20 in three of the included meta-analytic models after bias
> correction (Figure **3**b).

(The "three models" claim is correct as written — verified against
`results/revision/model_level_metrics.csv`: 3 of 48 under the Yang-2023 correction, none
before correction.)

## D2. ¶0076, primary-study level

> …estimated power among the studies represented in the included meta-analyses was
> **17.4%** using uncorrected mean effects and only **9.0%** using bias-corrected means
> (Table S1). […] Type M error also increased substantially, from **2.89** using
> uncorrected means to **7.88** using bias-corrected means, exceeding 20 for **13.9% of
> the 5,740 effect-size observations** (Figure **3**b, Table S1).

Note "averaged" should also go from this sentence, for the reason given in §3.

## D3. ¶0077, Type S — needs one added clause, not just new numbers

The paragraph compares our Type S values against Yang et al. Those published values are
**model-based** summaries. If our raw medians are substituted, the comparison stops being
like-for-like. Carry the comparison on the model-based values and say so.

> In contrast to our findings for power and Type M error, Type S error among the
> meta-analytic models included here was very low, with a median of 7 × 10⁻⁸ using
> uncorrected means and 0.22% using bias-corrected means (Table S1). Among the primary
> studies represented in these meta-analyses, the median rose from 1.8% to 13.6% after
> bias correction (Table S1). For comparability with published estimates, which are
> reported as model-based summaries, the corresponding model-based values here are 0.06%
> and 1.22% at the meta-analysis level and **2.76%** and **10.2%** at the primary-study
> level. On that basis, meta-analysis-level Type S error was below values reported for
> ecology and evolution [14], whereas the primary-study-level value after correction
> exceeded typical estimates in ecology and evolution (5–8% [14]) and global change
> biology (<5% [13]).

The rest of the paragraph — the caveats about what a sign reversal means substantively —
is unaffected and should be kept as written. It is one of the better passages in the
manuscript.

## D4. A framing sentence to add

The framing decision has changed: the meta-analysis-level results are now a descriptive
summary rather than a principal result, and the Discussion currently gives them equal
standing. Add at the end of ¶0075 or the start of ¶0076:

> We give more weight to the primary-study-level results throughout. At the
> meta-analysis level the assumed effect and the standard error it is judged against are
> outputs of the same fitted model, so the resulting metrics restate the model's own
> effect-to-uncertainty ratio; the summary is also sensitive to how meta-analyses are
> weighted and to the largest contributing meta-analysis.

## D5. ¶0081 already cites Yang et al. [61] as a recommendation

Since the bias-robust approach is now a reported sensitivity analysis rather than only a
recommendation, the sentence "Future work should consider regression-based alternatives
… such as the robust estimation approach recently proposed by Yang et al. [61]" reads
oddly beside a Results paragraph that applies it. Suggest: "…such as the robust
estimation approach of Yang et al. [61], which we applied here as a sensitivity analysis
and which has performed reliably across hundreds of meta-analyses in ecology and
evolution."

---

# 5. FIGURE 3 — the caption must be replaced entirely

The current caption (¶0149) describes a figure that no longer exists: it describes bar
rows per meta-analysis model grouped by effect-size metric, colour ramps in red, blue and
yellow-green, a "20+" category, two series, and an overall mean line. The replacement
figure has no tile rows, three series, no colour ramp, and summary bars read from the
canonical tables.

**Replace ¶0149 with:**

> **Figure 3. Statistical power, Type M error and Type S error under uncorrected and
> bias-corrected assumed effects, at the primary-study and meta-analysis levels.**
> (a) Statistical power, (b) Type M error and (c) Type S error. Each panel shows the
> primary-study level on the left, where every one of the 5,740 effect-size observations
> is evaluated against its own sampling standard error, and the meta-analysis level on
> the right, where each of the 48 meta-analytic models is evaluated against the standard
> error of its pooled estimate. Within each panel the three groups are the uncorrected
> pooled mean, the Yang et al. [14] bias-corrected mean used as our primary correction,
> and the Yang et al. [61] bias-robust estimate used as a sensitivity analysis. Violins
> show the distribution of the plotted values; each is scaled to the same maximum width,
> so a violin's shape shows relative density within that group and widths cannot be
> compared between groups. Points are the individual values, drawn smaller and more
> transparent at the primary-study level because there are 5,740 of them against 48.
> Horizontal bars are the summaries reported in the Results: at the primary-study level
> the back-transformed intercept of a linear mixed-effects model with a study-cluster
> random intercept, fitted over 1,415 study clusters; at the meta-analysis level the
> geometric mean across the 48 models, weighted by the number of effect-size observations
> each contributes. Meta-analysis-level Type S error is the one exception: there the bar
> is the median across models and the vertical line the interquartile range, because
> almost every value at that level falls below the offset of 0.025 that the log-scale
> model requires, so the model-based summary would reflect the offset more than the data.
> Those model-based estimates are retained in Table S1. In panel (c) the meta-analysis-level
> bars and intervals for the uncorrected and bias-robust effects are indistinguishable
> from zero at this scale (medians 7 × 10⁻⁸ and 8 × 10⁻⁷). Type M error is plotted on a
> logarithmic scale and is unbounded as the assumed effect approaches zero, because the
> assumed effect is its denominator. Type S error is bounded above at 0.5, the value it
> takes when the assumed effect is zero, and is not plotted on a logarithmic scale
> because values can reach or approach zero. Colour distinguishes the three assumed
> effects and repeats information already carried by position on the horizontal axis.

**Three things in this caption are load-bearing and should not be trimmed away.**

*"cannot be compared between groups"* — the violins use equal maximum width, not equal
area. Without this sentence a reader will read a wide violin as more data.

*"a study-cluster random intercept, fitted over 1,415 study clusters"* — this is the
estimand correction of §M-c. The bar is **not** an average over the 5,740 points drawn
beside it, and saying which model produced it is what stops that misreading.

*the panel (c) sentence* — two of the six meta-analysis-level Type S bars sit on the axis
and one has an interval too narrow to see. A reader who cannot tell a drawn zero from a
missing bar will assume the latter, which is exactly the failure mode of the submitted
figure's censored tiles.

If the caption must be shortened, cut the colour sentence first, then the sentence
beginning "Points are the individual values".

**A change already made to the code, for the record.** `R/revision/11_main_figure.R`
previously drew six lines of caption text into the image itself. Journals typeset
captions separately, so that text would have appeared twice in the article. It has been
removed; the caption above is now the only one, and the figure files have been
regenerated. Anything a reader needs in order to decode the panels must therefore be in
the caption above.

---

# 6. FOUR THINGS THAT ARE YOUR DECISION, NOT MINE

None of these can be settled from the analysis; each changes what gets written.

## 6.1 The meta-analysis-level confidence interval rests on an assumption that fails

`confint()` on the weighted fit treats the number of effect-size observations, *k*, as an
inverse-variance weight. That assumption does not hold here: the weighted residual
variance by *k* tertile is 20.2 / 15.8 / 43.4, so the largest meta-analyses are the most
variable, not the least. A nonparametric bootstrap over the 48 models (B = 20,000) gives
[65.2%, 91.4%] for uncorrected meta-analysis-level power, against the reported
[71.7%, 94.2%]. The reported interval is both narrower and shifted upward: its width on
the log scale is 0.273 against the bootstrap's 0.338, so it is about 19% narrower.

*(`docs/16` §4.4 records this as "about 28% narrower". Recomputing the two widths
directly from the intervals gives 19%, so **quote neither percentage without rerunning
the bootstrap** — the interval endpoints are what should be reported if this goes in at
all. The bootstrap itself remains single-derivation.)*

R1 above keeps the reported interval, because that is what Table S1 and the submitted
supplement contain. **Keeping it is a choice**, so: report as is; report as is with one
sentence in the Methods noting that the interval assumes *k* behaves as an
inverse-variance weight; or report the bootstrap interval instead. My recommendation is
the middle one — the point estimate is unaffected, the meta-analysis level is now a
descriptive summary anyway, and switching intervals would put Table S1 out of step with
the submitted supplement for no gain in the conclusion. This is single-derivation so far;
if you want the bootstrap in the paper it needs a second derivation first.

## 6.2 Table numbering — decide once and apply to every consumer

The evidence-base table (28 papers) is cited in the Results as "Table 1", but the only
Table 1 caption in the document is the recommendations table, which the Discussion cites
as "Table 2", and no Table 2 caption exists. Two coherent resolutions:

- **(b) Evidence base stays supplementary.** The Results citation at ¶0069 becomes
  "Table S2"; the recommendations table keeps its existing Table 1 caption; the
  Discussion's "Table 2" at ¶0082 becomes "Table 1". **Three words in the manuscript, no
  script changes, no change to your Google Sheet.**
- **(a) Evidence base becomes main-text Table 1**, recommendations become Table 2. The
  Results citation at ¶0069 is then already correct and the Discussion citation becomes
  correct — but everything downstream has to move.

**I recommend (b).** `docs/13` recommended (a) on the grounds that answering R2C11, R2C15
and R2C17 with a main-text table is stronger than answering them with a supplementary
one. That is true in itself, but it undercounted the cost: (a) is *not* the smaller edit.
It requires the same change, in the same commit, to the caption generated by
`R/revision/13_table_metadata.R`, the 14 column-dictionary rows labelled `Table S2` in
`metadata_columns.csv`, the file-index entry in `metadata_files.csv`, the filename
`TableS2_evidence_base.csv` — and to the Google Sheet tab you are part-way through
populating with 153 rows. That last one is the decisive argument: renumbering a sheet
mid-population is exactly the half-applied-decision failure we already wrote up as a
safeguard.

If you would still rather have it in the main text, say so and I will make all six
changes in one commit. But (b) costs three words.

## 6.3 Figure 2 will be cited nowhere

Once the eight `Figure 2a/b/c` callouts become `Figure 3a/b/c`, the conceptual workflow
diagram has no in-text citation at all. Journals object to that. Either add a citation —
the natural place is the first sentence of the Methods "Estimating proxies for underlying
effects" subsection, as "(Figure 2)" — or drop the figure. It is the weakest display item
in the paper and is the obvious candidate if the figure count needs reducing.

## 6.4 Two live comments are still in the file

They will be visible to the editor if the `.docx` is submitted as is:

- "Package and R version numbers will be updated from the final revision session
  information after all analyses are rerun." (14 Aug)
- A superseded Discussion paragraph held as a comment. (12 Aug)

Neither is a problem, but both should be resolved or deleted before submission.

---

# 7. Cross-references, to apply last

Per `docs/13_crossreference_fixes.md`, after the numbers are settled:

- eight `Figure 2a/b/c` → `Figure 3a/b/c` — **already handled** inside R1 and R2 above,
  plus two in the Discussion handled in D1 and D2
- `Figure S1` in the eligibility paragraph (¶0042) → `Figure 1`
- the Table 1 / Table 2 collision → §6.2
- `Table S1` is cited 11 times and correctly has no main-text caption

---

# 8. Numbers used above, with their sources

| quantity | value | source |
|---|---|---|
| MA power, uncorrected / Yang 2023 / Yang 2024 | 0.8221 / 0.3904 / 0.4789 | Table S1 part A |
| MA Type M | 1.1086 / 2.1771 / 1.5257 | Table S1 part A |
| MA Type S, raw median | 7.2e-08 / 0.0022 / 7.6e-07 | Table S1 part A |
| MA Type S, model-based | 0.00058 / 0.0122 / 0.0056 | Table S1 part A |
| MA models with Type M > 20 | 0 / 3 / 1 (of 48) | `model_level_metrics.csv` |
| MA power, equal weight, corrected | 0.2505 (25.0%) | Table S1 part A |
| primary power | 0.1735 / 0.0899 / 0.1339 | Table S1 part A |
| primary Type M | 2.8913 / 7.8787 / 3.8601 | Table S1 part A |
| primary Type S, model-based | 0.0276 / 0.1021 / 0.0484 | Table S1 part A |
| primary Type S, raw median [IQR] | 0.0183 [0.0020, 0.0573] / 0.1358 [0.0402, 0.2580] | Table S1 part A |
| primary observations with Type M > 20 | 70 (1.22%) / 796 (13.87%) of 5,740 | recomputed from `original_estimates.rds` |
| primary power, equal per meta-analysis | 0.2237 / 0.1054 | Table S1 part A |
| primary power, meta-analysis random effect | 0.2226 / 0.1052 | Table S1 part A |
| optimistic primary power / Type M | 0.2806 / 2.0166 | Table S1 part B |
| external primary, SMD 0.2/0.5/0.8 | 0.0767 / 0.2053 / 0.4144 | Table S1 part B |
| external primary, Zr 0.1/0.3/0.5 | 0.0865 / 0.3312 / 0.6531 | Table S1 part B |
| external primary, lnRR 10/25/50% | 0.1254 / 0.2726 / 0.4677 | Table S1 part B |
| external MA, medium effect, k-weighted | 0.8923 SMD / 0.9744 Zr / 0.9864 lnRR | Table S1 part C |
| leave-one-cluster-out power | 0.17354 → 0.17450; Type M 2.891 → 2.883 | `leave_one_cluster_out.csv` |
| sign reversals | 20/48 reported; 6/48 FE + VCV; 5/48 UWLS | `reversal_counts.csv` |
| effect sizes in reversing models | 1,932 of 5,740 | `reversal_counts.csv` |
| MA09 leave-one-out, FE + VCV, k-weighted | 0.4789 → 0.6931 (+44.7%), k = 1,297 | `loo_influence.csv` |
| largest equal-weighted paper influence | MA39, +8.4%; median |change| 1.4% | `leave_one_paper_out.csv` |
