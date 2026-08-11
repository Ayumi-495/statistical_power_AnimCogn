# Verification panel — falsification briefs (ready to fire)

Ayumi's instruction (2026-08-10): rerun the verification panel **from scratch**
once the agent spend limit is lifted. Reviewers must inspect the code and outputs
themselves and must **not** take the analyst's conclusions as given. Their task is
to **falsify** the six claims below, not to confirm them.

Panel, per `AYUMI/06_Agents/orchestration protocol.md`: **Lovelace** (code and
reproducibility) and **Wald** (adversarial re-derivation) and **Turing**
(statistics and formal reasoning) in parallel; then **Wald** receives Turing's
report to audit it; then **Nightingale** as the readiness gate; Pancho
consolidates. High-stakes quantitative output, so Turing plus Wald plus
Nightingale are all mandatory.

## Standing constraints for every reviewer

- Read-only on `S2_v2.R`, the manuscript, `R/`, `outputs/`, `docs/`. No commits.
- Do **not** re-run `R/01`, `R/03`, `R/04`, `R/05` in place — they overwrite the
  `outputs/` artifacts under audit. Copy the worktree elsewhere to run pipeline
  code.
- Own scratch directory per reviewer under `/private/tmp/claude-501/`.
- Report proposed fixes; do not apply them.
- **D3, D4, D7, D8 must stay open.** Supply evidence, never a decision.
- The analyst found two errors in their own work (claims 2 and 3 below). Treat the
  *corrected* versions as claims to falsify too — do not assume a self-corrected
  number is right.

## Paths

- MAIN: `/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/statistical_power_AnimCogn`
  — `S2_v2.R`, `docs/Yefeng0920-EcoEvo_PB-c6b1d62/EcoEvo_PB_script.Rmd`,
  `docs/AnimalCogn_statistical power_BiologyOpen.docx`, data in `lnRR/`, `SMD/`,
  `SMD/des_stat/`, `zr/`.
- WORKTREE under audit: `<MAIN>/.claude/worktrees/ayumi-methodological-audit-3061e3`
  — `R/00_setup.R`…`R/05_loo_study.R`, `R/README.md`, `outputs/` (29 files),
  `docs/02_audit_and_revision_plan.md`, `docs/03_manuscript_corrections.md`,
  `docs/04_verification_audit_provisional.md`.
- Sibling evidence repo: `Ayumi-495/systematic_mapping_AnimCogn`, file
  `data/papers/primary_paper.xlsx`, reachable with `gh api`.

---

## Claim 1 — the `~` defect and its downstream numbers

**Assert:** `S2_v2.R:544` (`mod = sei + year_pub.l`) and `S2_v2.R:554`
(`mod = var + year_pub.l`) omit the `~`, so `metafor` fits ONE composite
moderator equal to the arithmetic sum of the two variables rather than two
separate moderators. Only the lnRR scenario-1 block is affected; its sole member
is `MA09.csv` (1,297 of 5,740 rows, 126 `study_ID`). `beta0_c2` = 0.1060 as
written vs **0.0681** intended. Downstream, corrected side only: MA power
0.4486→0.3904 (−13.0%), MA Type M 2.040→2.180, MA Type S 0.01215→0.01223,
primary power 0.09057→0.08774, primary Type M 7.791→8.124, primary Type S
0.09854→0.10431; uncorrected side unchanged. The defect is **ours**, absent from
Yang et al.

**Try to falsify by:** proving `mod` does NOT partial-match `mods` (in which case
the moderator was silently *ignored* and `beta0_c2` would equal `beta0` — check
whether 0.1060 equals the intercept-only `beta0` of 0.2194; it does not, but
verify); showing more than one dataset sits in lnRR scenario 1; showing the same
defect elsewhere in `S2_v2.R` or present in Yang's Rmd; showing that `R/`'s
`as.formula(paste(...))` construction produces an unintended formula for some
scenario × moderator combination; showing the downstream deltas are wrong or that
"every corrected value moves in the direction that strengthens the paper" is
false or is not the whole truth (e.g. does MA-level Type S really move
meaningfully at +0.6%?).

---

## Claim 2 — the CORRECTED B1 concentration result

**Assert:** `R/04_sensitivity.R:54` sorts by `desc(shrinkage_pct)`, so the
reported `top5_share_of_kweighted_shrinkage` = **16.97%** answers the wrong
question. Ranked by `k`-weighted contribution the top five share is
**58.86%**, from **MA09, MA08, MA26, MA34, MA28_02**. Therefore the correction is
**concentrated**, not pervasive, and the earlier conclusion was wrong. The
statistic is well defined because 0 of 48 meta-analyses have negative shrinkage
(the one-directional gate guarantees it).

**Try to falsify by:** recomputing 58.86% independently; questioning whether
"top five by `k`-weighted contribution" is itself the right diagnostic for
Reviewer 1's question (Reviewer 1 asked whether "a few studies with heavily
corrected estimates are driving the main results" — is a share-of-total-shrinkage
statistic the right operationalisation at all, given the aggregation is a
`k`-weighted mean of *log* metrics, not a mean of shrinkages?); constructing a
better diagnostic (e.g. leave-those-five-out and report the change in the actual
reported medians, which is the quantity reviewers care about) and reporting what
it gives; checking whether 0 negative shrinkages is genuinely guaranteed or an
empirical accident; and checking whether MA09's contribution is inflated by the
`~` defect, which would make claims 1 and 2 mutually entangled.

---

## Claim 3 — B5 implementation, and aggregate versus individual

**Assert:** the exclusion unit is the whole `(study_ID x dataset)` cluster, not a
single row. All 1,415 clusters are present in
`outputs/B5_leave_one_study_out_rows.csv` (5,740 rows), nothing was silently
dropped, and no dataset has fewer than two distinct `study_ID` so neither `next`
branch in `R/05_loo_study.R` fired. Aggregate: median primary power
17.154%→17.224% (+0.07 pp). Individual: paired per-row power change ranges
−30.06 pp to +29.50 pp, with 21.88% of rows moving >1 pp, 3.24% >5 pp, 0.73%
>10 pp, mean |change| 0.873 pp. Two sign flips: MA15_04 `Wu_etal_2019`
(−0.1213→+0.5122) and MA39_2 `Diehl_2014` (+0.0520→−0.1365). B5 is a **lower
bound** on self-inclusion because `study_ID` is finer than "primary study" in 9
of 28 papers.

**Try to falsify by:** reading `R/05_loo_study.R` and proving a cluster was
dropped, or that the `keep` logic removes something other than the intended
cluster, or that `focal <- d[!keep, ]` and the refit are misaligned; checking
whether `sum(keep) < 2 || n_distinct(...) < 2` could ever fire and be
miscounted; verifying the ±30 pp tails are real rather than an artifact of
near-zero `mu` (if the leave-one-out `mu` approaches zero, power collapses — is
that a genuine self-inclusion effect or a numerical instability?); assessing
whether the analyst's revised claim ("negligible for the pooled estimate, not for
individual studies") is now *under*-stated or *over*-stated; and deciding whether
the two sign flips invalidate those clusters' Type S and Type M entirely.

---

## Claim 4 — the extended dependence model and the large CI change

**This is the claim Ayumi most wants stress-tested. It is D4 evidence only — do
not decide D4.**

**Assert:** for primary-study power under the uncorrected assumed effect,
`lmer(log(power) ~ 1 + (1|study_ID))` gives 17.15% [16.43, 17.91] with variance
components study 0.534 and residual 0.149; adding `(1|case)` gives **22.23%
[17.28, 28.60]** with study 0.160, **case 0.766**, residual 0.123 — a 5 pp shift
in the point estimate and a **7.6-fold** widening of the CI. Interpretation
offered: meta-analysis identity absorbs most of the variance, so the current
interval is materially overconfident.

**Try to falsify by:** checking whether `(1|study_ID) + (1|case)` is the correct
crossed specification given that `study_ID` labels recur across datasets and are
non-harmonised across four schemes — would `(1|case/study_ID)` be more
appropriate, and does the answer change the numbers?; testing whether the
variance reallocation is genuine structure or an artifact of the identifier
problem (e.g. refit on the subset of datasets with clean author-year labels, or
after collapsing `study_ID` to `case:study_ID` so it is unambiguously nested);
checking whether 48 levels of `case` support the term and whether the fit is
singular or boundary; asking whether the widened CI is *correct* rather than
merely wider (is the simple model's CI the wrong one, or is the extended model
over-parameterised for this data?); checking whether Wald CIs are appropriate here
at all versus profile or bootstrap; and reporting the same comparison for Type M,
Type S, and the corrected assumed effect, since only power was examined.

Also assess: does the extended model change the *substantive* conclusion, or only
the number? And which specification would Yang et al.'s framework imply?

---

## Claim 5 — joins, grouping, optimizers, object alignment

**Assert:** the new pipeline reproduces `S2_v2.R`'s specification exactly, with
all deviations deliberate and documented, and the legacy path reproduces the
submitted manuscript on all 13 checks in `outputs/reproduction_check.csv`.

**Try to falsify by:**

- **Object alignment.** `S2_v2.R` matches models to datasets *positionally* in
  many places (`model_est_lnRR$beta0[i]` with `lnRR[[i]]`;
  `model_SMD_pb <- append(model_SMD_sei.year, model_SMD_ess.sei.year)` against
  `case = names(SMD)`). Verify the new pipeline's canonical order (lnRR, SMD
  direct, SMD computed, Zr) really matches, and that `stopifnot(identical(...))`
  guards are sufficient rather than decorative. Hunt off-by-one.
- **Optimizers.** `S2_v2.R` uses `optim` everywhere except intercept-only Zr
  (`nlminb`, `:311-319`) and the SMD detection model with the `sei` moderator
  (`nlminb`, `:350-358`), while the 5 `ess.sei` detection models use `optim`
  (`:360-368`). Verify `optimizer_for()` encodes exactly that including the
  `uses_ess` argument, by tabulating every `rma.mv` call in `S2_v2.R`. Note this
  was already wrong once during development and corrected.
- **Joins, grouping, factor levels.** Every `left_join`, `merge`, `group_by`,
  `factor()`, `pivot_*`. Especially `R/04_sensitivity.R`'s B4 block, which builds
  a per-metric `dat_m` list — verify it selects the right datasets and that
  `primary_metrics()` receives a correctly shaped object with matching row counts.
- **The tibble data-masking hazard.** A column named `mu` shadows a vector `mu`
  in later arguments of the same `tibble()` call. This bug existed and was fixed;
  verify the fix is complete in `ma_metrics`, `primary_metrics`, and every
  `tibble()`/`mutate()` in the pipeline.
- **Reproduction-gate honesty.** Extract the 13 target values from the .docx
  yourself and judge each tolerance. Is any tolerance wide enough to hide a real
  discrepancy? Specifically, "MA type M, corrected" has tolerance 0.20 on 2.03
  and "primary type M, corrected" has 0.40 on 7.79.
- **Seed effectiveness.** Each script sources `00_setup.R` afresh, so the number
  and order of RNG draws differs between scripts. Test whether reported Type M
  values are reproducible across separate invocations, and quantify the residual
  Monte Carlo error on every reported Type M.

---

## Claim 6 — additional undocumented deviations or coding errors

**Assert:** the deviations from Yang et al. are limited to those documented — the
effective-sample-size moderator applied to SMD `des_stat` only (Yang also applied
it to lnRR; our lnRR datasets lack `C_n`/`T_n`); Yang's z-scaling computed at
`S2_v2.R:236-266` but never used in any model; no parallel scaled-data pipeline;
and `MMA_beta0_all_original` (`S2_v2.R:1131`) pooling folded `beta0` across three
effect-size scales without z-scaling.

**Try to falsify by:** systematic diff of the two implementations, looking for any
material difference not on that list; verifying the lnRR `C_n`/`T_n` absence
claim against the actual CSV columns; determining whether
`MMA_beta0_all_original` output reaches any reported figure or table; checking
`S2_v2.R:2927, 3081, 3119, 3157`, which assign `MA_case` into the Type M and
Type S figure objects **from the power object** (probably benign, unverified);
checking the 105 filtered rows for a missingness mechanism related to precision
(if low-precision rows were preferentially dropped, every power estimate is
biased upward); and searching the new pipeline for any error the analyst has not
already reported.

---

## Required output

Each reviewer uses their roster template. Wald additionally audits Turing's
review. Nightingale returns `pass` / `revise before Pancho` / `blocked`. Pancho
consolidates into four sections, per Ayumi's instruction:

1. confirmed coding facts;
2. remaining discrepancies;
3. methodological choices;
4. decisions that Ayumi and Shinichi need to make.
