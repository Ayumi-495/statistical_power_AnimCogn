# Revision analysis scripts

Implements the plan in `docs/02_audit_and_revision_plan.md`. **`S2_v2.R` is the
immutable baseline: it is never edited, sourced, or run by these scripts.**
Everything here re-implements its documented behaviour so the revision analyses
can run with the audit corrections applied.

Run in order:

```bash
Rscript R/01_estimates.R && Rscript R/03_baseline.R && Rscript R/04_sensitivity.R && Rscript R/05_loo_study.R
```

| Script | What it does | Runtime |
|---|---|---|
| `00_setup.R` | packages, estimators, data loading, hierarchy assertions, optimizer table. Sourced by the others; not run directly. | — |
| `01_estimates.R` | 48 × (intercept-only + bias-detection + 2 correction models). Produces `beta0`, `beta0_c`, `beta0_c2`, `beta0_c3`, the scenario, the gate, and the optimistic assumed effect. Builds both the legacy and corrected specifications. | ~3 min |
| `02_metrics.R` | power / Type M / Type S at both levels for any assumed effect, and the three candidate summary statistics. Sourced; not run directly. | — |
| `03_baseline.R` | headline results, the reproduction gate, and the effect of the composite-moderator fix. | ~4 min |
| `04_sensitivity.R` | B1 correction sensitivity, B2 leave-one-paper-out, B3 optimistic assumed effect, B4 external scenarios. | ~8 min |
| `05_loo_study.R` | B5 leave-one-study-out: ~1,415 intercept-only refits. | ~30 min |

## Audit corrections applied

1. **`zr` path** — `S2_v2.R:98,103` read `path = "Zr"` while the directory is
   `zr/`. That resolves only on a case-insensitive filesystem; on Linux it
   returns `character(0)` and the pipeline silently yields 37 models instead of
   48. Fixed, and backed by `stopifnot` assertions on 48 models / 28 papers /
   5,740 effect-size rows so a silent load failure cannot pass as a valid run.
2. **Seed** — `error_M()` is Monte Carlo with `N = 10000` and neither `S2_v2.R`
   nor Yang et al. set a seed. `SEED <- 20260810` in `00_setup.R`.
3. **The lognormal mean is not reported** — `exp(intercept + 0.5·var(log x))` has
   no upper bound and returns 1.137 for meta-analysis-level power. It is computed
   in `03_baseline.R` for diagnosis and reproduction only. The principal summary
   is the model-based median with its 95% CI; original-scale arithmetic means are
   reported alongside as descriptives.
4. **Corrected means use the corrected dispersion** — a consequence of (3): the
   defect where corrected "means" were built from the *uncorrected* `var(log x)`
   disappears with the quantity itself.
5. **Composite-moderator typo** — `S2_v2.R:544,554` write `mod = sei +
   year_pub.l` and `mod = var + year_pub.l`, omitting the `~`, which fits one
   composite moderator (the arithmetic sum of a sampling variance and a centred
   publication year) instead of two. Confined to the lnRR scenario-1 block, whose
   only member is `MA09.csv` — 1,297 of the 5,740 rows. Both paths are available:
   `build_estimates(legacy = TRUE)` reproduces the submitted manuscript,
   `legacy = FALSE` (default) fits the intended model. See
   `docs/03_manuscript_corrections.md` §8.

## The reproduction gate

`03_baseline.R` checks 13 values from the submitted manuscript against the
**legacy** specification. All 13 pass (`outputs/reproduction_check.csv`).
Nothing downstream should be trusted if this fails: a mismatch means the
pipeline has drifted from the baseline, not that the baseline was wrong.

Two residual deviations are expected and documented: Type M is Monte Carlo and
was unseeded in the baseline, and primary-study power reproduces as 17.15%
against 17.20% reported (the reported mean, 23.1%, reproduces exactly, and
17.15/23.10 are mutually consistent under one intercept — most likely an
lme4/R version difference, settled by regenerating Table S1 in the declared
environment: R 4.4.2, lme4 1.1.37, metafor 4.6.0).

## Not settled — do not treat as decisions

- **Extended dependence model.** `summarise_primary(grouping = "study+model")`
  adds `(1 | case)` alongside `(1 | study_ID)`. It is computed as a comparison
  only. The reported specification remains the exact Yang et al. port,
  `(1 | study_ID)`. Pending PI discussion (plan decision D4).
- **Meta-analysis-level power.** Retained and computed, but whether it stays in
  the main text, moves to supplementary, or is renamed is pending PI discussion
  (plan §3a). The code makes no editorial choice.
- **The `beta0_c3` gate** (D3), a **corrected optimistic scenario** on
  `beta0_c2` (D7), and the **Type S CI floor** (D8) are all unresolved. The
  ingredients are computed — `outputs/B1_shrinkage_by_meta_analysis.csv` reports
  gate frequency and SE inflation — but no choice is baked in.

## Known limits

- **`study_ID` is not a harmonised primary-study identifier.** Four schemes
  coexist across the 48 datasets, and in 9 source papers `study_ID` is finer than
  "primary study". `05_loo_study.R` therefore excludes a study *as the source
  dataset defines it*, which makes B5 a **lower bound** on the self-inclusion
  effect. Resolving this needs the identifier work in plan §2.
- **Corrected leave-one-out is not implemented.** It would need ~4 fits per
  exclusion and a dropped study can flip a scenario or the gate, mixing
  self-inclusion with scenario instability. Deferred pending the uncorrected
  result (plan §4b).
- **Study-level descriptives (B6)** are a data-extraction task, not implemented
  here. `Ayumi-495/systematic_mapping_AnimCogn`,
  `data/papers/primary_paper.xlsx` supplies author/title/journal/year/DOI for all
  1,299 rows covering the 28 papers, with zero missing DOIs, as a starting point.
