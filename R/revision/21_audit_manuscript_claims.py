#!/usr/bin/env python3
"""Step 21: check every number written into the manuscript text against the result files.

Every other check in this workflow verifies the pipeline against itself. This one checks
the last step, which is the one no gate can see: a human reading a CSV, rounding a value,
converting it to a percentage and typing it into a sentence. That is where the remaining
error will be, and it has already happened once - the Results paragraph carried 44.9%
for a value the pipeline had recomputed as 39.0%.

The claims below are enumerated BY HAND from `docs/17_results_and_figure_text.md` rather
than parsed out of the prose, because a parser would silently skip anything it failed to
recognise and report success. Each entry records the sentence the number appears in, the
figure exactly as written there, and where it comes from. A claim is correct when the
stored value, converted to the units the sentence uses and rounded to the number of
digits the sentence shows, equals what is written.

WHEN THIS DISAGREES, THE DOCUMENT IS WRONG, NOT THE CSV.

Run: python3 R/revision/21_audit_manuscript_claims.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
REV = os.path.join(ROOT, "results", "revision")
SUP = os.path.join(REV, "supplementary")


def load(path):
    with open(path, encoding="latin-1") as fh:
        return list(csv.DictReader(fh))


S1 = load(os.path.join(SUP, "TableS1_reported_metrics.csv"))
MLVL = [r for r in load(os.path.join(REV, "meta_analysis_level_sensitivity.csv"))
        if not r["role"].startswith("diagnostic")]
PLVL = load(os.path.join(REV, "primary_level_sensitivity.csv"))
SCEN = load(os.path.join(REV, "assumed_effect_scenarios.csv"))
LOCO = load(os.path.join(REV, "leave_one_cluster_out.csv"))
LOO = load(os.path.join(REV, "loo_influence.csv"))
LOPO = load(os.path.join(REV, "leave_one_paper_out.csv"))
REVR = load(os.path.join(REV, "reversal_counts.csv"))
MODL = load(os.path.join(REV, "model_level_metrics.csv"))


def one(rows, field, **keys):
    """The single row matching every key, or an error naming what went wrong."""
    hit = [r for r in rows if all(r[k] == v for k, v in keys.items())]
    if len(hit) != 1:
        raise KeyError(f"{len(hit)} rows match {keys} (expected exactly 1)")
    return float(hit[0][field])


def ma(est, metric, field, weighting="k_effect_sizes"):
    return one(MLVL, field, effect_estimator=est, metric=metric, weighting=weighting)


def pr(est, metric, field, weighting="unweighted_per_effect_size"):
    return one(PLVL, field, effect_estimator=est, metric=metric, weighting=weighting)


def sc(grouping, scenario, metric, field, agg="primary_study_level", weighting="unweighted"):
    return one(SCEN, field, aggregation=agg, grouping=grouping, scenario=scenario,
               metric=metric, weighting=weighting)


UNC = "uncorrected_beta0"
Y23 = "yang2023_gated_beta0_c3"
FEV = "yang2024_FE_VCV"
OPT = "optimistic (95% CI limit farther from zero)"

# (section, sentence fragment, value as written, decimals shown, scale, stored value)
# scale 100 means the sentence writes a percentage of a stored proportion.
CLAIMS = [
    # ---- R1, meta-analysis level -------------------------------------------
    ("R1", "power to detect the uncorrected mean effect was 82.2%", 82.2, 1, 100,
     lambda: ma(UNC, "power", "geometric_mean")),
    ("R1", "uncorrected power CI lower 71.7%", 71.7, 1, 100, lambda: ma(UNC, "power", "ci_lower")),
    ("R1", "uncorrected power CI upper 94.2%", 94.2, 1, 100, lambda: ma(UNC, "power", "ci_upper")),
    ("R1", "falling to 39.0% after bias correction", 39.0, 1, 100,
     lambda: ma(Y23, "power", "geometric_mean")),
    ("R1", "corrected power CI lower 29.6%", 29.6, 1, 100, lambda: ma(Y23, "power", "ci_lower")),
    ("R1", "corrected power CI upper 51.4%", 51.4, 1, 100, lambda: ma(Y23, "power", "ci_upper")),
    ("R1", "Type M error rose from 1.11", 1.11, 2, 1, lambda: ma(UNC, "type_M", "geometric_mean")),
    ("R1", "uncorrected Type M CI lower 1.02", 1.02, 2, 1, lambda: ma(UNC, "type_M", "ci_lower")),
    ("R1", "uncorrected Type M CI upper 1.20", 1.20, 2, 1, lambda: ma(UNC, "type_M", "ci_upper")),
    ("R1", "Type M error rose ... to 2.18", 2.18, 2, 1, lambda: ma(Y23, "type_M", "geometric_mean")),
    ("R1", "corrected Type M CI lower 1.48", 1.48, 2, 1, lambda: ma(Y23, "type_M", "ci_lower")),
    ("R1", "corrected Type M CI upper 3.21", 3.21, 2, 1, lambda: ma(Y23, "type_M", "ci_upper")),
    ("R1", "median across models was 7 x 10^-8 before correction", 7.2e-8, 8, 1,
     lambda: ma(UNC, "type_S", "raw_median")),
    ("R1", "0.22% after correction", 0.22, 2, 100, lambda: ma(Y23, "type_S", "raw_median")),
    ("R1", "corrected Type S IQR lower 0.001%", 0.001, 3, 100, lambda: ma(Y23, "type_S", "raw_q1")),
    ("R1", "corrected Type S IQR upper 5.5%", 5.5, 1, 100, lambda: ma(Y23, "type_S", "raw_q3")),
    ("R1", "model-based Type S 0.06% uncorrected", 0.06, 2, 100,
     lambda: ma(UNC, "type_S", "geometric_mean")),
    ("R1", "model-based Type S 1.22% corrected", 1.22, 2, 100,
     lambda: ma(Y23, "type_S", "geometric_mean")),
    ("R1", "corrected model-based Type S CI lower 0.35%", 0.35, 2, 100,
     lambda: ma(Y23, "type_S", "ci_lower")),
    ("R1", "corrected model-based Type S CI upper 2.36%", 2.36, 2, 100,
     lambda: ma(Y23, "type_S", "ci_upper")),
    # ---- R2, primary-study level -------------------------------------------
    ("R2", "power was low before bias correction, at 17.4%", 17.4, 1, 100,
     lambda: pr(UNC, "power", "geometric_mean")),
    ("R2", "uncorrected primary power CI lower 16.6%", 16.6, 1, 100, lambda: pr(UNC, "power", "ci_lower")),
    ("R2", "uncorrected primary power CI upper 18.1%", 18.1, 1, 100, lambda: pr(UNC, "power", "ci_upper")),
    ("R2", "fell to 9.0%", 9.0, 1, 100, lambda: pr(Y23, "power", "geometric_mean")),
    ("R2", "corrected primary power CI lower 8.7%", 8.7, 1, 100, lambda: pr(Y23, "power", "ci_lower")),
    ("R2", "corrected primary power CI upper 9.3%", 9.3, 1, 100, lambda: pr(Y23, "power", "ci_upper")),
    ("R2", "Type M error rose from 2.89", 2.89, 2, 1, lambda: pr(UNC, "type_M", "geometric_mean")),
    ("R2", "uncorrected primary Type M CI lower 2.79", 2.79, 2, 1, lambda: pr(UNC, "type_M", "ci_lower")),
    ("R2", "uncorrected primary Type M CI upper 2.99", 2.99, 2, 1, lambda: pr(UNC, "type_M", "ci_upper")),
    ("R2", "Type M error rose ... to 7.88", 7.88, 2, 1, lambda: pr(Y23, "type_M", "geometric_mean")),
    ("R2", "corrected primary Type M CI lower 7.25", 7.25, 2, 1, lambda: pr(Y23, "type_M", "ci_lower")),
    ("R2", "corrected primary Type M CI upper 8.56", 8.56, 2, 1, lambda: pr(Y23, "type_M", "ci_upper")),
    ("R2", "Type S error rose from 2.76%", 2.76, 2, 100, lambda: pr(UNC, "type_S", "geometric_mean")),
    ("R2", "uncorrected primary Type S CI lower 2.57%", 2.57, 2, 100, lambda: pr(UNC, "type_S", "ci_lower")),
    ("R2", "uncorrected primary Type S CI upper 2.97%", 2.97, 2, 100, lambda: pr(UNC, "type_S", "ci_upper")),
    ("R2", "Type S error rose ... to 10.2%", 10.2, 1, 100, lambda: pr(Y23, "type_S", "geometric_mean")),
    ("R2", "corrected primary Type S CI lower 9.6%", 9.6, 1, 100, lambda: pr(Y23, "type_S", "ci_lower")),
    ("R2", "corrected primary Type S CI upper 10.8%", 10.8, 1, 100, lambda: pr(Y23, "type_S", "ci_upper")),
    ("R2", "medians of 1.8% (uncorrected)", 1.8, 1, 100, lambda: pr(UNC, "type_S", "raw_median")),
    ("R2", "uncorrected Type S IQR lower 0.2%", 0.2, 1, 100, lambda: pr(UNC, "type_S", "raw_q1")),
    ("R2", "uncorrected Type S IQR upper 5.7%", 5.7, 1, 100, lambda: pr(UNC, "type_S", "raw_q3")),
    ("R2", "medians of 13.6% (corrected)", 13.6, 1, 100, lambda: pr(Y23, "type_S", "raw_median")),
    ("R2", "corrected Type S IQR lower 4.0%", 4.0, 1, 100, lambda: pr(Y23, "type_S", "raw_q1")),
    ("R2", "corrected Type S IQR upper 25.8%", 25.8, 1, 100, lambda: pr(Y23, "type_S", "raw_q3")),
    ("R2", "equal weighting gave 22.4%", 22.4, 1, 100,
     lambda: pr(UNC, "power", "geometric_mean", "equal_per_meta_analysis")),
    ("R2", "equal weighting CI lower 21.5%", 21.5, 1, 100,
     lambda: pr(UNC, "power", "ci_lower", "equal_per_meta_analysis")),
    ("R2", "equal weighting CI upper 23.3%", 23.3, 1, 100,
     lambda: pr(UNC, "power", "ci_upper", "equal_per_meta_analysis")),
    ("R2", "equal weighting after correction 10.5%", 10.5, 1, 100,
     lambda: pr(Y23, "power", "geometric_mean", "equal_per_meta_analysis")),
    ("R2", "equal weighting corrected CI lower 10.3%", 10.3, 1, 100,
     lambda: pr(Y23, "power", "ci_lower", "equal_per_meta_analysis")),
    ("R2", "equal weighting corrected CI upper 10.8%", 10.8, 1, 100,
     lambda: pr(Y23, "power", "ci_upper", "equal_per_meta_analysis")),
    ("R2", "random effect gave 22.3%", 22.3, 1, 100,
     lambda: pr(UNC, "power", "geometric_mean", "meta_analysis_random_effect")),
    ("R2", "random effect CI lower 17.3%", 17.3, 1, 100,
     lambda: pr(UNC, "power", "ci_lower", "meta_analysis_random_effect")),
    ("R2", "random effect CI upper 28.6%", 28.6, 1, 100,
     lambda: pr(UNC, "power", "ci_upper", "meta_analysis_random_effect")),
    ("R2", "random effect corrected 10.5%", 10.5, 1, 100,
     lambda: pr(Y23, "power", "geometric_mean", "meta_analysis_random_effect")),
    ("R2", "random effect corrected CI lower 8.4%", 8.4, 1, 100,
     lambda: pr(Y23, "power", "ci_lower", "meta_analysis_random_effect")),
    ("R2", "random effect corrected CI upper 13.2%", 13.2, 1, 100,
     lambda: pr(Y23, "power", "ci_upper", "meta_analysis_random_effect")),
    # ---- R3, sensitivity analyses ------------------------------------------
    ("R3", "optimistic primary power 28.1%", 28.1, 1, 100,
     lambda: sc("all metrics", OPT, "power", "geometric_mean")),
    ("R3", "optimistic power CI lower 26.9%", 26.9, 1, 100,
     lambda: sc("all metrics", OPT, "power", "ci_lower")),
    ("R3", "optimistic power CI upper 29.3%", 29.3, 1, 100,
     lambda: sc("all metrics", OPT, "power", "ci_upper")),
    ("R3", "optimistic Type M remained at 2.02", 2.02, 2, 1,
     lambda: sc("all metrics", OPT, "type_M", "geometric_mean")),
    ("R3", "external SMD d = 0.2 gives 7.7%", 7.7, 1, 100,
     lambda: sc("SMD", "small (d = 0.2)", "power", "geometric_mean")),
    ("R3", "external SMD d = 0.5 gives 20.5%", 20.5, 1, 100,
     lambda: sc("SMD", "medium (d = 0.5)", "power", "geometric_mean")),
    ("R3", "external SMD d = 0.8 gives 41.4%", 41.4, 1, 100,
     lambda: sc("SMD", "large (d = 0.8)", "power", "geometric_mean")),
    ("R3", "external r = 0.1 gives 8.7%", 8.7, 1, 100,
     lambda: sc("Zr", "small (r = 0.1)", "power", "geometric_mean")),
    ("R3", "external r = 0.3 gives 33.1%", 33.1, 1, 100,
     lambda: sc("Zr", "medium (r = 0.3)", "power", "geometric_mean")),
    ("R3", "external r = 0.5 gives 65.3%", 65.3, 1, 100,
     lambda: sc("Zr", "large (r = 0.5)", "power", "geometric_mean")),
    ("R3", "external lnRR 10% gives 12.5%", 12.5, 1, 100,
     lambda: sc("lnRR", "10% change", "power", "geometric_mean")),
    ("R3", "external lnRR 25% gives 27.3%", 27.3, 1, 100,
     lambda: sc("lnRR", "25% change", "power", "geometric_mean")),
    ("R3", "external lnRR 50% gives 46.8%", 46.8, 1, 100,
     lambda: sc("lnRR", "50% change", "power", "geometric_mean")),
    ("R3", "bias-robust primary power 13.4%", 13.4, 1, 100,
     lambda: pr(FEV, "power", "geometric_mean")),
    ("R3", "bias-robust primary power CI lower 12.8%", 12.8, 1, 100,
     lambda: pr(FEV, "power", "ci_lower")),
    ("R3", "bias-robust primary power CI upper 14.0%", 14.0, 1, 100,
     lambda: pr(FEV, "power", "ci_upper")),
    ("R3", "bias-robust primary Type M 3.86", 3.86, 2, 1,
     lambda: pr(FEV, "type_M", "geometric_mean")),
    ("R3", "leave-one-cluster-out 17.35% self-inclusive", 17.35, 2, 100,
     lambda: one(LOCO, "geometric_mean", effect_estimator=UNC, metric="power",
                 assumed_effect="self_inclusive", weighting="unweighted_per_effect_size")),
    ("R3", "leave-one-cluster-out 17.45% held out", 17.45, 2, 100,
     lambda: one(LOCO, "geometric_mean", effect_estimator=UNC, metric="power",
                 assumed_effect="leave_one_cluster_out", weighting="unweighted_per_effect_size")),
    ("R3", "leave-one-cluster-out Type M 2.89 self-inclusive", 2.89, 2, 1,
     lambda: one(LOCO, "geometric_mean", effect_estimator=UNC, metric="type_M",
                 assumed_effect="self_inclusive", weighting="unweighted_per_effect_size")),
    ("R3", "leave-one-cluster-out Type M 2.88 held out", 2.88, 2, 1,
     lambda: one(LOCO, "geometric_mean", effect_estimator=UNC, metric="type_M",
                 assumed_effect="leave_one_cluster_out", weighting="unweighted_per_effect_size")),
    ("R3", "MA-level external medium SMD 89.2%", 89.2, 1, 100,
     lambda: sc("SMD", "medium (d = 0.5)", "power", "geometric_mean",
                agg="meta_analysis_level", weighting="k_effect_sizes")),
    ("R3", "MA-level external medium Zr 97.4%", 97.4, 1, 100,
     lambda: sc("Zr", "medium (r = 0.3)", "power", "geometric_mean",
                agg="meta_analysis_level", weighting="k_effect_sizes")),
    ("R3", "MA-level external medium lnRR 98.6%", 98.6, 1, 100,
     lambda: sc("lnRR", "25% change", "power", "geometric_mean",
                agg="meta_analysis_level", weighting="k_effect_sizes")),
    ("R3", "equal weighting changed corrected power to 25.0%", 25.0, 1, 100,
     lambda: ma(Y23, "power", "geometric_mean", "equal_per_meta_analysis")),
    ("R3", "MA09 removal changes bias-robust summary from 47.9%", 47.9, 1, 100,
     lambda: ma(FEV, "power", "geometric_mean")),
    ("R3", "MA09 removal changes bias-robust summary to 69.3%", 69.3, 1, 100,
     lambda: one(LOO, "summary_without", effect_estimator=FEV, metric="power",
                 dropped_MA_model="MA09.csv", weighting="k_effect_sizes")),
    ("R3", "largest equal-weighted paper influence 8.4%", 8.4, 1, 1,
     lambda: one(LOPO, "pct_change", effect_estimator=FEV, metric="power",
                 dropped_source_paper="MA39", weighting="equal_per_meta_analysis")),
    # ---- Methods -----------------------------------------------------------
    ("M-b", "sign reversal occurred in 20 of the 48 models", 20, 0, 1,
     lambda: one(REVR, "n_reversal", effect_estimator=Y23)),
    ("M-b", "those models contribute 1,932 effect-size observations", 1932, 0, 1,
     lambda: one(REVR, "n_effect_size_in_reversals", effect_estimator=Y23)),
    ("R3", "reversals fall from 20 to 6 of 48 under the bias-robust estimator", 6, 0, 1,
     lambda: one(REVR, "n_reversal", effect_estimator=FEV)),
]

# Counts that are computed from a file rather than read from a single cell.
COUNTED = [
    ("R1", "Type M exceeded 20 in three of the 48 models after correction", 3,
     lambda: sum(1 for r in MODL if r["effect_estimator"] == Y23 and float(r["type_M"]) > 20)),
    ("D1", "Type M exceeded 20 in no model before correction", 0,
     lambda: sum(1 for r in MODL if r["effect_estimator"] == UNC and float(r["type_M"]) > 20)),
    ("R3", "MA09 contributes 1,297 effect-size observations", 1297,
     lambda: int(one(LOO, "dropped_k", effect_estimator=FEV, metric="power",
                     dropped_MA_model="MA09.csv", weighting="k_effect_sizes"))),
]


def main():
    fails = []
    print(f"{'sec':5s} {'claim':62s} {'written':>12s} {'stored':>14s}  ok")
    print("-" * 100)
    for sec, text, written, dec, scale, getter in CLAIMS:
        try:
            stored = getter()
        except KeyError as exc:
            fails.append((sec, text, written, f"LOOKUP FAILED: {exc}"))
            print(f"{sec:5s} {text[:62]:62s} {written:>12g} {'LOOKUP FAIL':>14s}  NO")
            continue
        shown = round(stored * scale, dec)
        ok = abs(shown - written) < 10 ** (-dec) / 2
        if not ok:
            fails.append((sec, text, written, shown))
        print(f"{sec:5s} {text[:62]:62s} {written:>12g} {stored * scale:>14.6g}  "
              f"{'ok' if ok else 'NO'}")
    for sec, text, written, getter in COUNTED:
        got = getter()
        ok = got == written
        if not ok:
            fails.append((sec, text, written, got))
        print(f"{sec:5s} {text[:62]:62s} {written:>12g} {got:>14g}  {'ok' if ok else 'NO'}")

    n = len(CLAIMS) + len(COUNTED)
    print("-" * 100)
    if fails:
        print(f"\n{len(fails)} of {n} claims DISAGREE with the result files:\n")
        for sec, text, written, got in fails:
            print(f"  [{sec}] {text}\n        written {written}, stored {got}")
        print("\nThe document is wrong, not the CSV. Fix docs/17_results_and_figure_text.md.")
        return 1
    print(f"\nPASS: all {n} manuscript claims match the result files "
          f"at the precision they are written to.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
