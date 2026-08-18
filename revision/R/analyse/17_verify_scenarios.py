#!/usr/bin/env python3
"""Step 16b: independent re-derivation of revision/results/assumed_effect_scenarios.csv.

Run `Rscript revision/R/analyse/16_export_scenario_inputs.R` first.

Nothing here is imported from the R side except the fitted meta-analytic inputs
(assumed effects, standard errors, cluster labels, effect-size metrics). The
scenario definitions, the three design-analysis metrics, the REML random-intercept
summary and the weighted least-squares summary are all re-implemented from the
mathematics, in another language, with no shared library.

The three metrics (Gelman & Carlin 2014), for t = |mu| / se and critical value c:

    power  = 2 - Phi(c - t) - Phi(c + t)
    Type S = Phi(-c - t) / [ (1 - Phi(c - t)) + Phi(-c - t) ]
    Type M = E[|X| | |X| > c*se] / |mu|,  X ~ N(mu, se)

The primary-study-level summary is exp(intercept) - offset from a REML fit of
log(value + offset) ~ 1 + (1 | cluster). For the intercept-only random-intercept
model the profile likelihood reduces to one dimension in lambda = sigma_b^2/sigma_e^2:

    mu_hat(lambda) = sum_i w_i ybar_i / sum_i w_i,      w_i = n_i / (1 + n_i lambda)
    -2 logL_REML   = (N-1) log(sigma_e^2) + sum_i log(1 + n_i lambda)
                     + log(sum_i w_i) + (N-1) + (N-1) log(2 pi)

with sigma_e^2 profiled out as Q/(N-1), Q the GLS residual sum of squares. That is
solved here by golden-section search on log(lambda), which shares no code with
lme4's penalised-least-squares route.

The meta-analysis-level summary is exp(intercept) - offset from a weighted least
squares fit of log(value + offset) ~ 1 with weights k, i.e. a k-weighted geometric
mean, with the interval from the usual weighted-residual variance.
"""
import csv
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))          # <repo>/revision/R/analyse
REVISION = os.path.dirname(os.path.dirname(HERE))          # <repo>/revision
TMP = os.path.join(REVISION, "results", "intermediate")
OUT = os.path.join(REVISION, "results", "assumed_effect_scenarios.csv")
# Fail loudly rather than silently reading nothing. These paths are derived from this
# file's own location, so moving the script without moving the results - or the reverse -
# would otherwise produce an empty comparison that reports success.
for _p in (TMP,):
    if not os.path.isdir(_p):
        raise SystemExit(
            "expected directory not found: " + _p +
            "\nThis script must live in <repo>/revision/R/analyse and the results in "
            "<repo>/revision/results.")

CRIT = 1.959963984540054          # qnorm(0.975)
TYPE_S_OFFSET = 0.025
# Student t 97.5% quantiles are needed for the optimistic scenario only as a check;
# the confidence limits themselves are read from the fitted models via the export.


# --- normal distribution, series/continued-fraction erf, no scipy -------------
def _erf(x):
    # Abramowitz & Stegun 7.1.26 is only ~1e-7; use the incomplete gamma series
    # instead so that Phi is accurate to full double precision in the tails.
    return math.erf(x)


def pnorm(x):
    return 0.5 * (1.0 + _erf(x / math.sqrt(2.0)))


def dnorm(x):
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


# --- the three design-analysis metrics ----------------------------------------
def power_cf(mu, se, crit=CRIT):
    t = abs(mu) / se
    return 2.0 - pnorm(crit - t) - pnorm(crit + t)


def type_s_cf(mu, se, crit=CRIT):
    t = abs(mu) / se
    pu = 1.0 - pnorm(crit - t)
    pl = pnorm(-crit - t)
    return pl / (pu + pl)


def type_m_cf(mu, se, crit=CRIT):
    t = mu / se
    num = (t * pnorm(t - crit) + dnorm(crit - t)
           - t * pnorm(-crit - t) + dnorm(crit + t))
    den = abs(t) * (pnorm(t - crit) + pnorm(-crit - t))
    return num / den


METRICS = {"power": power_cf, "type_M": type_m_cf, "type_S": type_s_cf}


# --- REML intercept-only random-intercept model -------------------------------
def _reml_objective(lam, groups, total_n):
    """-2 * profiled REML log-likelihood, and the GLS intercept, at this lambda."""
    sw = 0.0
    swy = 0.0
    logdet = 0.0
    for n_i, sum_i, _ in groups:
        w_i = n_i / (1.0 + n_i * lam)
        sw += w_i
        swy += w_i * (sum_i / n_i)
        logdet += math.log1p(n_i * lam)
    mu = swy / sw
    q = 0.0
    for n_i, sum_i, ss_i in groups:
        # sum_j (y_ij - mu)^2  -  lam * n_i^2 * (ybar_i - mu)^2 / (1 + n_i lam)
        dev = ss_i - 2.0 * mu * sum_i + n_i * mu * mu
        gbar = sum_i / n_i - mu
        q += dev - lam * n_i * n_i * gbar * gbar / (1.0 + n_i * lam)
    sigma2 = q / (total_n - 1)
    obj = (total_n - 1) * math.log(sigma2) + logdet + math.log(sw)
    return obj, mu, sigma2, sw


def fit_random_intercept(y, cluster):
    """REML fit of y ~ 1 + (1|cluster). Returns (intercept, se, sigma2_e, lambda)."""
    acc = defaultdict(lambda: [0, 0.0, 0.0])
    for value, g in zip(y, cluster):
        a = acc[g]
        a[0] += 1
        a[1] += value
        a[2] += value * value
    groups = [(a[0], a[1], a[2]) for a in acc.values()]
    total_n = len(y)

    # golden-section search on u = log(lambda) over a wide bracket, then check the
    # lambda -> 0 boundary explicitly (lme4 reports a singular fit there).
    lo, hi = -30.0, 12.0
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c = b - inv_phi * (b - a)
    d = a + inv_phi * (b - a)
    fc = _reml_objective(math.exp(c), groups, total_n)[0]
    fd = _reml_objective(math.exp(d), groups, total_n)[0]
    for _ in range(300):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = _reml_objective(math.exp(c), groups, total_n)[0]
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = _reml_objective(math.exp(d), groups, total_n)[0]
        if abs(b - a) < 1e-12:
            break
    best_lam = math.exp((a + b) / 2.0)
    best_obj = _reml_objective(best_lam, groups, total_n)[0]
    obj0 = _reml_objective(0.0, groups, total_n)[0]
    if obj0 <= best_obj:
        best_lam = 0.0
    obj, mu, sigma2, sw = _reml_objective(best_lam, groups, total_n)
    se = math.sqrt(sigma2 / sw)
    return mu, se, sigma2, best_lam


def aggregate_primary(values, cluster, metric):
    off = TYPE_S_OFFSET if metric == "type_S" else 0.0
    y = [math.log(v + off) for v in values]
    mu, se, _, _ = fit_random_intercept(y, cluster)
    return {
        "geometric_mean": math.exp(mu) - off,
        "ci_lower": math.exp(mu - CRIT * se) - off,
        "ci_upper": math.exp(mu + CRIT * se) - off,
    }


# --- weighted least squares, intercept only -----------------------------------
def aggregate_ma(values, weights, metric):
    off = TYPE_S_OFFSET if metric == "type_S" else 0.0
    y = [math.log(v + off) for v in values]
    sw = sum(weights)
    mu = sum(w * yi for w, yi in zip(weights, y)) / sw
    n = len(y)
    # weighted residual variance on n-1 degrees of freedom, as lm() reports it
    rss = sum(w * (yi - mu) ** 2 for w, yi in zip(weights, y))
    s2 = rss / (n - 1)
    se = math.sqrt(s2 / sw)
    tq = student_t_quantile(0.975, n - 1)
    return {
        "geometric_mean": math.exp(mu) - off,
        "ci_lower": math.exp(mu - tq * se) - off,
        "ci_upper": math.exp(mu + tq * se) - off,
    }


def student_t_quantile(p, df):
    """Two-sided t quantile by bisection on the regularised incomplete beta CDF."""
    def cdf(t):
        x = df / (df + t * t)
        ib = betainc_reg(df / 2.0, 0.5, x)
        return 1.0 - 0.5 * ib if t > 0 else 0.5 * ib
    lo, hi = 0.0, 1000.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def betainc_reg(a, b, x):
    """Regularised incomplete beta I_x(a,b) via the Lentz continued fraction."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1.0 - x) - lbeta)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _betacf(a, b, x):
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        if abs(d) < tiny:
            d = tiny
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        if abs(d) < tiny:
            d = tiny
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return h


# --- scenario definitions, re-derived -----------------------------------------
EXTERNAL = [
    ("SMD", "small (d = 0.2)", 0.2),
    ("SMD", "medium (d = 0.5)", 0.5),
    ("SMD", "large (d = 0.8)", 0.8),
    ("Zr", "small (r = 0.1)", math.atanh(0.1)),
    ("Zr", "medium (r = 0.3)", math.atanh(0.3)),
    ("Zr", "large (r = 0.5)", math.atanh(0.5)),
    ("lnRR", "10% change", math.log(1.10)),
    ("lnRR", "25% change", math.log(1.25)),
    ("lnRR", "50% change", math.log(1.50)),
]
OPTIMISTIC = "optimistic (95% CI limit farther from zero)"


def main():
    # A few source datasets carry non-UTF-8 bytes in their study labels. Cluster
    # labels are only ever compared for equality here, so decoding them
    # byte-preservingly with latin-1 keeps distinct labels distinct.
    with open(os.path.join(TMP, "verify_inputs_ma.csv"), encoding="latin-1") as fh:
        ma_rows = list(csv.DictReader(fh))
    with open(os.path.join(TMP, "verify_inputs_primary.csv"), encoding="latin-1") as fh:
        pr_rows = list(csv.DictReader(fh))
    assert len(ma_rows) == 48 and len(pr_rows) == 5740, "unexpected input sizes"

    for r in ma_rows:
        for c in ("k", "beta0", "se_beta0", "ci_lb_beta0", "ci_ub_beta0"):
            r[c] = float(r[c])
        # the optimistic assumed effect: the interval limit farther from zero
        r["mu_opt"] = r["ci_ub_beta0"] if r["beta0"] >= 0 else r["ci_lb_beta0"]
        assert (r["mu_opt"] >= 0) == (r["beta0"] >= 0), "optimistic scenario flipped sign"
    by_model = {r["MA_model"]: r for r in ma_rows}
    for r in pr_rows:
        r["sei"] = float(r["sei"])
        m = by_model[r["MA_model"]]
        r["mu_unc"] = m["beta0"]
        r["mu_opt"] = m["mu_opt"]

    computed = {}

    def add(agg, weighting, grouping, scenario, metric, res):
        computed[(agg, weighting, grouping, scenario, metric)] = res

    # primary-study level, all metrics pooled
    for label, key in (("uncorrected pooled mean", "mu_unc"), (OPTIMISTIC, "mu_opt")):
        for metric, fn in METRICS.items():
            vals = [fn(r[key], r["sei"]) for r in pr_rows]
            cl = [r["cluster"] for r in pr_rows]
            add("primary_study_level", "unweighted", "all metrics", label, metric,
                aggregate_primary(vals, cl, metric))

    # primary-study level, per effect-size metric: baseline, optimistic, external
    for est in ("lnRR", "SMD", "Zr"):
        sub = [r for r in pr_rows if r["effect_size_type"] == est]
        cl = [r["cluster"] for r in sub]
        for label, key in (("uncorrected pooled mean", "mu_unc"), (OPTIMISTIC, "mu_opt")):
            for metric, fn in METRICS.items():
                vals = [fn(r[key], r["sei"]) for r in sub]
                add("primary_study_level", "unweighted", est, label, metric,
                    aggregate_primary(vals, cl, metric))
        for e_type, label, mu in EXTERNAL:
            if e_type != est:
                continue
            for metric, fn in METRICS.items():
                vals = [fn(mu, r["sei"]) for r in sub]
                add("primary_study_level", "unweighted", est, label, metric,
                    aggregate_primary(vals, cl, metric))

    # meta-analysis level, per effect-size metric, both weightings
    for est in ("lnRR", "SMD", "Zr"):
        sub = [r for r in ma_rows if r["effect_size_type"] == est]
        ks = [r["k"] for r in sub]
        ones = [1.0] * len(sub)
        cases = [("uncorrected pooled mean", None)]
        cases += [(label, mu) for e_type, label, mu in EXTERNAL if e_type == est]
        for label, mu in cases:
            for metric, fn in METRICS.items():
                vals = [fn(r["beta0"] if mu is None else mu, r["se_beta0"]) for r in sub]
                add("meta_analysis_level", "k_effect_sizes", est, label, metric,
                    aggregate_ma(vals, ks, metric))
                add("meta_analysis_level", "equal_per_meta_analysis", est, label, metric,
                    aggregate_ma(vals, ones, metric))

    # --- compare against the R output -----------------------------------------
    with open(OUT) as fh:
        r_rows = list(csv.DictReader(fh))
    worst = {"rel": 0.0, "where": None}
    matched = 0
    unmatched_r = []
    for row in r_rows:
        key = (row["aggregation"], row["weighting"], row["grouping"],
               row["scenario"], row["metric"])
        if key not in computed:
            unmatched_r.append(key)
            continue
        matched += 1
        mine = computed[key]
        for col in ("geometric_mean", "ci_lower", "ci_upper"):
            a, b = float(row[col]), mine[col]
            denom = max(abs(a), abs(b), 1e-12)
            rel = abs(a - b) / denom
            if rel > worst["rel"]:
                worst = {"rel": rel, "where": (key, col, a, b)}

    print(f"rows in assumed_effect_scenarios.csv : {len(r_rows)}")
    print(f"rows independently re-derived here   : {len(computed)}")
    print(f"rows compared                        : {matched}")
    if unmatched_r:
        print(f"rows in the R output with no Python counterpart: {len(unmatched_r)}")
        for k in unmatched_r[:5]:
            print("   ", k)
    key, col, a, b = worst["where"]
    print(f"worst relative disagreement          : {worst['rel']:.3e}")
    print(f"   at {key} [{col}]  R={a:.10g}  Python={b:.10g}")

    if matched != len(r_rows):
        print("\nFAIL: not every row in the R output was checked")
        return 1
    # 1e-6 is generous for the point estimates and accommodates the different
    # optimiser paths to the REML variance components, which move the interval
    # limits in the last few digits.
    if worst["rel"] > 1e-6:
        print("\nFAIL: disagreement exceeds 1e-6")
        return 1
    print("\nPASS: every scenario row reproduces to better than 1e-6 relative.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
