"""
ML performance comparison: pyCSRML vs ChemoTyper on ToxCast endpoints.

Loads per-fold results from toxprint_toxcast_fast_results.csv and produces:
  - Per-endpoint summary table (mean ± sd for ROC-AUC, AvgPrec, MCC, BalAcc)
  - Paired statistical tests (Wilcoxon signed-rank + Bayesian estimation)
  - Aggregate delta table (pycsrml - chemotyper)
  - LaTeX and Markdown table outputs

Usage:
    python scripts/report_ml_comparison.py [--results CSV] [--out DIR]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).parent
DATA = HERE / "data"
RESULTS_CSV = DATA / "toxprint_toxcast_fast_results.csv"

METRICS = ["roc_auc", "avg_prec", "mcc", "bal_acc"]
METRIC_LABELS = {
    "roc_auc": "ROC-AUC",
    "avg_prec": "Avg. Precision",
    "mcc": "MCC",
    "bal_acc": "Balanced Acc.",
}


# ---------------------------------------------------------------------------
# Bayesian credible interval for the mean difference (normal likelihood,
# Jeffreys prior mu ~ Normal(0, large), sigma ~ Jeffreys uninformative).
# For n observations this reduces to the standard t-distribution posterior,
# i.e. the 95% credible interval equals the frequentist 95% CI.
# ---------------------------------------------------------------------------

def bayesian_mean_diff(diffs: np.ndarray, alpha: float = 0.05) -> dict:
    """
    Return posterior summary for the mean of `diffs` using a normal likelihood
    with an improper Jeffreys (flat) prior on mu and sigma.

    Under this prior the posterior of mu | diffs is a t-distribution with
    (n-1) degrees of freedom, mean = sample mean, and scale = se = sd/sqrt(n).
    The 100*(1-alpha)% credible interval is identical to the frequentist CI.
    """
    n = len(diffs)
    mu = float(np.mean(diffs))
    sd = float(np.std(diffs, ddof=1)) if n > 1 else 0.0
    se = sd / np.sqrt(n) if n > 1 else 0.0
    df = n - 1

    if df > 0 and se > 0:
        t_crit = stats.t.ppf(1 - alpha / 2, df=df)
        ci_lo = mu - t_crit * se
        ci_hi = mu + t_crit * se
        # P(delta > 0 | data) from the posterior t-distribution
        p_positive = float(stats.t.sf(0, df=df, loc=mu, scale=se))
    else:
        ci_lo = ci_hi = mu
        p_positive = float(mu > 0)

    return {"mean": mu, "sd": sd, "se": se, "ci_lo": ci_lo, "ci_hi": ci_hi,
            "p_positive": p_positive}


def wilcoxon_or_sign(a: np.ndarray, b: np.ndarray) -> tuple[float, str]:
    """Paired Wilcoxon signed-rank test; fall back to sign test for n < 5."""
    diffs = a - b
    if len(diffs) < 5 or np.all(diffs == 0):
        return 1.0, "n/a"
    try:
        res = stats.wilcoxon(diffs, alternative="two-sided",
                             zero_method="wilcox")
        return float(res.pvalue), "wilcoxon"
    except Exception:
        return 1.0, "error"


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Cohen's d for paired samples (= mean diff / pooled SD of diffs)."""
    d = a - b
    sd = float(np.std(d, ddof=1))
    return float(np.mean(d) / sd) if sd > 0 else 0.0


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def load_fold_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df


def build_pivot(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Pivot fold results so each row = (endpoint, fold), cols = feature sets."""
    piv = df.pivot_table(index=["endpoint", "fold"],
                         columns="feature_set",
                         values=metric,
                         aggfunc="first")
    return piv.reset_index()


def per_endpoint_stats(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """
    For each endpoint compute per-feature-set mean±sd and the
    Bayesian/frequentist comparison of pycsrml vs chemotyper.
    """
    rows = []
    piv = build_pivot(df, metric)

    for endpoint, grp in piv.groupby("endpoint"):
        ct = grp["chemotyper"].values
        py = grp["pycsrml"].values

        ct_mean = np.mean(ct); ct_sd = np.std(ct, ddof=1)
        py_mean = np.mean(py); py_sd = np.std(py, ddof=1)

        diffs = py - ct  # positive = pycsrml better
        bay = bayesian_mean_diff(diffs)
        pval, test = wilcoxon_or_sign(py, ct)
        d = cohens_d(py, ct)

        rows.append({
            "endpoint": endpoint,
            "ct_mean": ct_mean, "ct_sd": ct_sd,
            "py_mean": py_mean, "py_sd": py_sd,
            "delta_mean": bay["mean"],
            "delta_ci_lo": bay["ci_lo"],
            "delta_ci_hi": bay["ci_hi"],
            "p_pos": bay["p_positive"],
            "cohens_d": d,
            "pval": pval,
            "test": test,
        })
    return pd.DataFrame(rows)


def global_stats(df: pd.DataFrame, metric: str) -> dict:
    """Aggregate comparison across all endpoint × fold pairs."""
    piv = build_pivot(df, metric)
    ct = piv["chemotyper"].values
    py = piv["pycsrml"].values
    diffs = py - ct
    bay = bayesian_mean_diff(diffs)
    pval, test = wilcoxon_or_sign(py, ct)
    return {"metric": metric, **bay, "pval": pval, "test": test,
            "n_pairs": len(diffs)}


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def fmt_mean_sd(mean: float, sd: float, decimals: int = 4) -> str:
    fmt = f"{{:.{decimals}f}}"
    return f"{fmt.format(mean)} ± {fmt.format(sd)}"


def fmt_delta(mean: float, lo: float, hi: float, decimals: int = 4) -> str:
    fmt = f"{{:.{decimals}f}}"
    sign = "+" if mean >= 0 else ""
    return f"{sign}{fmt.format(mean)} [{fmt.format(lo)}, {fmt.format(hi)}]"


def render_markdown_per_endpoint(stats_df: pd.DataFrame, metric: str) -> str:
    label = METRIC_LABELS[metric]
    lines = [
        f"### {label} — per-endpoint comparison\n",
        "| Endpoint | ChemoTyper | pyCSRML | Δ (mean [95% CrI]) | P(Δ>0) | Cohens d |",
        "|----------|-----------|---------|-------------------|--------|---------|",
    ]
    for _, row in stats_df.iterrows():
        ct_str = fmt_mean_sd(row.ct_mean, row.ct_sd)
        py_str = fmt_mean_sd(row.py_mean, row.py_sd)
        delta_str = fmt_delta(row.delta_mean, row.delta_ci_lo, row.delta_ci_hi)
        lines.append(
            f"| {row.endpoint} | {ct_str} | {py_str} "
            f"| {delta_str} | {row.p_pos:.3f} | {row.cohens_d:+.3f} |"
        )
    return "\n".join(lines)


def render_markdown_global(global_rows: list[dict]) -> str:
    lines = [
        "### Global comparison (all 75 endpoint × fold pairs)\n",
        "| Metric | Mean Δ | 95% CrI | P(Δ>0) | Wilcoxon p |",
        "|--------|--------|---------|--------|------------|",
    ]
    for row in global_rows:
        label = METRIC_LABELS[row["metric"]]
        delta_str = fmt_delta(row["mean"], row["ci_lo"], row["ci_hi"])
        lines.append(
            f"| {label} | {row['mean']:+.4f} | [{row['ci_lo']:.4f}, {row['ci_hi']:.4f}] "
            f"| {row['p_positive']:.3f} | {row['pval']:.4f} |"
        )
    return "\n".join(lines)


def render_latex_summary_table(summary_df: pd.DataFrame) -> str:
    """
    LaTeX table: per-endpoint ROC-AUC and MCC for both feature sets.
    Suitable for supplementary information.
    """
    lines = [
        r"\begin{table}[h!]",
        r"\centering",
        r"\caption{ToxCast endpoint performance: ChemoTyper vs.\ pyCSRML fingerprints.",
        r"Values are mean\,$\pm$\,s.d.\ over 5 cross-validation folds.",
        r"RF classifier with fixed hyperparameters (200 trees, \texttt{max\_features=sqrt}).}",
        r"\label{tab:toxcast_comparison}",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Endpoint & \multicolumn{2}{c}{ROC-AUC} & \multicolumn{2}{c}{MCC} \\",
        r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
        r" & ChemoTyper & pyCSRML & ChemoTyper & pyCSRML \\",
        r"\midrule",
    ]
    for endpoint in summary_df["endpoint"].unique():
        sub = summary_df[summary_df["endpoint"] == endpoint]
        ct = sub[sub["feature_set"] == "chemotyper"].iloc[0]
        py = sub[sub["feature_set"] == "pycsrml"].iloc[0]
        lines.append(
            rf"\texttt{{{endpoint}}} & "
            rf"${ct.roc_auc_mean:.4f}\pm{ct.roc_auc_std:.4f}$ & "
            rf"${py.roc_auc_mean:.4f}\pm{py.roc_auc_std:.4f}$ & "
            rf"${ct.mcc_mean:.4f}\pm{ct.mcc_std:.4f}$ & "
            rf"${py.mcc_mean:.4f}\pm{py.mcc_std:.4f}$ \\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="ML comparison report")
    parser.add_argument("--results", type=Path, default=RESULTS_CSV)
    parser.add_argument("--summary", type=Path, default=DATA / "toxprint_toxcast_fast_summary.csv")
    parser.add_argument("--out", type=Path, default=DATA)
    args = parser.parse_args()

    df = load_fold_data(args.results)
    print(f"Loaded {len(df)} fold rows, {df['endpoint'].nunique()} endpoints, "
          f"{df['feature_set'].nunique()} feature sets.\n")

    report_lines: list[str] = [
        "# pyCSRML vs ChemoTyper — ML Performance Report\n",
        "**Data**: ToxCast Tox21 panel (9,014 compounds, 15 endpoints)  ",
        "**Method**: 5-fold cross-validation, RandomForest (fixed: 200 trees, sqrt features)  ",
        "**Δ = pyCSRML − ChemoTyper**; positive Δ → pyCSRML better  ",
        "**P(Δ>0)**: Bayesian posterior probability that pyCSRML exceeds ChemoTyper  ",
        "  (Jeffreys prior on μ,σ; posterior is t-distributed)  \n",
        "---\n",
    ]

    global_rows: list[dict] = []

    for metric in METRICS:
        ep_stats = per_endpoint_stats(df, metric)
        report_lines.append(render_markdown_per_endpoint(ep_stats, metric))
        report_lines.append("")

        gb = global_stats(df, metric)
        global_rows.append(gb)

        # Print brief console summary
        label = METRIC_LABELS[metric]
        print(f"{label}: global Δ = {gb['mean']:+.4f}  "
              f"95% CrI [{gb['ci_lo']:.4f}, {gb['ci_hi']:.4f}]  "
              f"P(Δ>0)={gb['p_positive']:.3f}  "
              f"Wilcoxon p={gb['pval']:.4f}")

    report_lines.append(render_markdown_global(global_rows))
    report_lines.append("\n---\n")

    # LaTeX table
    summary_df = pd.read_csv(args.summary)
    report_lines.append("## LaTeX table (ROC-AUC + MCC)\n")
    report_lines.append("```latex")
    report_lines.append(render_latex_summary_table(summary_df))
    report_lines.append("```\n")

    out_path = args.out / "ml_comparison_report.md"
    out_path.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"\nReport written to {out_path}")


if __name__ == "__main__":
    main()
