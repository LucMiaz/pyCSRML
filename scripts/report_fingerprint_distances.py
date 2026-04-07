"""
Fingerprint fidelity analysis: pyCSRML vs ChemoTyper on ToxCast and CLinventory.

Loads the two cached fingerprint matrices (chemotyper and pycsrml) and computes:
  - Per-compound Tanimoto similarity between matched pairs
  - Per-bit accuracy (fraction of compounds where both fingerprints agree)
  - Bits with accuracy < 99% flagged for investigation

Datasets covered:
  1. ToxCast (9 014 compounds) — ToxPrint v2 (729 bits)
     Source: benchmark_toxprint_toxcast.py NPY cache
  2. CLinventory (181 745 compounds) — ToxPrint v2 (729 bits)
     Source: tests/test_data/_toxprint_V2_vs_clinventory_molecules_cache.npz
  3. CLinventory (181 745 compounds) — TxP_PFAS v1 (129 bits)
     Source: tests/test_data/_TxP_PFAS_v1_vs_clinventory_molecules_cache.npz

Usage:
    python scripts/report_fingerprint_distances.py [--out DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
DATA = HERE / "data"
PFASGROUPS_BENCH = HERE.parent.parent / "PFASGroups" / "benchmark"
TESTS_DATA = HERE.parent / "tests" / "test_data"

CLINVENTORY_TOXPRINT_NPZ = TESTS_DATA / "_toxprint_V2_vs_clinventory_molecules_cache.npz"
CLINVENTORY_TXPPFAS_NPZ  = TESTS_DATA / "_TxP_PFAS_v1_vs_clinventory_molecules_cache.npz"


def find_cache_pair(data_dir: Path) -> tuple[Path, Path]:
    """Return (chemotyper_cache, pycsrml_cache) .npy paths."""
    ct_files = sorted(data_dir.glob("chemotyper_toxprint_cache_*.npy"))
    py_files = sorted(data_dir.glob("pycsrml_toxprint_cache_*.npy"))
    if not ct_files or not py_files:
        raise FileNotFoundError(
            f"Could not find fingerprint cache files in {data_dir}. "
            "Run benchmark_toxprint_toxcast.py first."
        )
    # Match by hash suffix
    for ct in ct_files:
        suffix = ct.stem.split("_")[-1]
        for py in py_files:
            if py.stem.endswith(suffix):
                return ct, py
    return ct_files[0], py_files[0]


def tanimoto_binary(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Row-wise Tanimoto (Jaccard) similarity for binary matrices.

    T(a,b) = |a & b| / |a | b|

    Returns shape (n,) float32. Compounds where both rows are all-zero get
    similarity 1.0 (identical empty vectors).
    """
    a = a.astype(bool)
    b = b.astype(bool)
    inter = (a & b).sum(axis=1).astype(np.float32)
    union = (a | b).sum(axis=1).astype(np.float32)
    sim = np.where(union > 0, inter / union, 1.0)
    return sim.astype(np.float32)


def per_bit_accuracy(ct: np.ndarray, py: np.ndarray) -> np.ndarray:
    """Fraction of compounds where bit i agrees between ct and py."""
    agree = (ct.astype(bool) == py.astype(bool))
    return agree.mean(axis=0)


def percentile_summary(arr: np.ndarray, label: str) -> str:
    pcts = [0, 1, 5, 25, 50, 75, 95, 99, 100]
    vals = np.percentile(arr, pcts)
    parts = ", ".join(f"p{p}={v:.4f}" for p, v in zip(pcts, vals))
    return f"{label}: n={len(arr)}, mean={arr.mean():.4f}, sd={arr.std():.4f}, {parts}"


def load_bit_names(bench_dir: Path) -> list[str] | None:
    """Try to load ToxPrint bit names from the benchmark TSV header."""
    tsv = bench_dir / "test_data" / "toxprint_V2.tsv"
    if not tsv.exists():
        return None
    with tsv.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
    # First column is SMILES/InChIKey, rest are bit names
    return header[1:] if len(header) > 1 else None


def load_clinventory_npz(npz_path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Load (pred, ref, bit_names) from a concordance-test .npz cache."""
    data = np.load(npz_path, allow_pickle=True)
    pred = data["pred"].astype(np.uint8)
    ref  = data["ref"].astype(np.uint8)
    bit_names = list(data["bit_names"])
    return pred, ref, bit_names


def analyse_pair(
    ct: np.ndarray,
    py: np.ndarray,
    bit_names: list[str],
    label: str,
    bad_threshold: float = 0.99,
) -> dict:
    """Compute Tanimoto and per-bit accuracy stats for a (ChemoTyper, pyCSRML) pair."""
    n_compounds, n_bits = ct.shape
    sim = tanimoto_binary(ct, py)
    total_agree = (ct.astype(bool) == py.astype(bool)).sum()
    overall_acc = total_agree / (n_compounds * n_bits)
    bit_acc = per_bit_accuracy(ct, py)

    pcts = np.percentile(sim, [0, 1, 5, 25, 50, 75, 95, 99, 100])
    return {
        "label": label,
        "n_compounds": n_compounds,
        "n_bits": n_bits,
        "sim": sim,
        "overall_acc": overall_acc,
        "bit_acc": bit_acc,
        "bit_names": bit_names,
        "pcts": pcts,
        "sim_eq1": int((sim == 1.0).sum()),
        "sim_lt99": int((sim < 0.99).sum()),
        "sim_lt95": int((sim < 0.95).sum()),
        "bad_threshold": bad_threshold,
    }


def _tanimoto_md_block(r: dict) -> list[str]:
    """Return markdown lines for the Tanimoto stats table of one dataset result."""
    p = r["pcts"]  # [p0, p1, p5, p25, p50, p75, p95, p99, p100]
    n = r["n_compounds"]
    sim = r["sim"]
    return [
        f"| Mean Tanimoto | {sim.mean():.4f} |",
        f"| Median Tanimoto | {float(np.median(sim)):.4f} |",
        f"| Std | {sim.std():.4f} |",
        f"| p1 / p5 / p25 | {p[1]:.4f} / {p[2]:.4f} / {p[3]:.4f} |",
        f"| p75 / p95 / p99 | {p[5]:.4f} / {p[6]:.4f} / {p[7]:.4f} |",
        f"| Maximum | {p[8]:.4f} |",
        f"| Identical (T=1.0) | {r['sim_eq1']:,d} / {n:,d} "
        f"({100*r['sim_eq1']/n:.1f}%) |",
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description="Fingerprint distance analysis")
    parser.add_argument("--out", type=Path, default=DATA)
    parser.add_argument("--bad-threshold", type=float, default=0.99,
                        help="Accuracy below this is flagged (default 0.99)")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Dataset 1 — ToxCast, ToxPrint v2
    # ------------------------------------------------------------------
    ct_path, py_path = find_cache_pair(DATA)
    print(f"ChemoTyper cache : {ct_path.name}")
    print(f"pyCSRML cache    : {py_path.name}")
    ct = np.load(ct_path)
    py = np.load(py_path)
    print(f"Shapes: CT={ct.shape}, PY={py.shape}")

    bit_names_toxprint = load_bit_names(PFASGROUPS_BENCH)
    if not bit_names_toxprint or len(bit_names_toxprint) != ct.shape[1]:
        bit_names_toxprint = [f"bit_{i}" for i in range(ct.shape[1])]
    else:
        print(f"  Loaded {len(bit_names_toxprint)} bit names from TSV header")

    r_tc = analyse_pair(ct, py, bit_names_toxprint, "ToxCast — ToxPrint v2", args.bad_threshold)

    print("\n" + percentile_summary(r_tc["sim"], "Per-compound Tanimoto similarity"))
    print(f"  Identical (T=1.0) : {r_tc['sim_eq1']:6d} / {r_tc['n_compounds']}"
          f" ({100*r_tc['sim_eq1']/r_tc['n_compounds']:.2f}%)")

    total_agree_tc = (ct.astype(bool) == py.astype(bool)).sum()
    print(f"\nOverall bit accuracy: {100*r_tc['overall_acc']:.4f}%")

    # Save per-bit CSV for ToxCast
    bit_df_tc = pd.DataFrame({
        "bit_index": range(r_tc["n_bits"]),
        "bit_name": r_tc["bit_names"],
        "accuracy": r_tc["bit_acc"],
        "ct_prevalence": ct.astype(bool).mean(axis=0),
        "py_prevalence": py.astype(bool).mean(axis=0),
    }).sort_values("accuracy")
    out_bits_tc = args.out / "toxcast_per_bit_accuracy.csv"
    bit_df_tc.to_csv(out_bits_tc, index=False)
    print(f"Per-bit accuracy saved to {out_bits_tc}")

    # Print worst bits for ToxCast
    n_show = min(20, int((r_tc["bit_acc"] < args.bad_threshold).sum()) + 5)
    if n_show > 0:
        print(f"\nLowest accuracy bits (top {n_show}):")
        print(f"{'Bit name':<60} {'Accuracy':>10} {'CT prev':>8} {'PY prev':>8}")
        print("-" * 90)
        for _, row in bit_df_tc.head(n_show).iterrows():
            print(f"{str(row.bit_name):<60} {row.accuracy:>10.4f} "
                  f"{row.ct_prevalence:>8.4f} {row.py_prevalence:>8.4f}")

    # Save per-compound Tanimoto for ToxCast
    pd.DataFrame({"tanimoto_sim": r_tc["sim"]}).to_csv(
        args.out / "toxcast_per_compound_tanimoto.csv", index=False
    )
    print(f"Per-compound Tanimoto saved to {args.out / 'toxcast_per_compound_tanimoto.csv'}")

    # ------------------------------------------------------------------
    # Dataset 2 — CLinventory, ToxPrint v2
    # ------------------------------------------------------------------
    results_clinv = []
    if CLINVENTORY_TOXPRINT_NPZ.exists():
        print(f"\nLoading CLinventory ToxPrint cache: {CLINVENTORY_TOXPRINT_NPZ.name}")
        pred_clinv_tp, ref_clinv_tp, bn_clinv_tp = load_clinventory_npz(CLINVENTORY_TOXPRINT_NPZ)
        r_clinv_tp = analyse_pair(
            ref_clinv_tp, pred_clinv_tp,
            bn_clinv_tp, "CLinventory — ToxPrint v2", args.bad_threshold,
        )
        results_clinv.append(r_clinv_tp)
        print(f"  Shapes: CT={ref_clinv_tp.shape}, PY={pred_clinv_tp.shape}")
        print(f"  Mean Tanimoto: {r_clinv_tp['sim'].mean():.4f}  "
              f"Median: {float(np.median(r_clinv_tp['sim'])):.4f}")
        print(f"  Overall bit accuracy: {100*r_clinv_tp['overall_acc']:.4f}%")
        pd.DataFrame({"tanimoto_sim": r_clinv_tp["sim"]}).to_csv(
            args.out / "clinventory_toxprint_per_compound_tanimoto.csv", index=False
        )
    else:
        print(f"\n[skip] CLinventory ToxPrint NPZ not found: {CLINVENTORY_TOXPRINT_NPZ}")

    # ------------------------------------------------------------------
    # Dataset 3 — CLinventory, TxP_PFAS v1
    # ------------------------------------------------------------------
    if CLINVENTORY_TXPPFAS_NPZ.exists():
        print(f"\nLoading CLinventory TxP_PFAS cache: {CLINVENTORY_TXPPFAS_NPZ.name}")
        pred_clinv_pf, ref_clinv_pf, bn_clinv_pf = load_clinventory_npz(CLINVENTORY_TXPPFAS_NPZ)
        r_clinv_pf = analyse_pair(
            ref_clinv_pf, pred_clinv_pf,
            bn_clinv_pf, "CLinventory — TxP_PFAS v1", args.bad_threshold,
        )
        results_clinv.append(r_clinv_pf)
        print(f"  Shapes: CT={ref_clinv_pf.shape}, PY={pred_clinv_pf.shape}")
        print(f"  Mean Tanimoto: {r_clinv_pf['sim'].mean():.4f}  "
              f"Median: {float(np.median(r_clinv_pf['sim'])):.4f}")
        print(f"  Overall bit accuracy: {100*r_clinv_pf['overall_acc']:.4f}%")
        pd.DataFrame({"tanimoto_sim": r_clinv_pf["sim"]}).to_csv(
            args.out / "clinventory_txppfas_per_compound_tanimoto.csv", index=False
        )
    else:
        print(f"\n[skip] CLinventory TxP_PFAS NPZ not found: {CLINVENTORY_TXPPFAS_NPZ}")

    # ------------------------------------------------------------------
    # Markdown report
    # ------------------------------------------------------------------
    n_tc = r_tc["n_compounds"]
    bit_acc_tc = r_tc["bit_acc"]
    n_bits_tc = r_tc["n_bits"]
    n_perfect_tc = int((bit_acc_tc == 1.0).sum())
    n_above99_tc = int((bit_acc_tc >= 0.99).sum())
    n_above95_tc = int((bit_acc_tc >= 0.95).sum())
    n_below_tc   = int((bit_acc_tc < args.bad_threshold).sum())

    report = [
        "# Fingerprint Fidelity: pyCSRML vs ChemoTyper\n",
        "---\n",
        "## ToxCast — ToxPrint v2\n",
        f"**Dataset**: ToxCast ({n_tc} compounds)  ",
        f"**Fingerprint**: ToxPrint v2 ({n_bits_tc} bits)  \n",
        "### Per-compound Tanimoto similarity\n",
        "| Statistic | Value |",
        "|-----------|-------|",
    ]
    report += _tanimoto_md_block(r_tc)
    report += [
        "",
        "### Overall bit-level accuracy\n",
        f"Overall concordance: **{100*r_tc['overall_acc']:.4f}%** "
        f"({total_agree_tc:,d} / {n_tc * n_bits_tc:,d} bit positions agree)\n",
        "### Per-bit accuracy summary\n",
        "| Category | Count | Fraction |",
        "|----------|-------|---------|",
        f"| Perfect (acc = 1.00) | {n_perfect_tc} | {100*n_perfect_tc/n_bits_tc:.1f}% |",
        f"| acc ≥ 0.99 | {n_above99_tc} | {100*n_above99_tc/n_bits_tc:.1f}% |",
        f"| acc ≥ 0.95 | {n_above95_tc} | {100*n_above95_tc/n_bits_tc:.1f}% |",
        f"| acc < {args.bad_threshold:.2f} | {n_below_tc} | {100*n_below_tc/n_bits_tc:.1f}% |",
        "",
        f"See `{out_bits_tc.name}` for the full per-bit table.\n",
    ]
    thresh_bad = min(args.bad_threshold, 0.99)
    bad_bits_tc = bit_df_tc[bit_df_tc["accuracy"] < thresh_bad].head(30)
    if not bad_bits_tc.empty:
        report += [
            f"### Bits with accuracy < {thresh_bad:.0%}\n",
            "| Bit name | Accuracy | CT prevalence | PY prevalence |",
            "|----------|----------|--------------|--------------|",
        ]
        for _, row in bad_bits_tc.iterrows():
            report.append(
                f"| `{row.bit_name}` | {row.accuracy:.4f} "
                f"| {row.ct_prevalence:.4f} | {row.py_prevalence:.4f} |"
            )

    # CLinventory sections
    for r in results_clinv:
        n_c = r["n_compounds"]
        n_b = r["n_bits"]
        ba  = r["bit_acc"]
        n_p99 = int((ba >= 0.99).sum())
        report += [
            f"\n---\n",
            f"## {r['label']}\n",
            f"**Dataset**: {n_c:,d} compounds  ",
            f"**Fingerprint**: {n_b} bits  \n",
            "### Per-compound Tanimoto similarity\n",
            "| Statistic | Value |",
            "|-----------|-------|",
        ]
        report += _tanimoto_md_block(r)
        tot_ag = int((r["bit_acc"] > 0).sum())  # placeholder; compute properly
        tot_agree_c = int(np.sum(
            (r["sim"] > 0).astype(int)  # real total agree computed below
        ))
        # recompute total agree properly
        # (sim alone doesn't give it; use overall_acc * n)
        total_agree_c = int(round(r["overall_acc"] * n_c * n_b))
        report += [
            "",
            "### Overall bit-level accuracy\n",
            f"Overall concordance: **{100*r['overall_acc']:.4f}%** "
            f"({total_agree_c:,d} / {n_c * n_b:,d} bit positions agree)\n",
            "### Per-bit accuracy summary\n",
            "| Category | Count | Fraction |",
            "|----------|-------|---------|",
            f"| acc ≥ 0.99 | {n_p99} | {100*n_p99/n_b:.1f}% |",
        ]

    out_md = args.out / "fingerprint_fidelity_report.md"
    out_md.write_text("\n".join(report), encoding="utf-8")
    print(f"\nMarkdown report written to {out_md}")


if __name__ == "__main__":
    main()
