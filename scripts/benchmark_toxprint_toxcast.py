#!/usr/bin/env python3
"""
Benchmark: ChemoTyper ToxPrint vs pyCSRML ToxPrint on ToxCast endpoints.

Compares the predictive power of 729-bit ToxPrint fingerprints computed by
two sources:
  * chemotyper  -- precomputed by the ChemoTyper reference tool
                   (PFASGroups/benchmark/test_data/toxprint_V2.tsv)
  * pycsrml     -- computed on-the-fly by this library (ToxPrintFingerprinter)

Cross-validation strategy
--------------------------
  Outer: StratifiedKFold(n_splits=5, shuffle=True)
  Inner: GridSearchCV(RandomForestClassifier, ..., cv=StratifiedKFold(3))
         RF param grid: n_estimators=[200,400], max_features=["sqrt",0.2]

This nested CV mirrors the strategy used in compare_fingerprints_toxcast.py
(PFASGroups benchmark) so results are directly comparable.

Usage
-----
    conda activate chem
    cd pyCSRML
    python scripts/benchmark_toxprint_toxcast.py              # full run
    python scripts/benchmark_toxprint_toxcast.py --fast       # skip inner CV
    python scripts/benchmark_toxprint_toxcast.py --dry-run    # data only
    python scripts/benchmark_toxprint_toxcast.py --db invitrodb_v4_3

Output
------
    scripts/data/toxprint_toxcast_results.csv
    scripts/data/toxprint_toxcast_summary.csv
    scripts/data/chemotyper_toxprint_cache_{hash}.npy
    scripts/data/pycsrml_toxprint_cache_{hash}.npy
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
import warnings

warnings.filterwarnings("ignore")

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent                               # pyCSRML/
DATA_OUT   = SCRIPT_DIR / "data"                             # output dir
PFASGROUPS = REPO_ROOT.parent / "PFASGroups" / "benchmark"  # sibling repo

TOXPRINT_TSV = PFASGROUPS / "test_data" / "toxprint_V2.tsv"
TXPSMILES    = PFASGROUPS / "test_data" / "toxcast_smiles.smi"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RANDOM_STATE = 42
OUTER_SPLITS = 5
INNER_SPLITS = 3
MIN_POS      = 30   # min. positives (and negatives) to run CV on an endpoint

LABEL_COLS = [
    "AR_antagonist",     "AR_agonist",
    "ERa_antagonist",    "ERa_agonist",
    "AhR_agonist",       "Aromatase_antagonist",
    "TR_antagonist",     "DT40_genotoxicity",
    "MMP_ratio",
    "CYP2D6_antagonist", "CYP2C19_antagonist", "CYP2C9_antagonist",
    "CYP3A4_antagonist", "p53_ratio",           "Caspase3_HEPG2",
]

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--db", "-d", default="toxcast",
        help=(
            "Dataset name: 'toxcast' (default) or 'invitrodb_v4_3'. "
            "Selects {db}_labels_summary.csv from PFASGroups/benchmark/data/."
        ),
    )
    p.add_argument(
        "--fast", action="store_true",
        help="Skip inner GridSearch; use fixed RF(n_estimators=200, max_features='sqrt').",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="Load and align data, then exit without running CV.",
    )
    p.add_argument(
        "--no-cache", action="store_true",
        help="Ignore existing .npy caches and recompute fingerprints.",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Caching helpers (mtime-keyed)
# ---------------------------------------------------------------------------

def _cache_hash(paths: list[Path]) -> str:
    """Stable 12-char hash from the mtimes of the given source files."""
    parts = [f"{p}:{int(p.stat().st_mtime)}" for p in sorted(paths)]
    return hashlib.md5("|".join(parts).encode()).hexdigest()[:12]


def _try_load_cache(cache_path: Path, shape_n: int) -> np.ndarray | None:
    """Load .npy if it exists and row count matches; otherwise return None."""
    if cache_path.exists():
        arr = np.load(cache_path)
        if arr.shape[0] == shape_n:
            return arr
    return None


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_labels(db: str) -> pd.DataFrame:
    """Return the labels CSV as a DataFrame with smiles + endpoint columns."""
    csv_path = PFASGROUPS / "data" / f"{db}_labels_summary.csv"
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Labels file not found: {csv_path}\n"
            "Make sure the PFASGroups repo is at ../PFASGroups relative to pyCSRML."
        )
    df = pd.read_csv(csv_path)
    present = [c for c in LABEL_COLS if c in df.columns]
    missing = [c for c in LABEL_COLS if c not in df.columns]
    if missing:
        print(f"  Warning: {len(missing)} LABEL_COLS absent from {csv_path.name}: {missing}")
    df = df[["smiles"] + present].dropna(subset=["smiles"]).reset_index(drop=True)
    print(f"  Loaded {db}: {len(df)} compounds, {len(present)} endpoints")
    return df


def _ik(smi: str) -> str | None:
    """Compute InChIKey for a SMILES string; return None on parse failure."""
    mol = Chem.MolFromSmiles(smi)
    return Chem.inchi.MolToInchiKey(mol) if mol else None


def load_chemotyper_aligned(
    ref_smiles: list[str],
    cache_path: Path | None = None,
) -> np.ndarray:
    """
    Load ChemoTyper ToxPrint TSV and re-index rows to match ref_smiles
    via InChIKey lookup.

    Returns uint8 array of shape (len(ref_smiles), 729).
    Compounds absent from the TSV receive all-zero rows.
    """
    if cache_path is not None:
        cached = _try_load_cache(cache_path, len(ref_smiles))
        if cached is not None:
            print(f"  [cache hit] ChemoTyper ToxPrint: {cached.shape}")
            return cached

    print(f"  Loading {TOXPRINT_TSV.name} ...")
    df_tsv = pd.read_csv(TOXPRINT_TSV, sep="\t", index_col=0)
    X_raw  = df_tsv.values.astype(np.uint8)
    print(f"  TSV shape: {X_raw.shape}  "
          f"nonzero rows: {int((X_raw.sum(1) > 0).sum())}")

    tsv_smiles = [
        ln.strip()
        for ln in TXPSMILES.read_text().splitlines()
        if ln.strip()
    ]
    if len(tsv_smiles) != X_raw.shape[0]:
        raise ValueError(
            f"{TXPSMILES.name} has {len(tsv_smiles)} entries but "
            f"{TOXPRINT_TSV.name} has {X_raw.shape[0]} rows -- cannot align"
        )

    print("  Building InChIKey index for TSV rows ...")
    tsv_key_to_idx: dict[str, int] = {}
    for i, smi in enumerate(tsv_smiles):
        key = _ik(smi)
        if key is not None:
            tsv_key_to_idx[key] = i

    print(f"  Aligning {len(ref_smiles)} reference SMILES ...")
    X_aligned = np.zeros((len(ref_smiles), X_raw.shape[1]), dtype=np.uint8)
    matched = 0
    for j, smi in enumerate(ref_smiles):
        idx = tsv_key_to_idx.get(_ik(smi))
        if idx is not None:
            X_aligned[j] = X_raw[idx]
            matched += 1

    print(
        f"  Alignment: {matched}/{len(ref_smiles)} matched, "
        f"{len(ref_smiles) - matched} filled with zeros"
    )

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache_path, X_aligned)
        print(f"  [cache saved] {cache_path.name}")

    return X_aligned


def compute_pycsrml(
    smiles_list: list[str],
    cache_path: Path | None = None,
) -> np.ndarray:
    """
    Compute ToxPrint fingerprints via pyCSRML for all SMILES.

    Returns uint8 array of shape (n, 729); invalid SMILES -> all-zero row.
    """
    if cache_path is not None:
        cached = _try_load_cache(cache_path, len(smiles_list))
        if cached is not None:
            print(f"  [cache hit] pyCSRML ToxPrint: {cached.shape}")
            return cached

    sys.path.insert(0, str(REPO_ROOT))
    from pyCSRML import ToxPrintFingerprinter  # noqa: PLC0415

    fp = ToxPrintFingerprinter()
    print(f"  Computing pyCSRML ToxPrint for {len(smiles_list)} SMILES ...")
    t0 = time.time()
    X = fp.fingerprint_batch(None, smiles_list=smiles_list).astype(np.uint8)
    dt = time.time() - t0
    print(
        f"  Done: {dt:.1f}s  ({dt / len(smiles_list) * 1000:.1f} ms/mol)  "
        f"shape={X.shape}  nonzero rows={int((X.sum(1) > 0).sum())}"
    )

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache_path, X)
        print(f"  [cache saved] {cache_path.name}")

    return X


# ---------------------------------------------------------------------------
# CV machinery
# ---------------------------------------------------------------------------

def _make_rf_fixed():
    from sklearn.ensemble import RandomForestClassifier  # noqa: PLC0415
    return RandomForestClassifier(
        n_estimators=200,
        max_features="sqrt",
        class_weight="balanced",
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )


def _make_rf_grid():
    from sklearn.ensemble import RandomForestClassifier  # noqa: PLC0415
    from sklearn.model_selection import GridSearchCV, StratifiedKFold  # noqa: PLC0415
    base = RandomForestClassifier(
        class_weight="balanced",
        random_state=RANDOM_STATE,
        n_jobs=1,
    )
    return GridSearchCV(
        base,
        {"n_estimators": [200, 400], "max_features": ["sqrt", 0.2]},
        cv=StratifiedKFold(INNER_SPLITS, shuffle=True, random_state=RANDOM_STATE),
        scoring="roc_auc",
        n_jobs=-1,
        refit=True,
    )


def run_nested_cv(
    X: np.ndarray,
    y: np.ndarray,
    endpoint: str,
    fset: str,
    fast: bool,
) -> list[dict]:
    """Nested (or simple) CV for one (feature_set, endpoint) pair."""
    from sklearn.metrics import (  # noqa: PLC0415
        average_precision_score,
        balanced_accuracy_score,
        matthews_corrcoef,
        roc_auc_score,
    )
    from sklearn.model_selection import StratifiedKFold  # noqa: PLC0415

    outer_cv = StratifiedKFold(OUTER_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    rows: list[dict] = []

    for fold_i, (tr, te) in enumerate(outer_cv.split(X, y)):
        X_tr, X_te = X[tr], X[te]
        y_tr, y_te = y[tr], y[te]

        if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
            continue

        model = _make_rf_fixed() if fast else _make_rf_grid()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.fit(X_tr, y_tr)
            prob = model.predict_proba(X_te)[:, 1]
            pred = (prob >= 0.5).astype(int)

        try:
            auc = float(roc_auc_score(y_te, prob))
        except ValueError:
            auc = float("nan")
        try:
            ap = float(average_precision_score(y_te, prob))
        except ValueError:
            ap = float("nan")

        best_params = (
            str(model.best_params_) if hasattr(model, "best_params_") else "fixed"
        )
        rows.append({
            "feature_set": fset,
            "endpoint":    endpoint,
            "model":       "RandomForest",
            "fold":        fold_i,
            "n_total":     int(len(y)),
            "n_pos":       int(y.sum()),
            "pos_rate":    round(float(y.mean()), 4),
            "roc_auc":     round(auc, 4),
            "avg_prec":    round(ap, 4),
            "mcc":         round(float(matthews_corrcoef(y_te, pred)), 4),
            "bal_acc":     round(float(balanced_accuracy_score(y_te, pred)), 4),
            "best_params": best_params,
        })

    return rows


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary(results: pd.DataFrame) -> None:
    """Print mean AUC per (endpoint, feature_set) side by side."""
    pivot = (
        results.groupby(["endpoint", "feature_set"])["roc_auc"]
        .mean()
        .unstack("feature_set")
    )
    if "chemotyper" in pivot.columns and "pycsrml" in pivot.columns:
        pivot["\u0394(pycsrml-chemotyper)"] = pivot["pycsrml"] - pivot["chemotyper"]

    w_ep  = max(len(ep) for ep in pivot.index) + 2
    cols  = list(pivot.columns)
    w_col = max(max(len(c) for c in cols) + 2, 10)
    hdr   = f"{'Endpoint':<{w_ep}}" + "".join(f"{c:>{w_col}}" for c in cols)

    print("\n" + "=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))
    for ep, row in pivot.iterrows():
        line = f"{ep:<{w_ep}}"
        for c in cols:
            v = row[c]
            line += f"{v:>{w_col}.4f}" if not np.isnan(v) else f"{'--':>{w_col}}"
        print(line)
    print("=" * len(hdr))
    print()

    overall = results.groupby("feature_set")["roc_auc"].mean()
    print("Overall mean AUC-ROC:")
    for fset, val in overall.items():
        print(f"  {fset:<15}  {val:.4f}")
    print()


def _print_endpoint_stats(df: pd.DataFrame, label_cols: list[str]) -> None:
    """Print positives/negatives per endpoint."""
    print(f"\n  {'Endpoint':<30}  {'n_pos':>6}  {'n_neg':>6}  "
          f"{'pos_rate':>9}  {'runnable':>8}")
    print("  " + "-" * 65)
    for ep in label_cols:
        y     = df[ep].fillna(0).values.astype(int)
        n_pos = int(y.sum())
        n_neg = int((y == 0).sum())
        ok    = "yes" if (n_pos >= MIN_POS and n_neg >= MIN_POS) else "skip"
        print(f"  {ep:<30}  {n_pos:>6}  {n_neg:>6}  "
              f"{y.mean():>9.3f}  {ok:>8}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = _parse_args()
    DATA_OUT.mkdir(parents=True, exist_ok=True)

    # -- 1. Load labels -------------------------------------------------------
    print("\n[1/4] Loading ToxCast labels ...")
    df = load_labels(args.db)
    smiles     = df["smiles"].tolist()
    label_cols = [c for c in LABEL_COLS if c in df.columns]

    # -- 2. Fingerprints ------------------------------------------------------
    print("\n[2/4] Loading / computing fingerprints ...")
    h         = _cache_hash([TOXPRINT_TSV, TXPSMILES])
    ct_cache  = DATA_OUT / f"chemotyper_toxprint_cache_{h}.npy"
    mqg_cache = DATA_OUT / f"pycsrml_toxprint_cache_{h}.npy"

    if args.no_cache:
        ct_cache.unlink(missing_ok=True)
        mqg_cache.unlink(missing_ok=True)

    X_ct  = load_chemotyper_aligned(smiles, cache_path=ct_cache)
    X_mqg = compute_pycsrml(smiles, cache_path=mqg_cache)

    feature_sets = {"chemotyper": X_ct, "pycsrml": X_mqg}

    if args.dry_run:
        print("\n[dry-run] Stopping after fingerprint loading.")
        print(f"  ChemoTyper : {X_ct.shape}")
        print(f"  pyCSRML    : {X_mqg.shape}")
        _print_endpoint_stats(df, label_cols)
        return

    # -- 3. Nested CV ---------------------------------------------------------
    mode = "fast (fixed params)" if args.fast else "nested CV (inner GridSearch)"
    print(f"\n[3/4] Running {mode} -- "
          f"{len(label_cols)} endpoints x 2 feature sets ...")

    all_rows: list[dict] = []
    t_start = time.time()

    for ep in label_cols:
        y     = df[ep].fillna(0).values.astype(int)
        n_pos = int(y.sum())
        n_neg = int((y == 0).sum())

        if n_pos < MIN_POS or n_neg < MIN_POS:
            print(f"  [skip] {ep}: n_pos={n_pos}, n_neg={n_neg}")
            continue

        print(f"  {ep}  (n_pos={n_pos}, n_neg={n_neg})")

        for fset_name, X in feature_sets.items():
            t0   = time.time()
            rows = run_nested_cv(X, y, ep, fset_name, args.fast)
            dt   = time.time() - t0
            if rows:
                mean_auc = np.nanmean([r["roc_auc"] for r in rows])
                print(f"    {fset_name:<15} AUC={mean_auc:.4f}  ({dt:.1f}s)")
            all_rows.extend(rows)

    print(f"\n  Total CV time: {time.time() - t_start:.1f}s")

    if not all_rows:
        print("Warning: no CV results (all endpoints skipped).")
        return

    # -- 4. Save & print -------------------------------------------------------
    print("\n[4/4] Saving results ...")
    results = pd.DataFrame(all_rows)

    tag          = f"{args.db}{'_fast' if args.fast else ''}"
    results_path = DATA_OUT / f"toxprint_{tag}_results.csv"
    results.to_csv(results_path, index=False)
    print(f"  {results_path}")

    summary = (
        results.groupby(["feature_set", "endpoint"])[
            ["roc_auc", "avg_prec", "mcc", "bal_acc"]
        ]
        .agg(["mean", "std"])
        .round(4)
    )
    summary.columns = ["_".join(c) for c in summary.columns]
    summary_path = DATA_OUT / f"toxprint_{tag}_summary.csv"
    summary.reset_index().to_csv(summary_path, index=False)
    print(f"  {summary_path}")

    print_summary(results)


if __name__ == "__main__":
    main()
