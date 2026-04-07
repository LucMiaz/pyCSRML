"""
test_benchmark_regression.py
============================
Regression tests for the 5 size-stratified CLinventory benchmark sets.

Two test types, both marked ``pytest.mark.slow``:

timing regression
    Reload the timing baseline written by
    ``python scripts/benchmark_pycsrml_timing.py``, re-run 3 reps, and assert
    that the new median is within 30 % of the stored baseline.
    Skips cleanly if no baseline JSON exists (first-run friendly).

accuracy regression
    Load ChemoTyper fingerprint zips placed by the user in
    ``tests/test_data/size_benchmarks/`` (e.g. ``bench_tiny_toxprint.zip``),
    compute concordance via helpers from test_chemotyper_concordance, and
    assert that overall accuracy has not dropped more than 0.1 percentage-point
    from the baseline stored in ``bench_accuracy_baseline.json``.
    Skips cleanly when the zip files have not yet been added.

Run all benchmark regression tests:
    pytest tests/test_benchmark_regression.py -v -m slow
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np
import pytest
from rdkit import Chem

pytestmark = pytest.mark.slow

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
DATA_DIR = os.path.join(TESTS_DIR, "test_data")
BENCH_DIR = os.path.join(DATA_DIR, "size_benchmarks")

TIMING_BASELINE_JSON  = os.path.join(BENCH_DIR, "pycsrml_timing_baseline.json")
ACCURACY_BASELINE_JSON = os.path.join(BENCH_DIR, "bench_accuracy_baseline.json")

SET_NAMES = [
    "bench_tiny",
    "bench_small",
    "bench_medium",
    "bench_large",
    "bench_xlarge",
]

FP_PARAMS = [
    ("ToxPrint_V2",  "toxprint"),
    ("TxP_PFAS_v1",  "txppfas"),
]

TIMING_TOLERANCE  = 1.30   # new_median <= baseline_median * 1.30
ACCURACY_TOLERANCE = 0.001  # new_acc    >= baseline_acc - 0.001

# ---------------------------------------------------------------------------
# Import concordance helpers from test_chemotyper_concordance
# ---------------------------------------------------------------------------
_conc_spec = importlib.util.spec_from_file_location(
    "test_chemotyper_concordance",
    os.path.join(TESTS_DIR, "test_chemotyper_concordance.py"),
)
_conc_mod = importlib.util.module_from_spec(_conc_spec)
_conc_spec.loader.exec_module(_conc_mod)

_compute_concordance_from_arrays = _conc_mod._compute_concordance_from_arrays


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _smiles_path(set_name: str) -> str:
    return os.path.join(BENCH_DIR, f"{set_name}.smiles")


def _zip_path(set_name: str, fp_short: str) -> str:
    return os.path.join(BENCH_DIR, f"{set_name}_{fp_short}.zip")


def _load_smiles(set_name: str) -> list[str]:
    path = _smiles_path(set_name)
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def _parse_mols(smiles_list: list[str]) -> list:
    return [m for m in (Chem.MolFromSmiles(s) for s in smiles_list) if m is not None]


def _load_timing_baseline() -> dict:
    if not os.path.exists(TIMING_BASELINE_JSON):
        return {}
    with open(TIMING_BASELINE_JSON, encoding="utf-8") as fh:
        return json.load(fh)


def _load_accuracy_baseline() -> dict:
    if not os.path.exists(ACCURACY_BASELINE_JSON):
        return {}
    with open(ACCURACY_BASELINE_JSON, encoding="utf-8") as fh:
        return json.load(fh)


def _save_accuracy_baseline(data: dict) -> None:
    os.makedirs(BENCH_DIR, exist_ok=True)
    with open(ACCURACY_BASELINE_JSON, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)


def _run_timed(fp_obj, mols: list, n_reps: int) -> float:
    """Return the median wall-clock time (seconds) over n_reps."""
    times = []
    for _ in range(n_reps):
        t0 = time.perf_counter()
        fp_obj.fingerprint_batch(mols)
        times.append(time.perf_counter() - t0)
    return sorted(times)[len(times) // 2]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def fp_toxprint():
    from pyCSRML import ToxPrintFingerprinter
    return ToxPrintFingerprinter()


@pytest.fixture(scope="module")
def fp_pfas():
    from pyCSRML import PFASFingerprinter
    return PFASFingerprinter()


def _get_fp(fp_name: str, fp_toxprint, fp_pfas):
    return fp_toxprint if fp_name == "ToxPrint_V2" else fp_pfas


# ---------------------------------------------------------------------------
# Timing regression tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("set_name", SET_NAMES)
@pytest.mark.parametrize("fp_name,_fp_short", FP_PARAMS)
def test_timing_does_not_regress(set_name, fp_name, _fp_short, fp_toxprint, fp_pfas):
    """Re-time fingerprint_batch and assert no >30 % degradation vs baseline.

    Skips if:
      - no baseline JSON (run benchmark_pycsrml_timing.py first)
      - the set's .smiles file does not exist (run create_size_benchmarks.py)
      - fewer than 2 valid molecules in the set
    """
    baseline = _load_timing_baseline()
    if not baseline:
        pytest.skip("No timing baseline found. Run: python scripts/benchmark_pycsrml_timing.py")

    if set_name not in baseline or fp_name not in baseline[set_name]:
        pytest.skip(f"No baseline entry for ({set_name}, {fp_name}).")

    smiles = _load_smiles(set_name)
    if not smiles:
        pytest.skip(f"{_smiles_path(set_name)} not found — run create_size_benchmarks.py.")

    mols = _parse_mols(smiles)
    if len(mols) < 2:
        pytest.skip(f"Too few valid molecules in {set_name} ({len(mols)}).")

    fp_obj = _get_fp(fp_name, fp_toxprint, fp_pfas)
    new_median = _run_timed(fp_obj, mols, n_reps=3)

    stored = baseline[set_name][fp_name]
    baseline_median = stored["median_s"]
    limit = baseline_median * TIMING_TOLERANCE

    assert new_median <= limit, (
        f"Timing regression for ({set_name}, {fp_name}): "
        f"new_median={new_median*1000:.1f} ms > "
        f"baseline_median={baseline_median*1000:.1f} ms × {TIMING_TOLERANCE} "
        f"= {limit*1000:.1f} ms."
    )


# ---------------------------------------------------------------------------
# Accuracy regression tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("set_name", SET_NAMES)
@pytest.mark.parametrize("fp_name,fp_short", FP_PARAMS)
def test_accuracy_does_not_regress(set_name, fp_name, fp_short, fp_toxprint, fp_pfas):
    """Compare pyCSRML accuracy vs ChemoTyper zip for each benchmark set.

    Skips if:
      - the ChemoTyper zip has not been placed in size_benchmarks/
      - the set's .smiles file does not exist

    First run with a new zip: accuracy is stored as the baseline (test passes).
    Subsequent runs: assert accuracy >= baseline - 0.001 (0.1 pp tolerance).
    """
    import pandas as pd

    zip_path = _zip_path(set_name, fp_short)
    if not os.path.exists(zip_path):
        pytest.skip(
            f"ChemoTyper zip not found: {zip_path}\n"
            f"Run ChemoTyper on {_smiles_path(set_name)} and place the result there."
        )

    smiles = _load_smiles(set_name)
    if not smiles:
        pytest.skip(f"{_smiles_path(set_name)} not found — run create_size_benchmarks.py.")

    ref_df = pd.read_csv(zip_path, sep="\t", index_col=0, compression="zip")
    fp_obj = _get_fp(fp_name, fp_toxprint, fp_pfas)

    pred, ref, bit_names, n_failed = _compute_concordance_from_arrays(
        fp_obj, smiles, ref_df
    )
    overall_acc = float((pred == ref).mean())

    baseline_key = f"{set_name}_{fp_name}"
    acc_baselines = _load_accuracy_baseline()

    if baseline_key not in acc_baselines:
        acc_baselines[baseline_key] = {
            "overall_acc":  overall_acc,
            "n_compounds":  int(pred.shape[0]),
            "n_failed":     n_failed,
        }
        _save_accuracy_baseline(acc_baselines)
        print(
            f"\n[baseline] Wrote accuracy baseline for '{baseline_key}': "
            f"overall_acc={overall_acc:.4f}"
        )
        return

    stored_acc = acc_baselines[baseline_key]["overall_acc"]
    assert overall_acc >= stored_acc - ACCURACY_TOLERANCE, (
        f"Accuracy regression for ({set_name}, {fp_name}): "
        f"new_acc={overall_acc:.4f} < "
        f"baseline_acc={stored_acc:.4f} - {ACCURACY_TOLERANCE} "
        f"= {stored_acc - ACCURACY_TOLERANCE:.4f}. "
        f"(SMILES failures: {n_failed})"
    )
