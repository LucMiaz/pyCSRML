"""
benchmark_pycsrml_timing.py
===========================
Time pyCSRML fingerprinting on the 5 size-stratified CLinventory benchmark sets.

For each (set × fingerprinter) pair the script:
  - Parses SMILES with RDKit (once, before the timed loop)
  - Runs ``fingerprint_batch`` N times
  - Reports median, mean, min, max, and ms/mol

Results are saved to size_benchmarks/pycsrml_timing_baseline.json (unless
``--no-save`` is passed), which is later read by the regression test suite.

Usage
-----
  python scripts/benchmark_pycsrml_timing.py
  python scripts/benchmark_pycsrml_timing.py --reps 3 --no-save
"""

import argparse
import json
import os
import sys
import time

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, REPO_ROOT)

BENCH_DIR = os.path.join(REPO_ROOT, "tests", "test_data", "size_benchmarks")
BASELINE_JSON = os.path.join(BENCH_DIR, "pycsrml_timing_baseline.json")

# Ordered so that the results table matches the progression from tiny → xlarge
SET_ORDER = [
    "bench_tiny",
    "bench_small",
    "bench_medium",
    "bench_large",
    "bench_xlarge",
]


def discover_sets(bench_dir: str) -> dict[str, str]:
    """Return {set_name: smiles_path} for sets that exist on disk."""
    found = {}
    for name in SET_ORDER:
        path = os.path.join(bench_dir, f"{name}.smiles")
        if os.path.exists(path):
            found[name] = path
    return found


def parse_smiles(smiles_path: str):
    """Parse SMILES with RDKit, return (mols, n_failed)."""
    from rdkit import Chem

    mols, n_failed = [], 0
    with open(smiles_path, encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    for smi in lines:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)
        else:
            n_failed += 1
    return mols, n_failed


def time_fingerprinter(fp_obj, mols: list, n_reps: int) -> dict:
    """Run fingerprint_batch *n_reps* times and return timing stats (seconds)."""
    times = []
    for _ in range(n_reps):
        t0 = time.perf_counter()
        fp_obj.fingerprint_batch(mols)
        times.append(time.perf_counter() - t0)
    times_sorted = sorted(times)
    median = times_sorted[len(times_sorted) // 2]
    mean = sum(times) / len(times)
    n_mols = len(mols)
    return {
        "n_mols":    n_mols,
        "reps":      n_reps,
        "median_s":  median,
        "mean_s":    mean,
        "min_s":     times_sorted[0],
        "max_s":     times_sorted[-1],
        "median_ms_per_mol": (median / n_mols * 1000) if n_mols else 0.0,
    }


def print_table(results: dict) -> None:
    """Print an aligned summary table."""
    col_w = 18
    fps = ["ToxPrint_V2", "TxP_PFAS_v1"]
    header = f"{'Set':<16} {'Mols':>6}  " + "  ".join(
        f"{'[' + fp + '] median ms/mol':>{col_w}}" for fp in fps
    )
    print()
    print(header)
    print("-" * len(header))
    for set_name, fp_data in results.items():
        n_mols = next(iter(fp_data.values()))["n_mols"]
        cells = []
        for fp in fps:
            if fp in fp_data:
                v = fp_data[fp]["median_ms_per_mol"]
                cells.append(f"{v:{col_w}.4f}")
            else:
                cells.append(f"{'N/A':>{col_w}}")
        print(f"{set_name:<16} {n_mols:>6}  " + "  ".join(cells))
    print()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reps", type=int, default=5,
        help="Number of timed repetitions per (set × fingerprinter) (default: 5).",
    )
    parser.add_argument(
        "--no-save", action="store_true",
        help="Do not write pycsrml_timing_baseline.json.",
    )
    parser.add_argument(
        "--bench-dir", default=BENCH_DIR,
        help="Directory containing bench_*.smiles files.",
    )
    args = parser.parse_args(argv)

    sets = discover_sets(args.bench_dir)
    if not sets:
        print(
            f"ERROR: No bench_*.smiles files found in {args.bench_dir}.\n"
            "Run `python scripts/create_size_benchmarks.py` first.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Found {len(sets)} benchmark set(s): {', '.join(sets)}")
    print(f"Timing with {args.reps} repetitions per (set × fingerprinter).\n")

    # Import once (loads bundled JSON definitions — cost paid once)
    from pyCSRML import Fingerprinter, TOXPRINT_PATH, TXPPFAS_PATH

    print("Initialising fingerprinters (loads bundled JSON definitions)…")
    fp_map = {
        "ToxPrint_V2":  Fingerprinter(TOXPRINT_PATH),
        "TxP_PFAS_v1":  Fingerprinter(TXPPFAS_PATH),
    }
    print("  Done.\n")

    results: dict[str, dict] = {}
    for set_name, smiles_path in sets.items():
        mols, n_failed = parse_smiles(smiles_path)
        print(f"{set_name}: {len(mols)} valid SMILES ({n_failed} failures)")
        if not mols:
            print(f"  Skipping — no valid molecules.")
            continue
        results[set_name] = {}
        for fp_name, fp_obj in fp_map.items():
            stats = time_fingerprinter(fp_obj, mols, args.reps)
            results[set_name][fp_name] = stats
            print(
                f"  {fp_name}: median {stats['median_s']*1000:.1f} ms total  "
                f"({stats['median_ms_per_mol']:.4f} ms/mol)"
            )

    print_table(results)

    if not args.no_save:
        os.makedirs(args.bench_dir, exist_ok=True)
        with open(BASELINE_JSON, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=2)
        print(f"Timing baseline saved to {BASELINE_JSON}")
        print(
            "Run `pytest tests/test_benchmark_regression.py -v -m slow` "
            "to execute timing regression tests."
        )
    else:
        print("(--no-save: baseline not written)")


if __name__ == "__main__":
    main()
