"""
create_size_benchmarks.py
=========================
Extract 5 molecule-size-stratified benchmark sets from the CLinventory CSV.

Heavy-atom bins
---------------
  bench_tiny    :   1 –  10 heavy atoms
  bench_small   :  11 –  20 heavy atoms
  bench_medium  :  21 –  35 heavy atoms
  bench_large   :  36 –  60 heavy atoms
  bench_xlarge  :  61+       heavy atoms

Each set is capped at 500 molecules (keeps ChemoTyper runs practical).

Outputs (all written to tests/test_data/size_benchmarks/)
----------------------------------------------------------
  bench_tiny.smiles, bench_small.smiles, ...   — one SMILES per line
  bench_metadata.csv                           — full metadata for each set
  chemotyper_timing_template.csv               — template for manual timing input

Usage
-----
  python scripts/create_size_benchmarks.py
  python scripts/create_size_benchmarks.py --max-per-set 200
"""

import argparse
import csv
import os
import re
import sys

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(REPO_ROOT, "tests", "test_data")
CLINVENTORY_CSV = os.path.join(DATA_DIR, "clinventory_molecules.csv")
OUT_DIR = os.path.join(DATA_DIR, "size_benchmarks")

# ---------------------------------------------------------------------------
# Size bins [(name, min_ha, max_ha)]
# ---------------------------------------------------------------------------
BINS = [
    ("bench_tiny",   1,  10),
    ("bench_small",  11, 20),
    ("bench_medium", 21, 35),
    ("bench_large",  36, 60),
    ("bench_xlarge", 61, None),
]

FP_TYPES = ["ToxPrint_V2", "TxP_PFAS_v1"]


def count_heavy_atoms(formula: str) -> int:
    """Sum all element counts in *formula*, excluding hydrogen.

    Examples
    --------
    >>> count_heavy_atoms("C9H10ClN")   # returns 9+1+1 = 11
    >>> count_heavy_atoms("HCl")        # returns 1
    >>> count_heavy_atoms("Na+")        # sodium ion — returns 1
    """
    total = 0
    for element, count_str in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if element == "H":
            continue
        total += int(count_str) if count_str else 1
    return total


def load_clinventory(csv_path: str) -> list[dict]:
    """Read CLinventory CSV and return list of row dicts with 'heavy_atoms' added."""
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            formula = row.get("formula", "")
            row["heavy_atoms"] = count_heavy_atoms(formula) if formula else 0
            rows.append(row)
    return rows


def assign_bins(rows: list[dict], max_per_set: int) -> dict[str, list[dict]]:
    """Partition rows into size bins, each capped at *max_per_set*."""
    sets: dict[str, list[dict]] = {name: [] for name, _, _ in BINS}
    for row in rows:
        ha = row["heavy_atoms"]
        for name, lo, hi in BINS:
            if ha >= lo and (hi is None or ha <= hi):
                if len(sets[name]) < max_per_set:
                    sets[name].append(row)
                break
    return sets


def write_outputs(sets: dict[str, list[dict]], out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    # --- per-set .smiles files ---
    for name, rows in sets.items():
        smiles_path = os.path.join(out_dir, f"{name}.smiles")
        with open(smiles_path, "w", encoding="utf-8") as fh:
            for row in rows:
                fh.write(row["smiles"].strip() + "\n")
        print(f"  {smiles_path}  ({len(rows)} molecules)")

    # --- bench_metadata.csv ---
    meta_path = os.path.join(out_dir, "bench_metadata.csv")
    fieldnames = ["set_name", "id", "inchikey", "formula", "heavy_atoms", "smiles"]
    with open(meta_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for name, rows in sets.items():
            for row in rows:
                writer.writerow({
                    "set_name":   name,
                    "id":         row.get("id", ""),
                    "inchikey":   row.get("inchikey", ""),
                    "formula":    row.get("formula", ""),
                    "heavy_atoms":row["heavy_atoms"],
                    "smiles":     row["smiles"].strip(),
                })
    print(f"  {meta_path}")

    # --- chemotyper_timing_template.csv ---
    template_path = os.path.join(out_dir, "chemotyper_timing_template.csv")
    t_fields = [
        "set_name", "fingerprint", "n_molecules",
        "rep1_seconds", "rep2_seconds", "rep3_seconds",
        "mean_seconds", "notes",
    ]
    with open(template_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=t_fields)
        writer.writeheader()
        for name, rows in sets.items():
            for fp in FP_TYPES:
                writer.writerow({
                    "set_name":    name,
                    "fingerprint": fp,
                    "n_molecules": len(rows),
                    "rep1_seconds": "",
                    "rep2_seconds": "",
                    "rep3_seconds": "",
                    "mean_seconds": "",
                    "notes":        "",
                })
    print(f"  {template_path}")


def print_summary(sets: dict[str, list[dict]]) -> None:
    print("\nMolecule-size distribution:")
    header = f"  {'Set':<16} {'HA range':<14} {'Count':>6}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for name, lo, hi in BINS:
        ha_range = f"{lo}–{hi}" if hi else f"{lo}+"
        print(f"  {name:<16} {ha_range:<14} {len(sets[name]):>6}")
    total = sum(len(v) for v in sets.values())
    print(f"  {'TOTAL':<16} {'':<14} {total:>6}")


def print_next_steps(out_dir: str) -> None:
    print("""
Next steps
----------
1. Time pyCSRML on these sets:
     python scripts/benchmark_pycsrml_timing.py

2. Run ChemoTyper on each .smiles file (both ToxPrint V2 and TxP_PFAS v1),
   export as TSV, place zips alongside the .smiles files:
     size_benchmarks/bench_tiny_toxprint.zip
     size_benchmarks/bench_tiny_txppfas.zip
     ...

3. Fill in the timing results from ChemoTyper:
     size_benchmarks/chemotyper_timing_template.csv

4. Run regression tests:
     pytest tests/test_benchmark_regression.py -v -m slow
""")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-per-set", type=int, default=500,
        help="Maximum molecules per size bin (default: 500).",
    )
    parser.add_argument(
        "--clinventory-csv", default=CLINVENTORY_CSV,
        help="Path to clinventory_molecules.csv.",
    )
    parser.add_argument(
        "--out-dir", default=OUT_DIR,
        help="Output directory for benchmark sets.",
    )
    args = parser.parse_args(argv)

    if not os.path.exists(args.clinventory_csv):
        print(f"ERROR: CLinventory CSV not found: {args.clinventory_csv}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading CLinventory from {args.clinventory_csv} …")
    rows = load_clinventory(args.clinventory_csv)
    print(f"  {len(rows)} molecules loaded.")

    sets = assign_bins(rows, args.max_per_set)
    print_summary(sets)

    print(f"\nWriting benchmark sets to {args.out_dir} …")
    write_outputs(sets, args.out_dir)

    print_next_steps(args.out_dir)


if __name__ == "__main__":
    main()
