"""
Concordance test: compare pyToxPrint PFASFingerprinter output against
ChemoTyper reference fingerprints from Richard et al. (2023).

Reference files (tests/test_data/):
  Richard2023_SI_TableS2.csv  — DTXSID + all 129 TxP_PFAS bits (ChemoTyper)
  Richard2023_SI_TableS5.csv  — DTXSID + SMILES + compound metadata

Run with:
    pytest tests/test_chemotyper_concordance.py -v -s
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyToxPrint import PFASFingerprinter

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
S2_PATH = os.path.join(DATA_DIR, "Richard2023_SI_TableS2.csv")  # reference fp
S5_PATH = os.path.join(DATA_DIR, "Richard2023_SI_TableS5.csv")  # SMILES


# ---------------------------------------------------------------------------
# Fixtures (module-scoped so data is loaded once)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def fp():
    return PFASFingerprinter()


@pytest.fixture(scope="module")
def reference_df():
    """Merge reference fingerprints (S2) with SMILES (S5) on DTXSID."""
    s2 = pd.read_csv(S2_PATH, index_col="DTXSID")
    s5 = pd.read_csv(S5_PATH, index_col="DTXSID")[["SMILES"]]
    merged = s2.join(s5, how="inner")
    # Drop rows with missing SMILES
    before = len(merged)
    merged = merged.dropna(subset=["SMILES"])
    after = len(merged)
    if before != after:
        print(f"\n[fixture] Dropped {before - after} rows with missing SMILES")
    return merged


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _compute_concordance(fp_obj, reference_df):
    """
    Returns
    -------
    pred_matrix  : (n_valid, n_bits) bool array  — our predictions
    ref_matrix   : (n_valid, n_bits) int array   — ChemoTyper reference
    valid_dtxsids: list of DTXSIDs that were successfully parsed
    bit_names    : list of bit name strings (length n_bits)
    n_failed     : number of compounds whose SMILES could not be parsed
    """
    bit_names = fp_obj.bit_names

    # Check that all our bit names are present in the reference DataFrame
    missing = [b for b in bit_names if b not in reference_df.columns]
    if missing:
        pytest.fail(
            f"{len(missing)} bit name(s) from PFASFingerprinter are absent in "
            f"the reference CSV, e.g. {missing[:3]}"
        )

    smiles_list = reference_df["SMILES"].tolist()
    dtxsids = reference_df.index.tolist()

    # Parse SMILES, keep track of valid indices
    mols, valid_idx = [], []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is not None:
            mols.append(mol)
            valid_idx.append(i)

    n_failed = len(smiles_list) - len(valid_idx)

    # Batch fingerprint (our implementation)
    pred_matrix = fp_obj.fingerprint_batch(mols).astype(int)          # (n, 129)

    # Reference fingerprint aligned to our bit order
    ref_matrix = reference_df.iloc[valid_idx][bit_names].to_numpy(dtype=int)  # (n, 129)

    valid_dtxsids = [dtxsids[i] for i in valid_idx]
    return pred_matrix, ref_matrix, valid_dtxsids, bit_names, n_failed


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_bit_names_match_reference_columns(fp, reference_df):
    """Every PFASFingerprinter bit name must exist as a column in TableS2."""
    missing = [b for b in fp.bit_names if b not in reference_df.columns]
    assert not missing, (
        f"{len(missing)} bit name(s) not found in reference CSV: {missing}"
    )


def test_n_bits(fp):
    """PFASFingerprinter should produce exactly 129 bits."""
    assert len(fp.bit_names) == 129, f"Expected 129 bits, got {len(fp.bit_names)}"


def test_chemotyper_concordance(fp, reference_df, capsys):
    """
    Compare pyToxPrint fingerprints against ChemoTyper reference.

    Prints a full concordance report.  Asserts that:
      - overall bit accuracy >= 90 %
      - per-bit accuracy >= 70 % for at least 90 % of bits
    """
    pred, ref, dtxsids, bit_names, n_failed = _compute_concordance(fp, reference_df)

    n_compounds, n_bits = pred.shape
    match = pred == ref  # (n, 129) bool

    overall_acc = match.mean()
    bit_acc = match.mean(axis=0)           # (129,)
    compound_acc = match.mean(axis=1)      # (n,)

    # -- TP / FP / FN per bit --
    tp = ((pred == 1) & (ref == 1)).sum(axis=0)
    fp_count = ((pred == 1) & (ref == 0)).sum(axis=0)
    fn = ((pred == 0) & (ref == 1)).sum(axis=0)
    precision = np.where(tp + fp_count > 0, tp / (tp + fp_count), np.nan)
    recall    = np.where(tp + fn > 0,        tp / (tp + fn),        np.nan)
    f1        = np.where(precision + recall > 0,
                         2 * precision * recall / (precision + recall), np.nan)

    # -- Print report --
    with capsys.disabled():
        print(f"\n{'='*65}")
        print(f"  ChemoTyper concordance report — TxP_PFAS (Richard 2023)")
        print(f"{'='*65}")
        print(f"  Compounds processed : {n_compounds}  (SMILES failures: {n_failed})")
        print(f"  Bits compared       : {n_bits}")
        print(f"  Overall accuracy    : {overall_acc:.4f}  ({overall_acc*100:.2f} %)")
        print(f"  Bits ≥ 90% acc.     : {(bit_acc >= 0.90).sum()} / {n_bits}")
        print(f"  Bits ≥ 70% acc.     : {(bit_acc >= 0.70).sum()} / {n_bits}")
        print(f"  Cmpds with 100% acc : {(compound_acc == 1).sum()} / {n_compounds}")
        print(f"  Mean cmpd accuracy  : {compound_acc.mean():.4f}")
        print()

        # Per-bit summary sorted by accuracy
        print(f"  {'Bit name':<60} {'acc':>5}  {'prec':>5}  {'rec':>5}  {'f1':>5}")
        print(f"  {'-'*60}  {'-----':>5}  {'-----':>5}  {'-----':>5}  {'-----':>5}")
        order = np.argsort(bit_acc)
        for i in order:
            p = f"{precision[i]:.3f}" if not np.isnan(precision[i]) else "  N/A"
            r = f"{recall[i]:.3f}"    if not np.isnan(recall[i])    else "  N/A"
            f = f"{f1[i]:.3f}"        if not np.isnan(f1[i])        else "  N/A"
            flag = " ← WORST" if i == order[0] else ""
            print(f"  {bit_names[i]:<60} {bit_acc[i]:.3f}  {p:>5}  {r:>5}  {f:>5}{flag}")

        print(f"{'='*65}\n")

    # -- Assertions --
    assert overall_acc >= 0.90, (
        f"Overall bit accuracy {overall_acc:.4f} is below the 0.90 threshold. "
        "Check the per-bit table printed above."
    )
    pct_bits_above_70 = (bit_acc >= 0.70).mean()
    assert pct_bits_above_70 >= 0.90, (
        f"Only {pct_bits_above_70*100:.1f}% of bits have ≥70% accuracy "
        f"(threshold: 90% of bits)."
    )


def test_known_compounds_spot_check(fp, reference_df):
    """
    Spot-check PFOA and PFOS against the ChemoTyper reference from S2.

    Fetches SMILES from S5 (via reference_df) and expected bits from S2.
    Asserts:
      - key structural bits are correctly set
      - per-compound accuracy >= 95 %
    Note: a handful of bits (e.g. pfas_atom:element_metal_metalloid_CF) are
    known false positives caused by SMARTS approximations of the G/Q pseudo-
    elements; 100% exact match is not expected.
    """
    known = {
        "DTXSID8031865": {           # PFOA, CAS 335-67-1
            "name": "PFOA",
            "must_have": [
                "pfas_bond:C(=O)O_carboxylicAcid_generic_CF",
                "pfas_chain:perF-linear_C7_plus",
            ],
            "must_not_have": [
                "pfas_bond:S(=O)O_sulfonicAcid_acyclic_(chain)_SCF",
            ],
        },
        "DTXSID3031864": {           # PFOS, CAS 1763-23-1
            "name": "PFOS",
            "must_have": [
                "pfas_bond:S(=O)O_sulfonicAcid_acyclic_(chain)_SCF",
                "pfas_chain:perF-linear_C8_plus",
            ],
            "must_not_have": [
                "pfas_bond:C(=O)O_carboxylicAcid_generic_CF",
            ],
        },
    }

    fp_cols = [c for c in fp.bit_names if c in reference_df.columns]

    for dtxsid, spec in known.items():
        name = spec["name"]
        if dtxsid not in reference_df.index:
            pytest.skip(f"{dtxsid} ({name}) not found in reference data")

        row = reference_df.loc[dtxsid]
        mol = Chem.MolFromSmiles(str(row["SMILES"]))
        assert mol is not None, f"Could not parse SMILES for {name} ({dtxsid})"

        arr, bit_names = fp.fingerprint(mol)
        on_bits = {n for n, v in zip(bit_names, arr) if v}

        for b in spec["must_have"]:
            assert b in on_bits, f"{name}: required bit '{b}' not set"
        for b in spec["must_not_have"]:
            assert b not in on_bits, f"{name}: bit '{b}' should NOT be set"

        # Per-compound accuracy against reference
        pred = arr.astype(int)
        ref = row[fp_cols].to_numpy(dtype=int)
        acc = (pred == ref).mean()
        assert acc >= 0.95, (
            f"{name} ({dtxsid}): accuracy {acc:.3f} < 0.95. "
            f"Mismatches: "
            + "; ".join(
                f"{fp_cols[i]}(ref={int(ref[i])},pred={int(pred[i])})"
                for i in range(len(fp_cols))
                if pred[i] != ref[i]
            )
        )

    # Non-PFAS (benzene) should always have zero bits
    benzene = Chem.MolFromSmiles("c1ccccc1")
    arr_benz, _ = fp.fingerprint(benzene)
    assert not arr_benz.any(), "Benzene should have no TxP_PFAS bits"
