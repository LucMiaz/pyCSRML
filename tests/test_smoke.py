"""
Fast smoke tests for pyCSRML — run in CI without any external test data.

These tests verify:
  - The package imports correctly
  - Fingerprinter loads the bundled ToxPrint and TxP_PFAS definitions
  - fingerprint() returns the expected shapes and types
  - Known PFAS compounds set characteristic bits
  - Non-PFAS compounds set zero PFAS bits
"""
from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem

from pyCSRML import Fingerprinter, TOXPRINT_PATH, TXPPFAS_PATH


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pfas_fp():
    return Fingerprinter(TXPPFAS_PATH)


@pytest.fixture(scope="module")
def toxprint_fp():
    return Fingerprinter(TOXPRINT_PATH)


PFOA_SMILES = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
PFOS_SMILES = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O"
BENZENE_SMILES = "c1ccccc1"
ETHANOL_SMILES = "CCO"


# ---------------------------------------------------------------------------
# Fingerprinter construction
# ---------------------------------------------------------------------------

def test_txppfas_n_bits(pfas_fp):
    assert pfas_fp.n_bits == 129


def test_toxprint_n_bits(toxprint_fp):
    assert toxprint_fp.n_bits == 729


def test_pfas_bit_names(pfas_fp):
    assert len(pfas_fp.bit_names) == 129
    assert all(isinstance(n, str) for n in pfas_fp.bit_names)


# ---------------------------------------------------------------------------
# fingerprint() output shape and dtype
# ---------------------------------------------------------------------------

def test_fingerprint_returns_array_and_names(pfas_fp):
    mol = Chem.MolFromSmiles(ETHANOL_SMILES)
    arr, names = pfas_fp.fingerprint(mol)
    assert isinstance(arr, np.ndarray)
    assert arr.dtype == bool
    assert arr.shape == (129,)
    assert len(names) == 129


def test_fingerprint_none_molecule(pfas_fp):
    """None mol (invalid SMILES) should return an all-zero array without raising."""
    arr, names = pfas_fp.fingerprint(None)
    assert arr.shape == (129,)
    assert arr.sum() == 0


# ---------------------------------------------------------------------------
# Known PFAS compounds produce positive bits
# ---------------------------------------------------------------------------

def test_pfoa_has_pfas_bits(pfas_fp):
    mol = Chem.MolFromSmiles(PFOA_SMILES)
    arr, names = pfas_fp.fingerprint(mol)
    assert arr.sum() >= 4, f"PFOA should set several PFAS bits, got {arr.sum()}"


def test_pfoa_perF_linear_bit(pfas_fp):
    """PFOA (CF3-(CF2)6-COOH, 7-carbon perF chain) must set the C7+ bit."""
    mol = Chem.MolFromSmiles(PFOA_SMILES)
    arr, names = pfas_fp.fingerprint(mol)
    name_list = list(names)
    idx = name_list.index("pfas_chain:perF-linear_C7_plus")
    assert arr[idx], "PFOA must hit the C7+ perfluorinated chain bit"


def test_pfos_has_pfas_bits(pfas_fp):
    mol = Chem.MolFromSmiles(PFOS_SMILES)
    arr, _ = pfas_fp.fingerprint(mol)
    assert arr.sum() > 5


# ---------------------------------------------------------------------------
# Non-PFAS compounds produce zero PFAS bits
# ---------------------------------------------------------------------------

def test_ethanol_no_pfas_bits(pfas_fp):
    mol = Chem.MolFromSmiles(ETHANOL_SMILES)
    arr, _ = pfas_fp.fingerprint(mol)
    assert arr.sum() == 0, "Ethanol should not match any PFAS chemotype"


def test_benzene_no_pfas_bits(pfas_fp):
    mol = Chem.MolFromSmiles(BENZENE_SMILES)
    arr, _ = pfas_fp.fingerprint(mol)
    assert arr.sum() == 0, "Benzene should not match any PFAS chemotype"


# ---------------------------------------------------------------------------
# ToxPrint sanity checks
# ---------------------------------------------------------------------------

def test_toxprint_benzene_ring_bit(toxprint_fp):
    """Benzene must set at least one aromatic ring chemotype."""
    mol = Chem.MolFromSmiles(BENZENE_SMILES)
    arr, names = toxprint_fp.fingerprint(mol)
    assert arr.sum() > 0, "Benzene should set ToxPrint ring bits"