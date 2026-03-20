"""
ToxPrint fingerprinter: compute binary chemotype fingerprints for molecules
using ToxPrint v2.0 (729 bits) or TxP_PFAS v1.0 (129 bits) definitions.

Usage::

    from pyToxPrint.fingerprinter import ToxPrintFingerprinter, PFASFingerprinter
    from rdkit import Chem

    fp = ToxPrintFingerprinter()          # loads bundled ToxPrint v2 XML
    mol = Chem.MolFromSmiles("c1ccccc1")
    arr, names = fp.fingerprint(mol)      # numpy bool array + list of bit names

    fp_pfas = PFASFingerprinter()         # loads bundled TxP_PFAS XML
    arr_pfas, names_pfas = fp_pfas.fingerprint(mol)

Pattern matching strategy
--------------------------
Each chemotype is defined by:
  1. A primary SMARTS pattern (substructureMatch molecule)
  2. Zero or more exception SMARTS patterns (substructureException molecules)

A fingerprint bit is set to 1 if:
  * The molecule contains a substructure match for the primary pattern, AND
  * The molecule does NOT contain a substructure match for any exception pattern
    (exception patterns are only applied when the exception molecule contains
     matchingQueryAtom cross-references to the main pattern; otherwise the
     exception acts as a global exclusion)

Note: The exception logic is a reasonable approximation; the original ChemoTyper
tool may produce slightly different results for edge cases.
"""
from __future__ import annotations

import json
import os
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Optional, Union

import numpy as np

_DATA_DIR = Path(__file__).parent / "data"
_TOXPRINT_XML = _DATA_DIR / "toxprint_V2.0_r711.xml"
_TXPPFAS_XML = _DATA_DIR / "TxP_PFAS_v1.0.4.xml"
_TOXPRINT_JSON = _DATA_DIR / "toxprint_V2.0_r711.json"
_TXPPFAS_JSON = _DATA_DIR / "TxP_PFAS_v1.0.4.json"


# ---------------------------------------------------------------------------
# Lazy RDKit import
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem
    from rdkit.Chem import MolFromSmarts

    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False

    class _FakeChem:  # type: ignore
        @staticmethod
        def MolFromSmarts(s):
            return None

        @staticmethod
        def MolFromSmiles(s):
            return None

    Chem = _FakeChem()  # type: ignore
    MolFromSmarts = Chem.MolFromSmarts


# ---------------------------------------------------------------------------
# Compiled pattern cache
# ---------------------------------------------------------------------------
_SMARTS_CACHE: dict[str, object] = {}


def _compile_smarts(smarts: str):
    """Return a cached compiled SMARTS query mol, or None on failure."""
    if smarts in _SMARTS_CACHE:
        return _SMARTS_CACHE[smarts]
    result = None
    if _HAS_RDKIT and smarts:
        try:
            result = Chem.MolFromSmarts(smarts)
        except Exception:
            result = None
    _SMARTS_CACHE[smarts] = result
    return result


# ---------------------------------------------------------------------------
# Load/build fingerprint spec
# ---------------------------------------------------------------------------

def _load_spec(json_path: Optional[Path], xml_path: Path) -> dict:
    """
    Load fingerprint spec from JSON cache if available and up-to-date,
    otherwise parse the XML and (optionally) save JSON.
    """
    if json_path and json_path.exists():
        # Check freshness
        if not xml_path.exists() or json_path.stat().st_mtime >= xml_path.stat().st_mtime:
            with open(json_path, encoding="utf-8") as f:
                return json.load(f)

    # Parse XML
    from pyToxPrint._csrml import parse_csrml_xml, ordered_bit_list

    parsed = parse_csrml_xml(str(xml_path))
    bit_order = ordered_bit_list(parsed)

    # Build flat spec: list of bit dicts
    bits = []
    for bit_id in bit_order:
        sg = parsed["subgraph_index"].get(bit_id)
        if sg is None:
            continue
        bits.append(
            {
                "id": sg["id"],
                "label": sg["label"],
                "smarts": sg["smarts"],
                "exception_smarts": sg["exception_smarts"],
            }
        )

    spec = {
        "id": parsed["id"],
        "title": parsed["title"],
        "csrml_version": parsed["csrml_version"],
        "n_bits": len(bits),
        "bits": bits,
    }

    # Save JSON for next time
    if json_path:
        try:
            json_path.parent.mkdir(parents=True, exist_ok=True)
            with open(json_path, "w", encoding="utf-8") as f:
                json.dump(spec, f, ensure_ascii=False, indent=2)
        except Exception:
            pass

    return spec


# ---------------------------------------------------------------------------
# Fingerprinter class
# ---------------------------------------------------------------------------

class Fingerprinter:
    """
    Compute binary chemotype fingerprints from a CSRML XML definition file.

    Parameters
    ----------
    xml_path : str or Path
        Path to the CSRML XML file.
    json_cache : str or Path, optional
        Path to a JSON cache file. If None, no caching is used.
    verbose : bool
        If True, warn about patterns that fail to compile.
    """

    def __init__(
        self,
        xml_path: Union[str, Path],
        json_cache: Optional[Union[str, Path]] = None,
        verbose: bool = False,
    ):
        if not _HAS_RDKIT:
            raise ImportError("RDKit is required for fingerprint computation.")

        self._xml_path = Path(xml_path)
        self._json_cache = Path(json_cache) if json_cache else None
        self._verbose = verbose
        self._spec = _load_spec(self._json_cache, self._xml_path)
        self._bits = self._spec["bits"]
        self._n_bits = len(self._bits)

        # Pre-compile SMARTS and warn on failures
        self._queries: list[tuple] = []  # (main_query, [exc_query, ...])
        self._valid: list[bool] = []
        n_failed = 0
        for bit in self._bits:
            main_q = _compile_smarts(bit["smarts"]) if bit.get("smarts") else None
            exc_qs = [
                q
                for s in bit.get("exception_smarts", [])
                for q in [_compile_smarts(s)]
                if q is not None
            ]
            ok = main_q is not None
            if not ok:
                n_failed += 1
            self._queries.append((main_q, exc_qs))
            self._valid.append(ok)

        if n_failed > 0 and self._verbose:
            warnings.warn(
                f"{n_failed}/{self._n_bits} patterns failed to compile and will "
                "always produce 0.",
                stacklevel=2,
            )

    # ------------------------------------------------------------------

    @property
    def n_bits(self) -> int:
        """Number of fingerprint bits."""
        return self._n_bits

    @property
    def bit_names(self) -> list[str]:
        """Ordered list of bit labels (one per bit)."""
        return [b["label"] for b in self._bits]

    @property
    def bit_ids(self) -> list[str]:
        """Ordered list of bit IDs (original subgraph IDs)."""
        return [b["id"] for b in self._bits]

    @property
    def title(self) -> str:
        return self._spec.get("title", "")

    # ------------------------------------------------------------------

    def fingerprint(self, mol) -> tuple[np.ndarray, list[str]]:
        """
        Compute the binary fingerprint for a molecule.

        Parameters
        ----------
        mol : rdkit.Chem.Mol
            An RDKit molecule object (must be pre-sanitized).

        Returns
        -------
        array : numpy.ndarray of dtype bool
            Binary fingerprint vector of length n_bits.
        names : list[str]
            Corresponding bit labels.
        """
        if mol is None:
            return np.zeros(self._n_bits, dtype=bool), self.bit_names

        arr = np.zeros(self._n_bits, dtype=bool)
        for i, (main_q, exc_qs) in enumerate(self._queries):
            if main_q is None:
                continue
            if mol.HasSubstructMatch(main_q):
                # Check exceptions
                hit = all(not mol.HasSubstructMatch(eq) for eq in exc_qs) if exc_qs else True
                arr[i] = hit
        return arr, self.bit_names

    def fingerprint_smiles(self, smiles: str) -> tuple[np.ndarray, list[str]]:
        """
        Compute fingerprint from a SMILES string.

        Returns all-zeros array if SMILES is invalid.
        """
        mol = Chem.MolFromSmiles(smiles)
        return self.fingerprint(mol)

    def fingerprint_batch(
        self,
        mols,
        smiles_list: Optional[list[str]] = None,
    ) -> np.ndarray:
        """
        Compute fingerprints for a list of molecules (or SMILES strings).

        Parameters
        ----------
        mols : iterable of rdkit.Chem.Mol or None
            If None is passed for a molecule, zeros are used.
        smiles_list : list of str, optional
            If provided, mols is ignored and this list of SMILES is used instead.

        Returns
        -------
        matrix : numpy.ndarray of shape (n_mols, n_bits), dtype bool
        """
        if smiles_list is not None:
            mols = [Chem.MolFromSmiles(s) for s in smiles_list]
        results = [self.fingerprint(m)[0] for m in mols]
        return np.vstack(results) if results else np.empty((0, self._n_bits), dtype=bool)


# ---------------------------------------------------------------------------
# Convenience subclasses
# ---------------------------------------------------------------------------

class ToxPrintFingerprinter(Fingerprinter):
    """
    ToxPrint v2.0 fingerprinter (729 chemotype bits).

    Parameters
    ----------
    xml_path : str or Path, optional
        Custom path to the ToxPrint XML. Defaults to the bundled file.
    """

    def __init__(
        self,
        xml_path: Optional[Union[str, Path]] = None,
        verbose: bool = False,
    ):
        xp = Path(xml_path) if xml_path else _TOXPRINT_XML
        jp = _TOXPRINT_JSON
        super().__init__(xml_path=xp, json_cache=jp, verbose=verbose)


class PFASFingerprinter(Fingerprinter):
    """
    TxP_PFAS v1.0 fingerprinter (129 PFAS-specific chemotype bits).

    Parameters
    ----------
    xml_path : str or Path, optional
        Custom path to the TxP_PFAS XML. Defaults to the bundled file.
    """

    def __init__(
        self,
        xml_path: Optional[Union[str, Path]] = None,
        verbose: bool = False,
    ):
        xp = Path(xml_path) if xml_path else _TXPPFAS_XML
        jp = _TXPPFAS_JSON
        super().__init__(xml_path=xp, json_cache=jp, verbose=verbose)
