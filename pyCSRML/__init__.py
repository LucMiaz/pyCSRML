"""
pyToxPrint — Python wrappers for ToxPrint v2.0 and TxP_PFAS v1.0 chemotype fingerprints.

Quick start::

    from pyToxPrint import ToxPrintFingerprinter, PFASFingerprinter, EmbeddingSet, from_fingerprinter
    from rdkit import Chem

    # --- Single compound ---
    fp = PFASFingerprinter()
    mol = Chem.MolFromSmiles("FCCCF")
    arr, names = fp.fingerprint(mol)
    print(f"PFAS bits set: {arr.sum()}/{fp.n_bits}")

    # --- Multiple compounds with analysis ---
    smiles = ["FC(F)(F)C(F)(F)C(F)(F)C(=O)O", "FCCCF", "c1ccccc1"]
    eset = from_fingerprinter(fp, smiles_list=smiles, names=["PFOA-like", "4F-butane", "benzene"])
    eset.plot(kind="heatmap")
    eset.plot(kind="umap")
"""

from pyToxPrint.fingerprinter import (
    Fingerprinter,
    ToxPrintFingerprinter,
    PFASFingerprinter,
)
from pyToxPrint.embedding import (
    Embedding,
    EmbeddingSet,
    from_fingerprinter,
)

__all__ = [
    "Fingerprinter",
    "ToxPrintFingerprinter",
    "PFASFingerprinter",
    "Embedding",
    "EmbeddingSet",
    "from_fingerprinter",
]
