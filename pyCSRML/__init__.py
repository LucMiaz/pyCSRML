"""
pyCSRML — Python implementation of CSRML chemotype fingerprints.

ToxPrint v2.0 (729 bits) and TxP_PFAS v1.0.4 (129 bits) definitions are bundled.

Quick start::

    from pyCSRML import ToxPrintFingerprinter, PFASFingerprinter, EmbeddingSet, from_fingerprinter
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

from pyCSRML.fingerprinter import (
    Fingerprinter,
    ToxPrintFingerprinter,
    PFASFingerprinter,
)
from pyCSRML.embedding import (
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
