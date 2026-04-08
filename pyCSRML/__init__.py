"""
pyCSRML — Python implementation of CSRML chemotype fingerprints.

ToxPrint v2.0 (729 bits) and TxP_PFAS v1.0.4 (129 bits) definitions are bundled.

Quick start::

    from pyCSRML import Fingerprinter, TOXPRINT_PATH, TXPPFAS_PATH
    from rdkit import Chem

    # ToxPrint v2.0 (729 bits)
    fp = Fingerprinter(TOXPRINT_PATH)
    mol = Chem.MolFromSmiles("c1ccccc1")
    arr, names = fp.fingerprint(mol)
    print(f"Benzene: {arr.sum()} bits set / {fp.n_bits}")

    # TxP_PFAS v1.0.4 (129 bits)
    fp_pfas = Fingerprinter(TXPPFAS_PATH)
    matrix = fp_pfas.fingerprint_batch(mols)   # shape (n_mols, 129)
"""

from pyCSRML.fingerprinter import (
    Fingerprinter,
    TOXPRINT_PATH,
    TXPPFAS_PATH,
)

__all__ = [
    "Fingerprinter",
    "TOXPRINT_PATH",
    "TXPPFAS_PATH",
]
