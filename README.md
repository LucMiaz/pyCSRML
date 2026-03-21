# pyCSRML

[![CI](https://github.com/luc-miaz/pyCSRML/actions/workflows/ci.yml/badge.svg)](https://github.com/luc-miaz/pyCSRML/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pyCSRML)](https://pypi.org/project/pyCSRML/)
[![Python](https://img.shields.io/pypi/pyversions/pyCSRML)](https://pypi.org/project/pyCSRML/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Docs](https://readthedocs.org/projects/pycsrml/badge/?version=latest)](https://pycsrml.readthedocs.io)

**pyCSRML** is a pure-Python implementation of the
[Chemical Subgraph Representation Markup Language (CSRML)](https://www.molecular-networks.com/products/chemotyper).
It parses CSRML XML files, converts the subgraph patterns to SMARTS, and computes
binary chemotype fingerprints for molecules using [RDKit](https://www.rdkit.org/).

Two industry-standard fingerprint definitions are bundled:

| Fingerprint | Bits | Description |
|---|---|---|
| **ToxPrint v2.0** | 729 | General toxicologically relevant substructures |
| **TxP_PFAS v1.0** | 129 | Per- and polyfluoroalkyl substance (PFAS) chemotypes |

Concordance with the reference ChemoTyper tool: **99.4 %** overall accuracy across
14 710 compounds from Richard *et al.* (2023).

---

## Installation

```bash
pip install pyCSRML
```

For the optional analysis and visualisation features (pandas, matplotlib):

```bash
pip install "pyCSRML[analysis]"
```

For UMAP-based dimensionality reduction in `EmbeddingSet`:

```bash
pip install "pyCSRML[full]"
```

> **Conda users:** RDKit integrates best when installed via conda.
> ```bash
> conda install -c conda-forge rdkit
> pip install pyCSRML
> ```

---

## Quick start

### Single molecule

```python
from pyCSRML import PFASFingerprinter
from rdkit import Chem

fp = PFASFingerprinter()

mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O")  # PFOA
arr, names = fp.fingerprint(mol)

print(f"Bits set: {arr.sum()} / {fp.n_bits}")
on_bits = [names[i] for i in range(len(arr)) if arr[i]]
print(on_bits[:5])
```

### Multiple molecules with analysis

```python
from pyCSRML import PFASFingerprinter, from_fingerprinter

fp = PFASFingerprinter()

smiles = [
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",   # PFOA
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
    "FCCCF",    # simple difluoro
    "CCO",      # negative control
]

eset = from_fingerprinter(fp, smiles_list=smiles, names=["PFOA", "PFOS", "4F-propane", "EtOH"])
eset.plot(kind="heatmap")
```

### ToxPrint fingerprints (729 bits)

```python
from pyCSRML import ToxPrintFingerprinter
from rdkit import Chem

fp = ToxPrintFingerprinter()
mol = Chem.MolFromSmiles("c1ccccc1")
arr, names = fp.fingerprint(mol)
print(f"Benzene: {arr.sum()} bits set")
```

### Low-level CSRML parsing

```python
from pyCSRML._csrml import parse_csrml_xml, ordered_bit_list

data = parse_csrml_xml("path/to/my_fingerprints.xml")
bits = ordered_bit_list(data)
for bit in bits[:3]:
    print(bit["id"], bit["smarts"])
```

---

## API overview

| Class / function | Module | Description |
|---|---|---|
| `PFASFingerprinter` | `pyCSRML` | 129-bit TxP_PFAS fingerprinter |
| `ToxPrintFingerprinter` | `pyCSRML` | 729-bit ToxPrint v2 fingerprinter |
| `Fingerprinter` | `pyCSRML` | Base class; load any CSRML XML or JSON |
| `Embedding` | `pyCSRML` | Single-compound fingerprint container with metadata |
| `EmbeddingSet` | `pyCSRML` | Multi-compound container with heatmap / UMAP / clustering |
| `from_fingerprinter` | `pyCSRML` | Convenience factory: list of SMILES → `EmbeddingSet` |
| `parse_csrml_xml` | `pyCSRML._csrml` | Parse raw CSRML XML → Python dict |
| `ordered_bit_list` | `pyCSRML._csrml` | Return all bits in order from a parsed dict |

Full API reference: **[pycsrml.readthedocs.io](https://pycsrml.readthedocs.io)**

---

## CSRML features supported

| Feature | Status |
|---|---|
| `substructureMatch` → SMARTS | ✅ Full |
| `substructureException` (global) | ✅ Full |
| `matchingQueryAtom` → `[!$(...)]` folding | ✅ Full |
| `combineAtomFeatures` (OR-of-AND trees) | ✅ Full |
| `atomList` with `negate="true"` | ✅ Full |
| `attachedHydrogenCount` ranges | ✅ Full |
| `ringCountAtom` / `ringAtom` / `chainAtom` | ✅ Full |
| Pseudo-elements G, Z, Q, X | ✅ Full |
| `mustMatch` / `mustNotMatch` (test cases) | parsed, not used for matching |

---

## Development

```bash
git clone https://github.com/luc-miaz/pyCSRML
cd pyCSRML
pip install -e ".[dev]"

# Run tests (fast)
pytest -m "not slow"

# Run concordance test (~45 s)
pytest tests/test_chemotyper_concordance.py -v -s

# Pylint
pylint pyCSRML/
```

---

## Citation

If you use pyCSRML in academic work, please cite the original ToxPrint / ChemoTyper paper and the TxP_PFAS reference:

- Yang, C., *et al.* (2015). New publicly available chemical query language, CSRML, to support chemotype representations for application to data mining and modelling. *J. Chem. Inf. Model.* **55**, 510–528.
- Richard, A.M., *et al.* (2023). ToxPrint chemotypes and ChemoTyper portal. *Chem. Res. Toxicol.* **36**, 488–510.

---

## License

MIT © 2026 Luc Miaz

