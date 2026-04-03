# pyCSRML

[![CI](https://github.com/luc-miaz/pyCSRML/actions/workflows/ci.yml/badge.svg)](https://github.com/luc-miaz/pyCSRML/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pyCSRML)](https://pypi.org/project/pyCSRML/)
[![Python](https://img.shields.io/pypi/pyversions/pyCSRML)](https://pypi.org/project/pyCSRML/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Docs](https://readthedocs.org/projects/pycsrml/badge/?version=latest)](https://pycsrml.readthedocs.io)

**pyCSRML** is a pure-Python re-implementation of the
[Chemical Subgraph Representation Markup Language (CSRML)](https://www.molecular-networks.com/products/chemotyper).
It parses CSRML XML files, converts the subgraph patterns to SMARTS, and computes
binary chemotype fingerprints for molecules using [RDKit](https://www.rdkit.org/).

The module is not an exact replicate of the original CSRML (see [performance section](#performance)). the original software should be preferred.

The module was implemented from two fingerprints descriptions:

| Fingerprint | Bits | Description | Sourcde |
|---|---|---|----|
| **ToxPrint v2.0** | 729 | General toxicologically relevant substructures | Yang *et al.* 2015  |
| **TxP_PFAS v1.0** | 129 | Per- and polyfluoroalkyl substance (PFAS) chemotypes | Richard *et al.* 2023



## Performance

Accuracy is measured by comparing pyCSRML bit vectors against the reference
[ChemoTyper](https://www.molecular-networks.com/products/chemotyper) tool output.  
Run `pytest tests/test_chemotyper_concordance.py -v -s` to reproduce; the full
per-bit breakdown is written to `tests/concordance_report.md`.

| Dataset | Compounds | Fingerprint | Overall accuracy | Bits ≥ 90 % acc | Macro MCC | Macro Bal Acc | Macro ROC-AUC |
|---|---|---|---|---|---|---|---|
| Richard *et al.* 2023 (PFAS set) | 14 710 | TxP_PFAS v1 | **99.99 %** | 129 / 129 | **0.9971** | **0.9989** | **0.9989** |
| ToxCast (full) | 9 014 | ToxPrint v2 | **98.17 %** | 711 / 729 | **0.9155** | **0.9634** | **0.9634** |
| ToxCast (CF-containing subset) | 808 | TxP_PFAS v1 | **99.98 %** | 129 / 129 | **0.9905** | **0.9924** | **0.9924** |

> **Reading the table:** "CF-containing subset" means only the 808 ToxCast compounds
> for which ChemoTyper sets at least one TxP_PFAS bit — the meaningful subset for
> PFAS accuracy benchmarking. Full-dataset TxP_PFAS accuracy appears inflated (100 %)
> because the vast majority of compounds are all-zero for every PFAS bit.

### Known discrepancies

The 18 bits below 90 % accuracy in ToxPrint v2 are all in the metal / inorganic
chemotype groups; TxP_PFAS v1 has 4 bits below 100 % (all above 98.9 %).  
Root causes (see `tests/_check_tsv_alignment.py` and `tests/concordance_report.md`):

| Bit / category | Fingerprint | Accuracy | Direction | Root cause |
|---|---|---|---|---|
| `atom:element_noble_gas` | ToxPrint | 0.0 % | False positives | Noble-gas SMARTS approximated as `[*]` — matches every atom |
| `atom:element_metal_group_III`, `atom:element_metal_poor_metal`, etc. | ToxPrint | 0.1 – 5 % | False positives | Metal / metalloid element-group patterns use G/Q pseudo-elements that are approximated as `[*]`, causing widespread false positives |
| `ring:hetero_[6]_N_tetrazine_generic`, `ring:hetero_[6]_N_triazine_generic` | ToxPrint | 30 – 32 % | False positives | Nitrogen-count constraints in 6-membered heteroaromatic rings use atom-count SMARTS that over-match similar rings |
| `pfas_chain:alkeneLinear_mono-ene_ethylene_generic_F` | TxP_PFAS | 98.9 % | False negatives (recall 40 %) | RDKit perceives the C=C of tautomeric fluoropyrimidines (5-fluorouracil) as aromatic; the SMARTS `[#9]-[#6;A]=[#6;A]` requires aliphatic atoms and misses them |
| `pfas_bond:C=N_imine_FCN` | TxP_PFAS | 99.5 % | False negatives (recall 33 %) | Same aromaticity issue: the C=N bond in fluorinated heterocycles is perceived as aromatic by RDKit, so the aliphatic imine SMARTS does not match |
| `pfas_bond:aromatic_FCc1c` | TxP_PFAS | 99.5 % | Slight false positives (precision 97.2 %) | Aromatic F-C pattern slightly over-matches due to SMARTS approximation of the exception clause |

---

## Installation

The module needs RDKit installed. If necessary, start by installing a environment manager first (e.g. Conda/Mamba, like [Miniforge3](https://github.com/conda-forge/miniforge)) and creating an environment, e.g.:

```bash
mamba create -n rdkit pytho
mamba activate rdkit
mamba install -y -c rdkit rdkit
```

Then install pyCSRML via PyPI:

```bash
pip install pyCSRML
```


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

## Licence
<a href="https://github.com/LucMiaz/pyCSRML">pyCSRML</a> © 1999 by <a href="https://cogitopia.dev">Luc T. Miaz</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">

## Acknowledgments
This project is part of the [ZeroPM project](https://zeropm.eu/) (WP2) and has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 101036756. This work was developed at the [Department of Environmental Science](https://aces.su.se) at Stockholm University.<br />


<img alt="EU logo" src="https://zeropm.eu/wp-content/uploads/2021/12/flag_yellow_low.jpg" width=100/>     <a rel='zeropm_web' href="https://zeropm.eu/"/><img alt="zeropm logo" src="https://zeropm.eu/wp-content/uploads/2022/01/ZeroPM-logo.png" width=250 /></a><a rel='zeropm_web' href="https://su.se/"/><img alt="zeropm logo" src="https://eu01web.zoom.us/account/branding/p/5065401a-9915-4baa-9c16-665dcd743470.png" width=200 /></a>

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
    
