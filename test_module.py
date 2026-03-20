"""
Quick smoke test for pyToxPrint.
"""
import sys
sys.path.insert(0, r'c:\Users\luc\git\ToxPrint')

from pyToxPrint import ToxPrintFingerprinter, PFASFingerprinter, from_fingerprinter
from rdkit import Chem

# --- Test PFAS fingerprinter ---
print("Loading PFASFingerprinter...", flush=True)
fp_pfas = PFASFingerprinter(verbose=True)
print(f"  n_bits = {fp_pfas.n_bits}", flush=True)

smiles_test = [
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFOA
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
    "FCCCF",                  # simple difluoro
    "c1ccccc1",               # benzene
    "CCO",                    # ethanol (negative control)
]

print("\nPFAS fingerprints:")
for smi in smiles_test:
    mol = Chem.MolFromSmiles(smi)
    arr, names = fp_pfas.fingerprint(mol)
    on_bits = int(arr.sum())
    on_names = [names[i] for i in range(len(arr)) if arr[i]][:5]
    print(f"  {smi[:40]:<40} → {on_bits:3d} bits  {on_names[:3]}")

# --- Test ToxPrint fingerprinter ---
print("\nLoading ToxPrintFingerprinter...", flush=True)
fp_tp = ToxPrintFingerprinter(verbose=True)
print(f"  n_bits = {fp_tp.n_bits}", flush=True)

for smi in smiles_test:
    mol = Chem.MolFromSmiles(smi)
    arr, names = fp_tp.fingerprint(mol)
    on_bits = int(arr.sum())
    print(f"  {smi[:40]:<40} → {on_bits:3d} bits")

# --- Test EmbeddingSet ---
print("\nBuilding EmbeddingSet...", flush=True)
eset = from_fingerprinter(
    fp_pfas,
    smiles_list=smiles_test,
    names=["PFOA", "PFOS", "4F-propane", "benzene", "ethanol"],
    dtxsids=["DTXSID1", "DTXSID2", "DTXSID3", "DTXSID4", "DTXSID5"],
)
print(f"  EmbeddingSet: {len(eset)} compounds, {eset.n_bits} bits")

# Tanimoto similarity
S = eset.similarity_matrix("tanimoto")
print(f"  PFOA-PFOS Tanimoto: {S[0,1]:.3f}")
print(f"  PFOA-benzene Tanimoto: {S[0,3]:.3f}")

# DataFrame
df = eset.to_dataframe()
print(f"  DataFrame shape: {df.shape}")

print("\nAll tests passed!", flush=True)
