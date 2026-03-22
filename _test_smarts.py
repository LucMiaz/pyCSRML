"""Test proposed SMARTS fixes against full Richard2023 dataset."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')
import numpy as np
import pandas as pd
from rdkit import Chem
import warnings; warnings.filterwarnings('ignore')

df_bits = pd.read_csv(r'tests/test_data/Richard2023_SI_TableS2.csv')
df_info = pd.read_csv(r'tests/test_data/Richard2023_SI_TableS5.csv', usecols=['DTXSID', 'SMILES'])
df = df_bits.merge(df_info, on='DTXSID')

# Preparse all mols
print("Parsing molecules...")
mol_h_list = []
for smi in df['SMILES']:
    mol = Chem.MolFromSmiles(smi)
    mol_h_list.append(Chem.AddHs(mol) if mol else None)
print(f"Parsed {sum(1 for m in mol_h_list if m is not None)} / {len(mol_h_list)}")

def count_fp_fn(smarts_str, ref_col):
    q = Chem.MolFromSmarts(smarts_str)
    if q is None:
        print(f"  INVALID SMARTS: {smarts_str}")
        return None
    pred = np.array([int(m.HasSubstructMatch(q)) if m else -1 for m in mol_h_list])
    ref = ref_col.astype(int).values
    valid = pred >= 0
    tp  = ((pred==1)&(ref==1)&valid).sum()
    fp_ = ((pred==1)&(ref==0)&valid).sum()
    fn  = ((pred==0)&(ref==1)&valid).sum()
    tn  = ((pred==0)&(ref==0)&valid).sum()
    prec = tp/(tp+fp_) if (tp+fp_) else float('nan')
    rec  = tp/(tp+fn)  if (tp+fn)  else float('nan')
    f1   = 2*prec*rec/(prec+rec) if (prec+rec) else float('nan')
    return tp, fp_, fn, tn, prec, rec, f1

print()
# === 1. CN_amine ===
print("=== pfas_bond:CN_amine_ter-N_generic_CF ===")
bit = 'pfas_bond:CN_amine_ter-N_generic_CF'
ref_col = df[bit]
tests = [
    ("[#6](-[#6;!$(*=[!#6;!#1])]-[#7](-[#6])-[#6])(-[#9])-*",
     "original (pre-fix)"),
    ("[#6](-[#6;!$(*=[!#6;!#1])]-[#7](-[#6;!$(*=[!#6;!#1])])-[#6;!$(*=[!#6;!#1])])(-[#9])-*",
     "after _csrml fix"),
    ("[#6](-[#6;!$(*=[!#6;!#1])]-[#7;!+](-[#6;!$(*=[!#6;!#1])])-[#6;!$(*=[!#6;!#1])])(-[#9])-*",
     "+!+ on N"),
]
for smar, label in tests:
    r = count_fp_fn(smar, ref_col)
    if r:
        tp, fp_, fn, tn, prec, rec, f1 = r
        print(f"  {label}: TP={tp} FP={fp_} FN={fn} prec={prec:.4f} rec={rec:.4f} f1={f1:.4f}")

print()
# === 2. aromatic_FCc1c ===
print("=== pfas_bond:aromatic_FCc1c ===")
bit = 'pfas_bond:aromatic_FCc1c'
ref_col = df[bit]
tests = [
    ("[#9]-[#6;A]-[#6;a]", "current"),
    ("[#9]-[#6;A]-[c;!$([c]~[#7,#8,#16;a])]", "no heteroatom neighbor"),
    ("[#9]-[#6;A]-[$([c]1:[c]:[c]:[c]:[c]:[c]:1)]", "benzene ring only"),
]
for smar, label in tests:
    r = count_fp_fn(smar, ref_col)
    if r:
        tp, fp_, fn, tn, prec, rec, f1 = r
        print(f"  {label}: TP={tp} FP={fp_} FN={fn} prec={prec:.4f} rec={rec:.4f} f1={f1:.4f}")

print()
# === 3. C=N imine ===
print("=== pfas_bond:C=N_imine_FCN ===")
bit = 'pfas_bond:C=N_imine_FCN'
ref_col = df[bit]
tests = [
    ("[#9]-[#6;A]-[#6;A]=[#7;A;H0,H1]", "current"),
    ("[#9]-[#6;A]-[#6]!-[#7;H0,H1]", "!- bond, no A"),
    ("[#9]-[#6;A]-[#6]!-[#7;!a;H0,H1]", "!- bond, N not aromatic"),
]
for smar, label in tests:
    r = count_fp_fn(smar, ref_col)
    if r:
        tp, fp_, fn, tn, prec, rec, f1 = r
        print(f"  {label}: TP={tp} FP={fp_} FN={fn} prec={prec:.4f} rec={rec:.4f} f1={f1:.4f}")

print()
# === 4. alkeneLinear_F ===
print("=== pfas_chain:alkeneLinear_mono-ene_ethylene_generic_F ===")
bit = 'pfas_chain:alkeneLinear_mono-ene_ethylene_generic_F'
ref_col = df[bit]
tests = [
    ("[#9]-[#6;A]=[#6;A]", "current"),
    ("[#9]-[#6]!-[#6]", "!- bond, no A"),
    ("[#9]-[#6;!a]!-[#6;!$([#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)]", "not on benzene"),
]
for smar, label in tests:
    r = count_fp_fn(smar, ref_col)
    if r:
        tp, fp_, fn, tn, prec, rec, f1 = r
        print(f"  {label}: TP={tp} FP={fp_} FN={fn} prec={prec:.4f} rec={rec:.4f} f1={f1:.4f}")
