"""Diagnose FP/FN molecules for 4 low-metric bits."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')
import pandas as pd
from rdkit import Chem
from pyCSRML.fingerprinter import PFASFingerprinter

fp = PFASFingerprinter()
bit_names = fp.bit_names

df_bits = pd.read_csv(r'tests/test_data/Richard2023_SI_TableS2.csv')
df_info = pd.read_csv(r'tests/test_data/Richard2023_SI_TableS5.csv', usecols=['DTXSID', 'SMILES'])
df = df_bits.merge(df_info, on='DTXSID')

TARGET_BITS = [
    'pfas_bond:CN_amine_ter-N_generic_CF',
    'pfas_chain:alkeneLinear_mono-ene_ethylene_generic_F',
    'pfas_bond:aromatic_FCc1c',
    'pfas_bond:C=N_imine_FCN',
]

# subsets of target bits found in our fp
valid_targets = [b for b in TARGET_BITS if b in bit_names]
print(f"Analyzing {len(valid_targets)} bits\n")

for target in valid_targets:
    if target not in df.columns:
        print(f"\n{target}: NOT IN REFERENCE\n"); continue

    ref_col = df[target].values
    bit_idx = bit_names.index(target)

    fps_pred = []
    n_err = 0
    for smi in df['SMILES']:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                fps_pred.append(-1)
            else:
                arr, _ = fp.fingerprint(mol)
                fps_pred.append(int(arr[bit_idx]))
        except Exception:
            fps_pred.append(-1)
            n_err += 1

    import numpy as np
    pred = np.array(fps_pred)
    ref = ref_col.astype(int)
    valid = pred >= 0

    tp = ((pred == 1) & (ref == 1) & valid).sum()
    fp_ = ((pred == 1) & (ref == 0) & valid).sum()
    fn = ((pred == 0) & (ref == 1) & valid).sum()
    tn = ((pred == 0) & (ref == 0) & valid).sum()

    prec = tp/(tp+fp_) if (tp+fp_) else float('nan')
    rec  = tp/(tp+fn) if (tp+fn) else float('nan')
    print(f"\n=== {target} ===")
    print(f"  TP={tp}  FP={fp_}  FN={fn}  TN={tn}")
    print(f"  prec={prec:.4f}  rec={rec:.4f}")

    # Show FP molecules (pred=1, ref=0) up to 10
    if fp_ > 0:
        fp_mask = (pred == 1) & (ref == 0) & valid
        fp_rows = df[fp_mask][['DTXSID', 'SMILES']].head(10)
        print(f"  False Positives ({fp_} total, showing ≤10):")
        for _, row in fp_rows.iterrows():
            print(f"    {row['DTXSID']}  {row['SMILES']}")

    # Show FN molecules (pred=0, ref=1) up to 10
    if fn > 0:
        fn_mask = (pred == 0) & (ref == 1) & valid
        fn_rows = df[fn_mask][['DTXSID', 'SMILES']].head(10)
        print(f"  False Negatives ({fn} total, showing ≤10):")
        for _, row in fn_rows.iterrows():
            print(f"    {row['DTXSID']}  {row['SMILES']}")
