"""Verify JSON regeneration worked and fixes are correct."""
import json
from pathlib import Path
from rdkit import Chem
import rdkit.RDLogger as RL
RL.DisableLog('rdApp.*')

data = json.loads(Path(r'pyCSRML/data/TxP_PFAS_v1.0.4.json').read_text())
bits = {b['label']: b for b in data['bits']}

# Check the fixed patterns
for label in ['pfas_chain:perF-linear_cap_C6_excl_mod', 'pfas_chain:FT_n1_hetero']:
    b = bits[label]
    print(f'{label}:')
    smarts = b['smarts']
    print(f'  smarts = {smarts[:120]}')
    q = Chem.MolFromSmarts(smarts)
    print(f'  valid = {q is not None}')
    print()

# Test PFHpA on cap_C6
pfhpa = Chem.MolFromSmiles('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F')
q = Chem.MolFromSmarts(bits['pfas_chain:perF-linear_cap_C6_excl_mod']['smarts'])
pfhpa_h = Chem.AddHs(pfhpa)
print('PFHpA matches perF-linear_cap_C6_excl_mod (orig):', pfhpa.HasSubstructMatch(q))
print('PFHpA matches perF-linear_cap_C6_excl_mod (addHs):', pfhpa_h.HasSubstructMatch(q))

# Test FT_n1_hetero on CHF2-CF2-CH2-NH2
mol = Chem.MolFromSmiles('FC(F)C(F)(F)CN')
mol_h = Chem.AddHs(mol)
q2 = Chem.MolFromSmarts(bits['pfas_chain:FT_n1_hetero']['smarts'])
print('CHF2-CF2-CH2-NH2 matches FT_n1_hetero (orig):', mol.HasSubstructMatch(q2))
print('CHF2-CF2-CH2-NH2 matches FT_n1_hetero (addHs):', mol_h.HasSubstructMatch(q2))

# Also verify actual reference compounds for perF-linear_cap_C6_excl_mod
import pandas as pd
try:
    s2 = pd.read_csv(r'tests/test_data/Richard2023_SI_TableS2.csv', index_col='DTXSID')
    s5 = pd.read_csv(r'tests/test_data/Richard2023_SI_TableS5.csv', index_col='DTXSID')
    df = s2.join(s5[['SMILES']])
    pos = df[df['pfas_chain:perF-linear_cap_C6_excl_mod'] == 1].head(3)
    print()
    print('Reference positives for perF-linear_cap_C6_excl_mod:')
    for dtxsid, row in pos.iterrows():
        mol = Chem.MolFromSmiles(str(row['SMILES']))
        if mol:
            mol_h = Chem.AddHs(mol)
            matches = mol_h.HasSubstructMatch(q)
            print(f'  {dtxsid}: {row["SMILES"][:50]} -> {matches}')
except Exception as e:
    print(f'Could not load reference data: {e}')
