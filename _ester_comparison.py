"""Compare ester SMARTS and test the negate fix."""
import json, re
from pathlib import Path
from rdkit import Chem
import rdkit.RDLogger as RL
RL.DisableLog('rdApp.*')

data = json.loads(Path('pyCSRML/data/TxP_PFAS_v1.0.4.json').read_text())
bits = {b['label']: b for b in data['bits']}

labels = [
    'pfas_bond:C(=O)O_carboxylicEster_acyclic_C(=O)CF',
    'pfas_bond:C(=O)O_carboxylicEster_acyclic_OCCF',
    'pfas_bond:S(=O)O_sulfonicEster_acyclic_S-C_(chain)_SCF',
]

print("=== Current SMARTS ===")
for label in labels:
    b = bits[label]
    print(f'{label}:')
    print(f'  {b["smarts"]}')
    print()

# Test the sulfonyl ester on reference compounds
print("=== Sulfonyl ester on reference compounds (current SMARTS) ===")
sulfester = Chem.MolFromSmarts(bits[labels[2]]['smarts'])
smarts_fixed = bits[labels[2]]['smarts'].replace('[#6;R]', '[#6;!R]')
sulfester_fixed = Chem.MolFromSmarts(smarts_fixed)
print(f'Fixed SMARTS would be: {smarts_fixed}')

test_smiles = [
    'FC(F)(F)S(=O)(=O)OC',
    'FC(F)(F)S(=O)(=O)OCC',
    'O=S(=O)(OCCCF)C(F)(F)F',
    'FC(F)OS(=O)(=O)C(F)F',
    'FC(F)(F)S(=O)(=O)OC=C',
]
for smi in test_smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None: continue
    mol_h = Chem.AddHs(mol)
    orig = mol_h.HasSubstructMatch(sulfester)
    fixed = mol_h.HasSubstructMatch(sulfester_fixed)
    print(f'  {smi}: orig={orig}, fixed={fixed}')
