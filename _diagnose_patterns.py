"""Diagnose the [#1] (explicit H) issue and perF-linear_cap pattern issue."""
import json
from pathlib import Path
import sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import rdkit.RDLogger as RDLogger
RDLogger.DisableLog('rdApp.*')

# Load patterns
data = json.loads(Path(r'C:\Users\luc\git\pyCSRML\pyCSRML\data\TxP_PFAS_v1.0.4.json').read_text())
bits = {b['label']: b for b in data['bits']}

def test_smarts(smarts_str, mol, label=""):
    """Test a SMARTS pattern on both original and AddHs mol."""
    q = Chem.MolFromSmarts(smarts_str)
    if q is None:
        print(f"  INVALID SMARTS: {smarts_str[:60]}")
        return
    mol_h = Chem.AddHs(mol)
    match_orig = mol.HasSubstructMatch(q)
    match_hs   = mol_h.HasSubstructMatch(q)
    print(f"  {label or 'pattern'}: orig={match_orig}, with_Hs={match_hs}")

# ----- Test 1: FT_n1_hetero on a simple fluorotelomer -----
print("=== Test 1: FT_n1_hetero ===")
# CF2-CH2-OH pattern (2:1 fluorotelomer alcohol, FTOH): FC(F)CCO
# FT_n1 means 1 CH2 group between CF2 and the heteroatom
# Structure: (F)(F)C-CH2-OH  where the CF2 is attached to a CF2CF fragment
# 2:1 FTOH = F5C2-CH2-OH = FC(F)(F)C(F)(F)CCO ... hmm
# Simpler: 2:1 FTOH = FC(F)CCO (2:1 ratio = 2 CF2 : 1 CH2)
smiles_list = [
    ("2:1 FT alcohol", "OCCC(F)(F)C(F)(F)F"),          # short fluorotelomer alcohol 
    ("4:2 FT alcohol", "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),  # 4:2 FTOH
    ("6:2 FTOH", "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("8:2 FTOH", "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
]
b = bits['pfas_chain:FT_n1_hetero']
for name, smi in smiles_list:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  Could not parse: {smi}")
        continue
    print(f"\n{name} ({smi}):")
    test_smarts(b['smarts'], mol, 'FT_n1_hetero primary')

# ----- Test 2: FT_n1_O on specific compounds -----
print("\n=== Test 2: FT_n1_O ===")
b = bits['pfas_chain:FT_n1_O']
for name, smi in smiles_list[:2]:
    mol = Chem.MolFromSmiles(smi)
    if mol is None: continue
    print(f"\n{name} ({smi}):")
    test_smarts(b['smarts'], mol, 'FT_n1_O primary')

# ----- Test 3: perF-linear_cap_C6_excl_mod -----
print("\n=== Test 3: perF-linear_cap_C6_excl_mod ===")
b = bits['pfas_chain:perF-linear_cap_C6_excl_mod']
perfluoro_smiles = [
    ("PFHpA (C7 PFCA)", "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),  # CF3-(CF2)5-COOH
    ("PFOA (C8 PFCA)", "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),  # CF3-(CF2)6-COOH
    ("PFHxA (C6 PFCA)", "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),  # CF3-(CF2)4-COOH
    ("PFHxS (C6 PFSA)", "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),  # CF3-(CF2)5-SO3H
]
print(f"  SMARTS: {b['smarts'][:120]}")
print(f"  (cont): {b['smarts'][120:] if len(b['smarts']) > 120 else ''}")
for name, smi in perfluoro_smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None: continue
    print(f"\n{name} ({smi}):")
    test_smarts(b['smarts'], mol, 'cap_C6_excl_mod')

# ----- Test 4: COH_alcohol_generic_OCCF -----
print("\n=== Test 4: COH_alcohol_generic_OCCF ===")
b = bits['pfas_bond:COH_alcohol_generic_OCCF']
print(f"  SMARTS: {b['smarts']}")
alcohol_smiles = [
    ("FTOH 2:1", "OCCC(F)(F)C(F)(F)F"),
    ("fluorinated alcohol", "OC(C(F)(F)F)CF"),
    ("trifluoroethanol", "OCC(F)(F)F"),
]
for name, smi in alcohol_smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None: continue
    print(f"\n{name} ({smi}):")
    test_smarts(b['smarts'], mol, 'COH_alcohol')

# ----- Test 5: polyF_cap_CHF2CF -----
print("\n=== Test 5: polyF_cap_CHF2CF ===")
b = bits['pfas_chain:polyF_cap_CHF2CF']
print(f"  SMARTS: {b['smarts']}")
chf_smiles = [
    ("CHF2-CF2 fragment", "FC(F)C(F)F"),  # CHF2CF2 (difluoromethyl cap)
    ("CHF2CF2CF2", "FC(F)C(F)(F)C(F)F"),
    ("CHF2CF2CF2CF2F", "FC(F)C(F)(F)C(F)(F)CF"),
]
for name, smi in chf_smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None: continue
    print(f"\n{name} ({smi}):")
    test_smarts(b['smarts'], mol, 'polyF_cap_CHF2CF')
