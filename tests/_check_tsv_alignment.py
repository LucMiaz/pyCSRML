"""Diagnose aromaticity & SMARTS matching for problematic FN cases."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')
from rdkit import Chem

# Check how RDKit perceives the FN examples
cases = {
    'alkene FN': '[H]c1c(F)c(=O)n([H])c(=O)n1[H]',   # 5-fluorouracil
    'imine FN1': '[H]c1nc(=O)n([H])c(N([H])[H])c1F',   # 5-fluorocytosine
    'imine FN2': '[H]c1c([H])c([H])c(C([H])([H])[H])c(-n2c(C([H])([H])F)nc3c([H])c([H])c(N([H])[H])nc3=O)c1[H]',
    'benzene FP': '[H]OC([H])([H])C(=O)[C@@]1(O[H])[C@]([H])(C([H])([H])[H])C([H])([H])[C@@]2([H])[C@]([H])(C([H])([H])[H])C([H])([H])C([H])(C([H])([H])[H])[C@@]([H])([C@@]1([H])[H])[C@@]2([H])F',
}

smarts_tests = [
    ('[#9]-[#6;A]=[#6;A]', 'current alkene'),
    ('[#9]-[#6]=[#6]', 'no-A alkene'),
    ('[#9]-c:c', 'arene F-c:c'),
    ('[#9]-[#6;A]-[#6;A]=[#7;A;H0,H1]', 'current imine'),
    ('[#9]-[#6;A]-[#6]=[#7;H0,H1]', 'imine no-A on N'),
    ('[#9]-c1c(=O)n', 'uracil-like F-c(=O)n'),
]

for name, smi in cases.items():
    mol = Chem.MolFromSmiles(smi)
    print(f'\n{name}: {smi[:60]}...')
    if mol is None:
        print('  ERROR: Could not parse SMILES')
        continue
    # Show first 10 atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in (9, 6, 7):  # F, C, N only
            print(f'  atom {atom.GetIdx()}: {atom.GetSymbol()} aromatic={atom.GetIsAromatic()} '
                  f'charge={atom.GetFormalCharge()} Hs={atom.GetTotalNumHs()}')
    for pat_smi, pat_name in smarts_tests:
        pat = Chem.MolFromSmarts(pat_smi)
        matches = mol.GetSubstructMatches(pat) if pat else []
        print(f'  {pat_name}: {"MATCH " + str(len(matches)) if matches else "no match"}')
