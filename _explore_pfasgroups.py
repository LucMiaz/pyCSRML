"""Explore how PFASGroups handles sulfonyl ester fragmentation."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\PFASGroups')
import rdkit.RDLogger as RL
RL.DisableLog('rdApp.*')

from PFASGroups import parse_smiles

smiles_list = [
    'FC(F)(F)S(=O)(=O)OC',              # trifluoromethylsulfonyl methyl ester
    'OC(=O)C(F)(F)C(F)(F)F',            # PFBA (carboxylic acid)
    'OC(=O)C(F)(F)C(F)(F)C(F)(F)F',    # PFPeA (carboxylic acid)
    'FC(F)(F)S(=O)(=O)OCC',             # triflate ester
]

for smi in smiles_list:
    mols = parse_smiles([smi])
    m = mols[0]
    print(f'SMILES: {smi}')
    for comp in m.components:
        print(f'  Component: group_type={comp.group_type}, smiles={comp.smiles}')
    print()
