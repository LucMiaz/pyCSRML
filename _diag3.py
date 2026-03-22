"""Detailed bond-type analysis for FN molecules."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')
from rdkit import Chem

def analyze_mol(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  CANNOT PARSE: {smi}")
        return
    mol_h = Chem.AddHs(mol)
    print(f"  SMILES: {smi}")
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        btype = bond.GetBondTypeAsDouble()
        barom = bond.GetIsAromatic()
        print(f"    bond {a1.GetIdx()}({a1.GetSymbol()},arom={a1.GetIsAromatic()}) "
              f"- {a2.GetIdx()}({a2.GetSymbol()},arom={a2.GetIsAromatic()}) "
              f"order={btype} arom={barom}")

# Test C=N imine matching with different SMARTS
print("=== C=N imine SMARTS tests ===")
fns = [
    "FC(F)(F)C1=NC(=N)NN1",
    "FC(F)(F)C1=NNC(=O)O1",
    "O=C1C=C(NC(=N1)C(F)(F)F)C(F)(F)F",
    "FC(F)(F)C=1C=C(NC(=S)N=1)C(F)(F)F",
]
smarts_tests = [
    ("[#9]-[#6;A]-[#6;A]=[#7;A;H0,H1]", "current"),
    ("[#9]-[#6;A]-[#6]=[#7;H0,H1]", "no-A on C=N"),
    ("[#9]-[#6;A]-[#6]=[#7;!a;H0,H1]", "no-A on C=N, N not aromatic"),
    ("[#9]-[#6;A]-[#6;A]=[#7;H0,H1]", "no-A on N only"),
]

for smi in fns:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"Can't parse {smi}")
        continue
    mol_h = Chem.AddHs(mol)
    results = []
    for smarts, label in smarts_tests:
        q = Chem.MolFromSmarts(smarts)
        m = mol_h.HasSubstructMatch(q)
        results.append(f"{label}={m}")
    print(f"  {smi}")
    print(f"    {', '.join(results)}")

print()
# Check bond types for a ring-containing FN
print("=== bond analysis for ring FN ===")
analyze_mol("FC(F)(F)C1=NC(=N)NN1")
analyze_mol("FC(F)(F)C1=NNC(=O)O1")

print()
# Test alkeneLinear_F with no-A constraint
print("=== alkeneLinear_F SMARTS tests ===")
fns_alkene = [
    "FC1=C(F)C1=O",
    "FC=1C=C(F)C(=O)NC=1F",
    "O=C1C(F)=C(F)NC(F)=C1F",
]
smarts_alkene_tests = [
    ("[#9]-[#6;A]=[#6;A]", "current"),
    ("[#9]-[#6]=[#6]", "no-A"),
    ("[#9]-[#6]=[#6;!a]", "C2 not aromatic"),
]
for smi in fns_alkene:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"Can't parse {smi}")
        continue
    mol_h = Chem.AddHs(mol)
    results = []
    for smarts, label in smarts_alkene_tests:
        q = Chem.MolFromSmarts(smarts)
        m = mol_h.HasSubstructMatch(q)
        results.append(f"{label}={m}")
    print(f"  {smi}")
    print(f"    {', '.join(results)}")

print()
# Test fix for aromatic_FCc1c FPs
print("=== aromatic_FCc1c SMARTS tests ===")
fps_arom = [
    "FC(F)(F)C1=CC(=N)NN1",  # FP
    "FC(F)(F)C1=CC(=O)NN1",  # FP
    "FC(F)(F)C1=COC(=O)O1",  # FP
    "FC1=C(C=CNC1=O)C(F)(F)F",  # FP
    "FC(F)(F)c1ccccc1",  # should be TP
    "FC(F)(F)c1ccncc1",  # CF3 on pyridine, should probably be TP
]
smarts_arom_tests = [
    ("[#9]-[#6;A]-[#6;a]", "current"),
    ("[#9]-[#6;A]-[c;!$([c]~[#7,#8,#16;a])]", "no heteroatom neighbor"),
    ("[#9]-[#6;A]-[c;$([c]1[c][c][c][c][c]1)]", "in benzene ring"),
]
for smi in fps_arom:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"Can't parse {smi}")
        continue
    mol_h = Chem.AddHs(mol)
    results = []
    for smarts, label in smarts_arom_tests:
        q = Chem.MolFromSmarts(smarts)
        if q is None:
            results.append(f"{label}=INVALID_SMARTS")
            continue
        m = mol_h.HasSubstructMatch(q)
        results.append(f"{label}={m}")
    print(f"  {smi}")
    print(f"    {', '.join(results)}")
