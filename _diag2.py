"""Diagnose specific problem molecules."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')
from rdkit import Chem
from rdkit.Chem import AllChem

def try_match(smi, smarts, label=""):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  Can't parse: {smi}")
        return
    mol_h = Chem.AddHs(mol)
    q = Chem.MolFromSmarts(smarts)
    match = mol_h.HasSubstructMatch(q)
    print(f"  {label}: {match}  | {smi}")

# === CN_amine FPs: quaternary ammonium ===
print("=== CN_amine FPs: quaternary N+ should NOT match ===")
smarts_old = "[#6](-[#6;!$(*=[!#6;!#1])]-[#7](-[#6])-[#6])(-[#9])-*"
smarts_new = "[#6](-[#6;!$(*=[!#6;!#1])]-[#7](-[#6;!$(*=[!#6;!#1])])-[#6;!$(*=[!#6;!#1])])(-[#9])-*"
smarts_nocharge = "[#6](-[#6;!$(*=[!#6;!#1])]-[#7;!+](-[#6;!$(*=[!#6;!#1])])-[#6;!$(*=[!#6;!#1])])(-[#9])-*"
fps = ["C[N+](C)(C)CC(F)(F)C(F)F", "[I-].C[N+](C)(C)CC(F)(F)C(F)F"]
for smi in fps:
    try_match(smi, smarts_new, "new (FP?)")
    try_match(smi, smarts_nocharge, "!+ fix")

# CN_amine FN
print("\n=== CN_amine FN: should match ===")
fn_smi = "FC(F)(F)C=1C=C(N(C)C2=CC(=O)C=CC2=1)C(F)(F)F"
try_match(fn_smi, smarts_new, "new (FN?)")
try_match(fn_smi, smarts_nocharge, "!+ fix")

# === alkeneLinear_F FNs ===
print("\n=== alkeneLinear_F FNs: should match [#9]-[#6;A]=[#6;A] ===")
smarts_alkene = "[#9]-[#6;A]=[#6;A]"
fns = [
    "FC1=C(F)C1=O",
    "FC=1C=C(F)C(=O)NC=1F",
    "O=C1C(F)=C(F)NC(F)=C1F",
    "FC1=C(C=CNC1=O)C(F)(F)F",
    "FC1=CC(=CNC1=O)C(F)(F)F",
    "S=C1C(F)=C(F)NC(F)=C1F",
    "FC(F)(F)C=1NC(=O)NC(=O)C=1F",
]
for smi in fns:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  Can't parse: {smi}")
        continue
    mol_h = Chem.AddHs(mol)
    q = Chem.MolFromSmarts(smarts_alkene)
    match = mol_h.HasSubstructMatch(q)
    # Inspect aromaticity of carbons
    are_aromatic = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            are_aromatic.append(atom.GetIdx())
    print(f"  match={match}, aromatic_C_idxs={are_aromatic}, smi={smi}")

# === aromatic_FCc1c FPs: rings that are NOT aromatic ===
print("\n=== aromatic_FCc1c FPs: non-aromatic rings matching [#9]-[#6;A]-[#6;a] ===")
smarts_arom = "[#9]-[#6;A]-[#6;a]"
fps_arom = ["FC(F)(F)C1=CC(=N)NN1", "FC(F)(F)C1=CC(=O)NN1", "FC(F)(F)C1=COC(=O)O1"]
for smi in fps_arom:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  Can't parse: {smi}")
        continue
    mol_h = Chem.AddHs(mol)
    q = Chem.MolFromSmarts(smarts_arom)
    match = mol_h.HasSubstructMatch(q)
    are_aromatic_c = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIsAromatic()]
    are_aromatic_any = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()]
    print(f"  match={match}, aromatic_all={are_aromatic_any}, smi={smi}")
    # Try what atoms are aromatic
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            print(f"    atom {atom.GetIdx()} ({atom.GetSymbol()}) is aromatic")

# === C=N_imine FNs ===
print("\n=== C=N_imine FNs: should match [#9]-[#6;A]-[#6;A]=[#7;A;H0,H1] ===")
smarts_imine = "[#9]-[#6;A]-[#6;A]=[#7;A;H0,H1]"
fns_imine = [
    "FC(F)(F)C1=NC(=N)NN1",
    "FC(F)(F)C1=NNC(=O)O1",
    "FC(F)(F)C1=NNC(=N)S1",
    "FC(F)(F)C1=NNC(=O)S1",
    "O=C1C=C(NC(=N1)C(F)(F)F)C(F)(F)F",
    "FC(F)(F)N1C(=S)NN=C1C(F)(F)F",
]
for smi in fns_imine:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  Can't parse: {smi}")
        continue
    mol_h = Chem.AddHs(mol)
    q = Chem.MolFromSmarts(smarts_imine)
    match = mol_h.HasSubstructMatch(q)
    are_aromatic = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()]
    # Check N atoms
    n_atoms = [(a.GetIdx(), a.GetIsAromatic(), a.GetTotalNumHs()) for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    print(f"  match={match}, N_atoms(idx,aro,H)={n_atoms}, smi={smi}")
