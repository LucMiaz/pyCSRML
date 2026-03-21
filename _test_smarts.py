"""Test SMARTS parsing for X3,#6 issue."""
from rdkit import Chem
import rdkit.RDLogger as RL
RL.DisableLog('rdApp.*')

pfhpa = Chem.MolFromSmiles('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F')
# COOH carbon is atom 1 (X3, aliphatic, no H)

tests = [
    ('[#6;A;X3]', 'baseline - should match atom 1'),
    ('[#6;A;X3,#6;A]', 'X3 OR A-carbon'),
    ('[#6;A,#6;A;X3]', 'reversed order'),
    ('[#6;A;X2,#6;A;X3]', 'X2 OR X3'),
    ('[#6;A;X3,#8;A]', 'X3-C OR A-O - should match atom 1'),
    ('[#6;A;X3,#6;R1]', 'X3 OR in-ring - BUG?'),
    ('[$([#6;A;X3]),$([#6;R1])]', 'recursive OR X3 or ring'),
    ('[$([#6;A;X3]),$([#6;A;!H0]),$([#6;R1]),$([#8;A]),$([#16;A]),$([#7;A]),$([#15;A]),$([#1]),$([#17]),$([#35])]',
     'full recursive OR'),
]

for smarts, comment in tests:
    q = Chem.MolFromSmarts(smarts)
    if q is None:
        print(f'INVALID SMARTS: {smarts[:60]} | {comment}')
    else:
        matches = pfhpa.GetSubstructMatches(q)
        atoms = [m[0] for m in matches]
        print(f'{smarts[:60]} | {comment}: atoms={atoms}')

print()
# Now test the full perF chain pattern with recursive terminal
core = '[#6](-[#9])(-[#9])(-[#9])-[#6](-[#9])(-[#9])-[#6](-[#9])(-[#9])-[#6](-[#9])(-[#9])-[#6](-[#9])(-[#9])-[#6](-[#9])(-[#9])'
term_recursive = '[$([#6;A;X3]),$([#6;A;!H0]),$([#6;R1]),$([#8;A]),$([#16;A]),$([#7;A]),$([#15;A]),$([#1]),$([#17]),$([#35])]'
full_recursive = core + '-' + term_recursive

q_full = Chem.MolFromSmarts(full_recursive)
print('Full pattern with recursive terminal valid:', q_full is not None)
if q_full:
    print('Matches PFHpA:', pfhpa.HasSubstructMatch(q_full))
