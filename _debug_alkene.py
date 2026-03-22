"""Debug alkeneLinear SMARTS generation."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')
from pyCSRML._csrml import _parse_molecule, _graph_to_smarts, _atom_to_smarts, NS
import xml.etree.ElementTree as ET

tree = ET.parse(r'pyCSRML/data/TxP_PFAS_v1.0.4.xml')
root = tree.getroot()

for sg in root.iter():
    if not sg.tag.endswith('subgraph'):
        continue
    lbl = ''
    for child in sg:
        if child.tag.endswith('label'):
            lbl = (child.text or '').strip()
            break
    if 'alkeneLinear' not in lbl:
        continue

    print('Found subgraph:', lbl)
    for mol_el in sg:
        if not mol_el.tag.endswith('molecule'):
            continue
        is_match = False
        for mif in mol_el:
            if mif.tag.endswith('matchIf') and mif.get('feature') == 'substructureMatch':
                is_match = True
                break
        if not is_match:
            continue

        mol = _parse_molecule(mol_el)
        if not mol:
            print('  No atoms parsed!')
            continue

        print(f'\nMatch molecule id={mol_el.get("id")}:')
        for aid, ainfo in mol['atoms'].items():
            s = _atom_to_smarts(ainfo['element'], ainfo['matchifs'])
            print(f'  atom {aid}({ainfo["element"]}): matchifs={ainfo["matchifs"]} -> SMARTS={s}')
        print('Bonds:', mol['bonds'])
        print('DFS start:', list(mol['atoms'].keys())[0])
        smarts = _graph_to_smarts(mol['atoms'], mol['bonds'])
        print('Generated SMARTS:', smarts)
    break

print('\nDone')
