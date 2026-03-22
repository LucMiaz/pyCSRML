"""Check matchingQueryAtom exception values in XML for CN_amine."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')
import json
import xml.etree.ElementTree as ET

with open(r'C:\Users\luc\git\pyCSRML\pyCSRML\data\TxP_PFAS_v1.0.4.json') as f:
    j = json.load(f)

TARGET_BITS = {
    'pfas_bond:aromatic_FCc1c',
    'pfas_bond:CN_amine_ter-N_generic_CF',
    'pfas_bond:C=N_imine_FCN',
    'pfas_chain:alkeneLinear_mono-ene_ethylene_generic_F',
}
print("=== Current SMARTS ===")
for b in j['bits']:
    if b['label'] in TARGET_BITS:
        print(f"{b['label']}")
        print(f"  smarts:      {b['smarts']}")
        print(f"  exceptions:  {b['exception_smarts']}")
        print()

# Inspect CN_amine exception XML
NS1 = 'http://www.openrisknet.org/schema/csRML/v1'
NS = {'csrml': NS1}
tree = ET.parse(r'C:\Users\luc\git\pyCSRML\pyCSRML\data\TxP_PFAS_v1.0.4.xml')
root = tree.getroot()

print("=== CN_amine exception molecules in XML ===")
for sg in root.iter():
    if sg.tag.endswith('subgraph'):
        # find label text
        label = ''
        for child in sg:
            if child.tag.endswith('label'):
                label = child.text or ''
                break
        if 'CN_amine_ter-N_generic' in label:
            for mol_el in sg:
                if not mol_el.tag.endswith('molecule'):
                    continue
                # check matchIf features
                for atom_el in mol_el.iter():
                    if not atom_el.tag.endswith('atom'):
                        continue
                    for mif in atom_el:
                        if not mif.tag.endswith('matchIf'):
                            continue
                        feat = mif.get('feature','')
                        if feat == 'matchingQueryAtom':
                            vals = [v.text for v in mif if v.tag.endswith('value')]
                            print(f"  mol id={mol_el.get('id')}, "
                                  f"atom id={atom_el.get('id')}, "
                                  f"matchingQueryAtom values={vals}")


