"""Dump raw XML for problem subgraphs."""
import xml.etree.ElementTree as ET, sys
sys.path.insert(0, r'C:\Users\luc\git\pyCSRML')

XML = r'C:\Users\luc\git\pyCSRML\pyCSRML\data\TxP_PFAS_v1.0.4.xml'
ns = 'http://www.molecular-networks.com/schema/csrml'
ids = [
    'txp-pfas-001',  # element_metal_metalloid_CF
    'txp-pfas-038',  # OZ_oxide_hydroxy_CF
    'txp-pfas-100',  # perF-linear_cap_C1_excl
    'txp-pfas-101',  # perF-linear_cap_C2_excl
    'txp-pfas-108',  # perF-linear_cap_C6_excl_mod
    'txp-pfas-109',  # perF-linear_cap_C7_excl_mod
    'txp-pfas-113',  # perF-linear_nocap_C1_excl
]

tree = ET.parse(XML)
root = tree.getroot()
for sg in root.iter(f'{{{ns}}}subgraph'):
    if sg.get('id') not in ids:
        continue
    print(f"\n{'='*70}")
    ET.indent(sg)
    print(ET.tostring(sg, encoding='unicode')[:4000])
    print()
