import json
from pathlib import Path

data = json.loads(Path(r'C:\Users\luc\git\pyCSRML\pyCSRML\data\TxP_PFAS_v1.0.4.json').read_text())
bits = {b['label']: b for b in data['bits']}

targets = [
    'pfas_chain:FT_n1_hetero', 'pfas_chain:FT_n2_hetero', 'pfas_chain:FT_n3_hetero',
    'pfas_chain:FT_n1_O', 'pfas_chain:FT_n2_O', 'pfas_chain:FT_n1_S', 'pfas_chain:FT_n2_S',
    'pfas_chain:FT_n1_N', 'pfas_chain:FT_n2_N',
    'pfas_chain:FT_n1_X', 'pfas_chain:FT_n2_X',
    'pfas_chain:FT_n1_OP', 'pfas_chain:FT_n2_OP',
    'pfas_chain:FT_n1_OS', 'pfas_chain:FT_n1_C=O', 'pfas_chain:FT_n2_C=O',
    'pfas_chain:FT_n2_Si',
    'pfas_chain:polyF_cap_CHF2CF', 'pfas_chain:polyF_cap_CH2FCF',
    'pfas_chain:polyF_nocap_CFCHFCF', 'pfas_chain:polyF_nocap_CFCH2CF',
    'pfas_chain:perF-linear_cap_C6_excl_mod', 'pfas_chain:perF-linear_cap_C7_excl_mod',
    'pfas_chain:perF-linear_cap_C8_excl_mod',
    'pfas_bond:COH_alcohol_generic_OCCF',
    'pfas_bond:C(=O)N_carboxamide_(NHR)_C(=O)CF',
    'pfas_bond:CN_amine_sec-NH_alkyl_CF', 'pfas_bond:CN_amine_pri-NH2_alkyl_CF',
    'pfas_bond:C=O_aldehyde_alkyl_CF',
]

for name in targets:
    b = bits.get(name, {})
    if not b:
        print(f'--- {name}: NOT FOUND ---')
        continue
    ps = b.get('smarts', 'N/A')
    excs = b.get('exception_smarts', [])
    print(f'--- {name} ---')
    print(f'  primary:    {ps[:120]}')
    if len(ps) > 120:
        print(f'  (cont...)   {ps[120:240]}')
    if excs:
        for e in excs:
            print(f'  exception:  {e[:100]}')
    print()
