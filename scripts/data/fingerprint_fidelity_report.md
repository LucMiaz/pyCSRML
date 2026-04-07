# Fingerprint Fidelity: pyCSRML vs ChemoTyper

---

## ToxCast — ToxPrint v2

**Dataset**: ToxCast (9014 compounds)  
**Fingerprint**: ToxPrint v2 (729 bits)  

### Per-compound Tanimoto similarity

| Statistic | Value |
|-----------|-------|
| Mean Tanimoto | 0.8602 |
| Median Tanimoto | 0.8846 |
| Std | 0.1240 |
| p1 / p5 / p25 | 0.4124 / 0.6190 / 0.8000 |
| p75 / p95 / p99 | 0.9412 / 1.0000 / 1.0000 |
| Maximum | 1.0000 |
| Identical (T=1.0) | 1,660 / 9,014 (18.4%) |

### Overall bit-level accuracy

Overall concordance: **99.6993%** (6,551,449 / 6,571,206 bit positions agree)

### Per-bit accuracy summary

| Category | Count | Fraction |
|----------|-------|---------|
| Perfect (acc = 1.00) | 514 | 70.5% |
| acc ≥ 0.99 | 697 | 95.6% |
| acc ≥ 0.95 | 721 | 98.9% |
| acc < 0.99 | 32 | 4.4% |

See `toxcast_per_bit_accuracy.csv` for the full per-bit table.

### Bits with accuracy < 99%

| Bit name | Accuracy | CT prevalence | PY prevalence |
|----------|----------|--------------|--------------|
| `ring:hetero_[6]_Z_generic` | 0.5462 | 0.2425 | 0.6964 |
| `chain:alkaneBranch_isopropyl_C3` | 0.7414 | 0.1162 | 0.3748 |
| `chain:alkaneCyclic_ethyl_C2_(connect_noZ)` | 0.7592 | 0.1766 | 0.4175 |
| `chain:alkeneCyclic_ethene_generic` | 0.8695 | 0.1030 | 0.1735 |
| `chain:alkeneCyclic_ethene_C_(connect_noZ)` | 0.9049 | 0.0576 | 0.1389 |
| `bond:quatN_b-carbonyl` | 0.9334 | 0.0007 | 0.0672 |
| `ring:hetero_[6_6]_Z_generic` | 0.9343 | 0.0714 | 0.1371 |
| `bond:CC(=O)C_ketone_alkene_cyclic_2-en-1-one_generic` | 0.9396 | 0.0476 | 0.0363 |
| `bond:CC(=O)C_ketone_alkane_cyclic` | 0.9586 | 0.0489 | 0.0721 |
| `ring:hetero_[5_6]_Z_generic` | 0.9612 | 0.0839 | 0.1227 |
| `bond:CC(=O)C_ketone_alkene_cyclic_2-en-1-one` | 0.9665 | 0.0230 | 0.0258 |
| `chain:alkeneCyclic_diene_cyclohexene` | 0.9700 | 0.0443 | 0.0150 |
| `chain:alkeneLinear_mono-ene_ethylene` | 0.9718 | 0.0413 | 0.0694 |
| `bond:C(=O)N_carboxamide_generic` | 0.9719 | 0.2060 | 0.1797 |
| `bond:NC=O_aminocarbonyl_generic` | 0.9727 | 0.2068 | 0.1797 |
| `bond:CC(=O)C_ketone_aliphatic_acyclic` | 0.9727 | 0.0470 | 0.0743 |
| `bond:C=O_carbonyl_ab-unsaturated_generic` | 0.9758 | 0.0605 | 0.0363 |
| `chain:alkeneLinear_mono-ene_allyl` | 0.9774 | 0.0371 | 0.0597 |
| `bond:CN_amine_ter-N_generic` | 0.9785 | 0.1305 | 0.1487 |
| `bond:CN_amine_ter-N_aliphatic` | 0.9787 | 0.1302 | 0.1487 |
| `chain:oxy-alkaneLinear_ethyleneOxide_EO10` | 0.9825 | 0.0001 | 0.0176 |
| `chain:oxy-alkaneLinear_ethyleneOxide_EO1` | 0.9827 | 0.0159 | 0.0014 |
| `chain:aromaticAlkane_Ph-C1_cyclic` | 0.9827 | 0.1076 | 0.0903 |
| `ring:hetero_[6]_N_pyrimidine` | 0.9849 | 0.0193 | 0.0344 |
| `bond:CC(=O)C_ketone_alkene_cyclic_(C6)` | 0.9868 | 0.0145 | 0.0013 |
| `bond:CC(=O)C_ketone_generic` | 0.9877 | 0.1047 | 0.0924 |
| `chain:oxy-alkaneLinear_ethylenOxide_EO1(O)` | 0.9878 | 0.0116 | 0.0006 |
| `bond:CN_amine_aliphatic_generic` | 0.9881 | 0.2080 | 0.2199 |
| `bond:NC=O_urea_generic` | 0.9884 | 0.0294 | 0.0178 |
| `bond:CC(=O)C_ketone_alkene_generic` | 0.9884 | 0.0416 | 0.0300 |

---

## CLinventory — ToxPrint v2

**Dataset**: 181,745 compounds  
**Fingerprint**: 729 bits  

### Per-compound Tanimoto similarity

| Statistic | Value |
|-----------|-------|
| Mean Tanimoto | 0.8849 |
| Median Tanimoto | 0.9091 |
| Std | 0.1145 |
| p1 / p5 / p25 | 0.4348 / 0.6667 / 0.8462 |
| p75 / p95 / p99 | 0.9565 / 1.0000 / 1.0000 |
| Maximum | 1.0000 |
| Identical (T=1.0) | 44,284 / 181,745 (24.4%) |

### Overall bit-level accuracy

Overall concordance: **99.7744%** (132,193,244 / 132,492,105 bit positions agree)

### Per-bit accuracy summary

| Category | Count | Fraction |
|----------|-------|---------|
| acc ≥ 0.99 | 706 | 96.8% |

---

## CLinventory — TxP_PFAS v1

**Dataset**: 181,745 compounds  
**Fingerprint**: 129 bits  

### Per-compound Tanimoto similarity

| Statistic | Value |
|-----------|-------|
| Mean Tanimoto | 0.9991 |
| Median Tanimoto | 1.0000 |
| Std | 0.0232 |
| p1 / p5 / p25 | 1.0000 / 1.0000 / 1.0000 |
| p75 / p95 / p99 | 1.0000 / 1.0000 / 1.0000 |
| Maximum | 1.0000 |
| Identical (T=1.0) | 181,425 / 181,745 (99.8%) |

### Overall bit-level accuracy

Overall concordance: **99.9983%** (23,444,695 / 23,445,105 bit positions agree)

### Per-bit accuracy summary

| Category | Count | Fraction |
|----------|-------|---------|
| acc ≥ 0.99 | 129 | 100.0% |