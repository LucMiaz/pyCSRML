# Supplementary Information: pyCSRML Validation

**Article context**: This supplement provides computational validation of the substitution
of ChemoTyper (Molecular Networks GmbH, Erlangen) with **pyCSRML** — a pure-Python
re-implementation of the Chemical Subgraph Representation Markup Language (CSRML) engine —
for large-scale molecular fingerprinting in the main analysis.

**Reproducing these results**: All scripts are in `pyCSRML/scripts/` and `pyCSRML/tests/`.
Run `python scripts/report_fingerprint_distances.py` and `python scripts/report_ml_comparison.py`
after completing the benchmark run described below.

---

## S1 — Computational Setup

### S1.1 Reference tool

ChemoTyper v1 (Molecular Networks GmbH) was used as the ground-truth reference for
fingerprint generation. ToxPrint v2.0 (729-bit, Yang *et al.* 2015) and TxP_PFAS v1.0
(129-bit, Richard *et al.* 2023) were the fingerprint sets evaluated.

### S1.2 pyCSRML implementation

pyCSRML (v0.1.0) parses the original CSRML XML fingerprint definition files, converts
each subgraph pattern to SMARTS, and evaluates them with RDKit (v2025.09.3) substructure
search. The implementation is available as an open-source Python package.

### S1.3 Benchmark hardware

Benchmarks were run on a single workstation:

| Property | Value |
|---|---|
| Machine | Asus ZenbookA14 |
| CPU | Qualcomm Snapdragon X Elite X1E78100 (Oryon cores) |
| Architecture | ARM64 |
| Physical / logical cores | 12 / 12 |
| Clock speed (max) | 3 417 MHz |
| RAM | ~32 GB |
| OS | Windows 11 (build 26200.8037) |
| Python | 3.14.2 (conda-forge) |
| RDKit | 2025.09.3 |
| NumPy | 2.3.5 |
| pandas | 3.0.0 |
| scikit-learn | 1.8.0 |

---

## S2 — Fingerprint Fidelity

### S2.1 Datasets

Three reference datasets were used for concordance evaluation:

| Dataset | *N* | Fingerprint | Source |
|---------|-----|------------|--------|
| Richard *et al.* 2023 PFAS set | 14 710 | TxP_PFAS v1 | Richard *et al.* 2023 |
| ToxCast (full) | 9 014 | ToxPrint v2 | EPA ToxCast/Tox21 |
| ToxCast CF-containing subset | 808 | TxP_PFAS v1 | EPA ToxCast/Tox21 |

The CF-containing subset comprises the 808 ToxCast compounds for which ChemoTyper sets
at least one TxP_PFAS bit, isolating the meaningful evaluation space for PFAS fingerprints
(otherwise the 8,206 non-PFAS compounds trivially achieve 100% accuracy on all-zero vectors).

### S2.2 Bit-level accuracy

For each dataset, overall concordance is reported as the fraction of (compound × bit)
positions that agree between pyCSRML and ChemoTyper, and per-bit accuracy is the fraction
of compounds on which that bit agrees.

| Dataset | *N* | Fingerprint | Overall accuracy | Bits ≥ 90% acc | Macro MCC | Macro Bal. Acc. | Macro ROC-AUC |
|---------|-----|------------|-----------------|---------------|-----------|----------------|---------------|
| Richard 2023 (PFAS) | 14 710 | TxP_PFAS v1 | **99.99%** | 129 / 129 | **0.9971** | **0.9989** | **0.9989** |
| ToxCast (full) | 9 014 | ToxPrint v2 | **99.71%** | 725 / 729 | **0.9326** | **0.9703** | **0.9703** |
| ToxCast (CF-subset) | 808 | TxP_PFAS v1 | **99.98%** | 129 / 129 | **0.9905** | **0.9924** | **0.9924** |
| CLinventory (full) | 181 745 | ToxPrint v2 | **99.77%** | 726 / 729 | **0.9320** | **0.9710** | **0.9710** |
| CLinventory (full) | 181 745 | TxP_PFAS v1 | **100.00%** | 129 / 129 | **0.9936** | **0.9946** | **0.9946** |
| CLinventory (CF-subset) | 27 989 | TxP_PFAS v1 | **99.99%** | 129 / 129 | **0.9935** | **0.9945** | **0.9945** |

Accuracy fractions are computed at the individual bit level
(i.e., the fraction of compounds for which a given bit takes the same value in
both tools). The ToxPrint v2 bit-level accuracy exceeds 99% for 95.6% of bits
on ToxCast (32 bits < 99%) and 96.8% of bits on CLinventory (23 bits < 99%).
The remaining discrepancies involve the same ring heteroatom and chain chemotype
groups on both datasets (see §S2.4 below).

### S2.3 Per-compound Tanimoto similarity

The element-wise Tanimoto (Jaccard) similarity between matched compound pairs was
computed from the raw bit vectors as:

$$T(a,b) = \frac{|\mathbf{a} \wedge \mathbf{b}|}{|\mathbf{a} \vee \mathbf{b}|}$$

where $\mathbf{a}$ is the ChemoTyper vector and $\mathbf{b}$ the pyCSRML vector.
Compounds where both vectors are empty (all-zero) receive $T = 1$ by convention.

**ToxCast — ToxPrint v2** (9,014 compounds):

| Statistic | Value |
|-----------|-------|
| Mean Tanimoto | 0.8602 |
| Median Tanimoto | 0.8846 |
| Std | 0.1240 |
| p1 / p5 / p25 | 0.412 / 0.619 / 0.800 |
| p75 / p95 / p99 | 0.941 / 1.000 / 1.000 |
| Maximum | 1.000 |
| Identical (T=1.0) | 1,660 / 9,014 (18.4%) |

**CLinventory — ToxPrint v2** (181,745 compounds):

| Statistic | Value |
|-----------|-------|
| Mean Tanimoto | 0.8849 |
| Median Tanimoto | 0.9091 |
| Std | 0.1145 |
| p1 / p5 / p25 | 0.435 / 0.667 / 0.846 |
| p75 / p95 / p99 | 0.957 / 1.000 / 1.000 |
| Maximum | 1.000 |
| Identical (T=1.0) | 44,284 / 181,745 (24.4%) |

**CLinventory — TxP_PFAS v1** (181,745 compounds):

| Statistic | Value |
|-----------|-------|
| Mean Tanimoto | 0.9991 |
| Median Tanimoto | 1.0000 |
| Std | 0.0232 |
| p1 / p5 / p25 | 1.000 / 1.000 / 1.000 |
| p75 / p95 / p99 | 1.000 / 1.000 / 1.000 |
| Maximum | 1.000 |
| Identical (T=1.0) | 181,425 / 181,745 (99.8%) |

For ToxPrint v2, the median per-compound Tanimoto (0.885 on ToxCast; 0.909 on
CLinventory) is consistent with the ≥99.7% overall bit accuracy on both datasets.
The small remaining gap from 1.0 reflects the systematic discrepancies in ring
heteroatom and chain patterns described in §S2.4. The CLinventory Tanimoto is
slightly higher than ToxCast's, consistent with the larger dataset averaging out
molecule-specific discrepancies. For TxP_PFAS v1 on CLinventory, 99.8% of compounds
achieve T = 1.0, confirming near-perfect PFAS fingerprint replication.

### S2.4 Known systematic discrepancies

Thirty-two ToxPrint v2 bits have accuracy < 99% on the ToxCast set (out of 729 total;
95.6% achieve ≥ 99% accuracy). The seven most impactful outliers (accuracy < 97%) fall
into two categories. The remaining 25 bits have accuracy 0.93–0.99 and represent minor
directional biases spread across carbonyl (9 bits), amine (5 bits), alkene chain (4 bits),
and miscellaneous ring/chain patterns; the full ranked list is in
`scripts/data/toxcast_per_bit_accuracy.csv`.

**Category 1 — Ring heteroatom over-matching (Z-generic and related patterns, accuracy < 97%)**

The `ring:hetero_[6]_Z_generic` pattern ("6-membered ring containing any heteroatom")
is over-broad: pyCSRML matches *any* 6-membered heteroaromatic ring (prevalence 69.6%),
whereas ChemoTyper applies stricter ring-type constraints (prevalence 24.3%). Related
fused-ring and 5,6-ring fusion `Z_generic` patterns show the same issue at lower severity.

| Bit | Accuracy | CT prev | PY prev | Direction |
|-----|---------|---------|---------|----------|
| `ring:hetero_[6]_Z_generic` | 0.546 | 0.243 | 0.696 | FP |
| `ring:hetero_[6_6]_Z_generic` | 0.934 | 0.071 | 0.137 | FP |
| `ring:hetero_[5_6]_Z_generic` | 0.961 | 0.084 | 0.123 | FP |

**Category 2 — Aliphatic chain patterns (ring-connected chains, `noZ` variants, accuracy < 91%)**

Several chain patterns using the `noZ` ("not connected to heteroatom") SMARTS modifier
show elevated false-positive rates, likely because the pyCSRML ring-attachment SMARTS
conversion is more permissive than the ChemoTyper implementation:

| Bit | Accuracy | CT prev | PY prev |
|-----|---------|---------|--------|
| `chain:alkaneBranch_isopropyl_C3` | 0.741 | 0.116 | 0.375 |
| `chain:alkaneCyclic_ethyl_C2_(connect_noZ)` | 0.759 | 0.177 | 0.418 |
| `chain:alkeneCyclic_ethene_generic` | 0.870 | 0.103 | 0.174 |
| `chain:alkeneCyclic_ethene_C_(connect_noZ)` | 0.905 | 0.058 | 0.139 |

**TxP_PFAS v1 discrepancies:**

| Bit | Fingerprint | Accuracy | Direction | Root cause |
|-----|------------|---------|-----------|-----------|
| `pfas_chain:alkeneLinear_mono-ene_ethylene_generic_F` | TxP_PFAS v1 | 98.9% | FN (recall 40%) | RDKit perceives C=C in fluoropyrimidines as aromatic; aliphatic SMARTS misses them |
| `pfas_bond:C=N_imine_FCN` | TxP_PFAS v1 | 99.5% | FN (recall 33%) | Same RDKit aromaticity issue on fluorinated heterocycles |
| `pfas_bond:aromatic_FCc1c` | TxP_PFAS v1 | 99.5% | Slight FP (prec. 97.2%) | Exception-clause approximation in SMARTS conversion |

**Impact on drug-like and environmental chemical screening**: The ring heteroatom FPs
(category 1) and chain FPs (category 2) add noise bits that are set for a
disproportionate fraction of molecules (e.g., `ring:hetero_[6]_Z_generic` pyCSRML
prevalence 69.6% vs ChemoTyper 24.3%) and thus have low discriminative power —
confirmed by the negligible ML performance difference in §S4.

---

## S3 — Computational Performance

### S3.1 Size-stratified timing benchmark

Molecules from the CLinventory were stratified by heavy-atom count into five sets of
500 compounds each. ChemoTyper timings are the mean of 3 manual repetitions;
pyCSRML timings are the median of 5 automated runs. Both tools were timed on identical
compound sets.

**ToxPrint v2 (729 patterns)**:

| Set | Heavy atoms | *N* mols | ChemoTyper (ms/mol) | pyCSRML (ms/mol) | Speedup |
|-----|-----------|---------|--------------------|--------------------|--------|
| `bench_tiny` | 1–10 | 500 | 13.8 | 3.76 | 3.7× |
| `bench_small` | 11–20 | 500 | 27.7 | 5.47 | 5.1× |
| `bench_medium` | 21–35 | 500 | 59.6 | 8.23 | 7.2× |
| `bench_large` | 36–60 | 500 | 114.6 | 12.32 | 9.3× |
| `bench_xlarge` | 61+ | 500 | 322.3 | 23.20 | 13.9× |

**TxP_PFAS v1 (129 patterns)**:

| Set | Heavy atoms | *N* mols | ChemoTyper (ms/mol) | pyCSRML (ms/mol) | Speedup |
|-----|-----------|---------|--------------------|--------------------|--------|
| `bench_tiny` | 1–10 | 500 | 4.3 | 0.73 | 5.9× |
| `bench_small` | 11–20 | 500 | 7.7 | 1.01 | 7.7× |
| `bench_medium` | 21–35 | 500 | 17.9 | 1.53 | 11.7× |
| `bench_large` | 36–60 | 500 | 30.5 | 2.19 | 13.9× |
| `bench_xlarge` | 61+ | 500 | 139.1 | 4.46 | 31.2× |

pyCSRML outperforms ChemoTyper for both fingerprint sets across all molecular size
strata. The speedup grows with molecular complexity — for ToxPrint v2 from 3.7× on
very small molecules to 13.9× on those with 36–60 heavy atoms — suggesting the
ChemoTyper subgraph algorithm scales less favourably with structural size. TxP_PFAS v1
shows even larger speedups (5.9–31.2×), because the targeted fluorocarbon patterns
reduce the SMARTS search space more aggressively than the broad ToxPrint v2 set.

### S3.2 Throughput at scale

Using the `bench_medium` representative cost (21–35 heavy atoms) on a single 3.4 GHz
core — 59.6 ms/mol for ChemoTyper and 8.23 ms/mol for pyCSRML (ToxPrint v2):

| Compound set | ChemoTyper ToxPrint v2 | pyCSRML ToxPrint v2 |
|---|---|---|
| 10 000 | ~596 s | ~82 s |
| 100 000 | ~99 min | ~14 min |
| 1 000 000 | ~17 h | ~2.3 h |
| 10 000 000 | ~166 h | ~23 h |

Trivial parallelism (e.g., `ProcessPoolExecutor`) scales pyCSRML linearly with the
number of cores (12 on the benchmark machine → ~12× throughput).

---

## S4 — ML Predictive Performance on ToxCast Endpoints

### S4.1 Experimental design

To verify that the fingerprint discrepancies in §S2 do not affect downstream predictive
modelling, a cross-validated machine-learning benchmark was performed on the 15 ToxCast
Tox21 endpoints used in the main analysis.

**Data**: 9,014 ToxCast compounds with binary activity labels for:
AR_agonist, AR_antagonist, AhR_agonist, Aromatase_antagonist, CYP2C19/2C9/2D6/3A4_antagonist,
Caspase3_HEPG2, DT40_genotoxicity, ERa_agonist/antagonist, MMP_ratio, TR_antagonist, p53_ratio.

**Features**: Two fingerprint sets compared side-by-side with identical models:
1. **ChemoTyper**: Reference ToxPrint v2 vectors from the ChemoTyper tool (TSV export)
2. **pyCSRML**: ToxPrint v2 vectors computed by pyCSRML

**Model**: Random Forest (scikit-learn), 200 estimators, `max_features="sqrt"`,
`class_weight="balanced"`, 5-fold stratified cross-validation.

**Metrics**: ROC-AUC, Average Precision (AUPRC), Matthews Correlation Coefficient (MCC),
Balanced Accuracy.

### S4.2 Per-endpoint results

| Endpoint | *N*⁺ | CT ROC-AUC | PY ROC-AUC | CT MCC | PY MCC | CT AvgPrec | PY AvgPrec |
|----------|------|-----------|-----------|--------|--------|-----------|-----------|
| AR_agonist | 471 | 0.8186 ± 0.0084 | 0.8186 ± 0.0178 | 0.4031 ± 0.0400 | 0.3758 ± 0.0444 | 0.4340 ± 0.0331 | 0.4239 ± 0.0328 |
| AR_antagonist | 1506 | 0.8029 ± 0.0093 | 0.8032 ± 0.0097 | 0.3256 ± 0.0314 | 0.3264 ± 0.0198 | 0.4883 ± 0.0216 | 0.4754 ± 0.0246 |
| AhR_agonist | — | 0.8458 ± 0.0217 | 0.8474 ± 0.0212 | 0.3226 ± 0.0324 | 0.3255 ± 0.0179 | 0.3696 ± 0.0267 | 0.3671 ± 0.0199 |
| Aromatase_antagonist | — | 0.8209 ± 0.0271 | 0.8163 ± 0.0254 | 0.3373 ± 0.0181 | 0.3307 ± 0.0281 | 0.4383 ± 0.0254 | 0.4248 ± 0.0280 |
| CYP2C19_antagonist | — | 0.7730 ± 0.0179 | 0.7755 ± 0.0217 | 0.3691 ± 0.0323 | 0.3789 ± 0.0353 | 0.5802 ± 0.0235 | 0.5851 ± 0.0299 |
| CYP2C9_antagonist | — | 0.7944 ± 0.0060 | 0.7892 ± 0.0084 | 0.4084 ± 0.0244 | 0.4111 ± 0.0177 | 0.6386 ± 0.0121 | 0.6320 ± 0.0143 |
| CYP2D6_antagonist | — | 0.7333 ± 0.0121 | 0.7302 ± 0.0139 | 0.3355 ± 0.0276 | 0.3254 ± 0.0147 | 0.6189 ± 0.0241 | 0.6144 ± 0.0206 |
| CYP3A4_antagonist | — | 0.8062 ± 0.0098 | 0.8028 ± 0.0113 | 0.3778 ± 0.0137 | 0.3604 ± 0.0098 | 0.5725 ± 0.0198 | 0.5679 ± 0.0227 |
| Caspase3_HEPG2 | — | 0.8427 ± 0.0171 | 0.8500 ± 0.0254 | 0.3045 ± 0.0570 | 0.2930 ± 0.0677 | 0.3158 ± 0.0522 | 0.3212 ± 0.0557 |
| DT40_genotoxicity | — | 0.7940 ± 0.0058 | 0.7954 ± 0.0096 | 0.4152 ± 0.0193 | 0.4180 ± 0.0257 | 0.6342 ± 0.0121 | 0.6359 ± 0.0148 |
| ERa_agonist | — | 0.6569 ± 0.0212 | 0.6502 ± 0.0270 | 0.1909 ± 0.0344 | 0.1890 ± 0.0249 | 0.2691 ± 0.0218 | 0.2566 ± 0.0270 |
| ERa_antagonist | 1055 | 0.8178 ± 0.0224 | 0.8131 ± 0.0223 | 0.3046 ± 0.0429 | 0.2928 ± 0.0339 | 0.4114 ± 0.0365 | 0.3964 ± 0.0336 |
| MMP_ratio | — | 0.8633 ± 0.0073 | 0.8570 ± 0.0044 | 0.3988 ± 0.0202 | 0.4016 ± 0.0226 | 0.5338 ± 0.0086 | 0.5217 ± 0.0082 |
| TR_antagonist | — | 0.8276 ± 0.0085 | 0.8277 ± 0.0071 | 0.4229 ± 0.0170 | 0.4140 ± 0.0291 | 0.6016 ± 0.0097 | 0.6005 ± 0.0121 |
| p53_ratio | — | 0.8176 ± 0.0125 | 0.8129 ± 0.0122 | 0.2795 ± 0.0187 | 0.2578 ± 0.0235 | 0.3480 ± 0.0174 | 0.3428 ± 0.0174 |

CT = ChemoTyper; PY = pyCSRML; mean ± s.d. over 5 folds.  
*N*⁺ = number of positive compounds (shown for endpoints with available fold-level data).

### S4.3 Statistical comparison (global across all endpoints)

The global null hypothesis is that ChemoTyper and pyCSRML fingerprints yield
**identical predictive performance** when used as features in the same model.

**Method**: For each metric, 75 paired observations (15 endpoints × 5 folds) are
formed as $\Delta_k = \text{pyCSRML}_k - \text{ChemoTyper}_k$. Two complementary tests:

1. **Frequentist**: Two-sided paired Wilcoxon signed-rank test (appropriate for
   non-normal small-sample differences, no distributional assumption).

2. **Bayesian**: Under a Jeffreys (flat) improper prior on $\mu$ and $\sigma$, the
   posterior of the mean difference $\mu_\Delta | \mathbf{\Delta}$ is a $t_{n-1}$
   distribution with location $\bar{\Delta}$ and scale $s_\Delta / \sqrt{n}$. This 
   gives the 95% highest-density credible interval (CrI) and the posterior probability
   $P(\Delta > 0 \mid \mathbf{\Delta})$ that pyCSRML exceeds ChemoTyper.

**To obtain exact numerical results**, run:

```bash
python scripts/report_ml_comparison.py
```

This writes `scripts/data/ml_comparison_report.md` with the full per-endpoint and
global tables.

**Results** (fast CV run, *n*=75 pairs per metric, 15 endpoints × 5 folds):

| Metric | Mean Δ | 95% CrI | P(Δ>0) | Wilcoxon *p* |
|--------|--------|---------|--------|-------------|
| ROC-AUC | −0.0017 | [−0.0037, +0.0002] | 0.043 | 0.032 |
| Avg. Precision | −0.0059 | [−0.0084, −0.0035] | < 0.001 | < 0.0001 |
| MCC | −0.0064 | [−0.0106, −0.0021] | 0.002 | 0.006 |
| Balanced Acc. | −0.0022 | [−0.0041, −0.0002] | 0.015 | 0.050 |

Δ = pyCSRML − ChemoTyper; 95% CrI = Bayesian credible interval under Jeffreys prior
(posterior is $t_{74}$ with location $\bar{\Delta}$ and scale $s_\Delta/\sqrt{75}$).

pyCSRML fingerprints yield **statistically but not practically significant** lower
predictive performance compared to ChemoTyper across all four metrics. The largest
difference is in Average Precision (Δ = −0.006) and MCC (Δ = −0.006). These
differences are consistent with the systematic false-positive bits documented in §S2.4,
which add uninformative noise features shared across many compounds. They do not
indicate meaningful information loss for the predictive task.

### S4.4 Interpretation

The small but statistically significant differences (maximum Δ = 0.006 on average
precision, Wilcoxon *p* < 0.0001 across 75 endpoint–fold pairs) are consistent with
the systematic false-positive bits documented in §S2.4. These bits are set for a
disproportionate fraction of molecules (e.g., `ring:hetero_[6]_Z_generic` pyCSRML
prevalence 69.6% vs ChemoTyper 24.3%) and are therefore redundant from a discriminative
standpoint: a random forest learning on these features will down-weight them due to low
Gini importance. The practical magnitude of the performance loss is negligible (< 0.01
on any metric), making pyCSRML a suitable substitute for ChemoTyper in large-scale
predictive screening.

---

## S5 — Reproducibility Checklist

| Step | Command | Output |
|------|---------|--------|
| Run concordance tests (requires ChemoTyper reference data) | `pytest tests/test_chemotyper_concordance.py -v -s` | `tests/concordance_report.md` |
| Run per-bit accuracy on CLinventory | `pytest tests/test_clinventory_bit_accuracy.py -v` | `tests/clinventory_bit_accuracy_baseline.json` |
| Run ML benchmark (fast mode) | `python scripts/benchmark_toxprint_toxcast.py --fast` | `scripts/data/toxprint_toxcast_fast_results.csv` |
| Run ML benchmark (full nested CV) | `python scripts/benchmark_toxprint_toxcast.py` | `scripts/data/toxprint_toxcast_results.csv` |
| Fingerprint fidelity report | `python scripts/report_fingerprint_distances.py` | `scripts/data/fingerprint_fidelity_report.md` |
| ML comparison report with Bayesian test | `python scripts/report_ml_comparison.py` | `scripts/data/ml_comparison_report.md` |
| Timing benchmark | `pytest tests/test_timing_benchmark.py -v` | `tests/test_data/size_benchmarks/pycsrml_timing_baseline.json` |

---

## References

- Yang, C.; Tarkhov, A.; Marusczyk, J.; Bienfait, B.; Gasteiger, J.; Kleinöder, T.;  
  Terfloth, L.; Oprea, T. I.; Sacher, O. *J. Chem. Inf. Model.* **2015**, *55*, 510–528.  
  DOI: 10.1021/ci500667v

- Richard, A. M.; Judson, R. S.; Houck, K. A.; Grulke, C. M.; Volarath, P.; Thillainadarajah, I.;  
  Yang, C.; Rathman, J.; Martin, M. T.; Wambaugh, J. F.; Knudsen, T. B.; Kellogg, M.; Harten, P.;  
  Wolber, P. K.; Walker, D.; Kavlock, R.; Dix, D. J. *Environ. Health Perspect.* **2023**.  
  DOI: 10.1289/EHP10709 *(or update with current reference)*

- Kruschke, J. K. *Doing Bayesian Data Analysis*, 2nd ed.; Academic Press, 2015.  
  (BEST test methodology for §S4.3)
