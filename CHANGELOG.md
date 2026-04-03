# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - 2026-04-03

Updated documentation and preparation for publication.

## [0.1.0] — 2026-03-20

Initial module.

### Added
- Parse CSRML XML (ToxPrint v2 and TxP_PFAS v1) to SMARTS
- Bundled fingerprint definitions: ToxPrint v2.0 (729 bits), TxP_PFAS v1.0.4 (129 bits)
- `PFASFingerprinter` and `ToxPrintFingerprinter` classes using RDKit
- `Embedding` and `EmbeddingSet` containers for multi-compound analysis
- Full support for all key CSRML atom features (pseudo-elements G/Z/Q/X, negated atomList,
  H-count ranges, combineAtomFeatures OR-of-AND trees, matchingQueryAtom folding,
  atomHeteroAttachedCount via recursive SMARTS)
- JSON and YAML input support in `Fingerprinter` (in addition to XML)
- 99.4 % concordance with ChemoTyper reference (Richard *et al.*, 2023); all 129 bits now ≥ 90 % accuracy
