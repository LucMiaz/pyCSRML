# Changelog

All notable changes to this project will be documented in this file.

## [0.3.0] - 2026-04-07

Improved accuracy on worst bits for both ToxPrint and TxP_PFAS. Renewed performance statistics. Added comparison of endpoint prediction.

## [0.2.0] - 2026-04-03

Updated documentation and preparation for publication.

## [0.1.0] — 2026-03-20

Initial module.

### Added
- Parse CSRML XML (ToxPrint v2 and TxP_PFAS v1) to SMARTS
- Bundled fingerprint definitions: ToxPrint v2.0 (729 bits) and TxP\_PFAS v1.0.4 (129 bits), accessible via `TOXPRINT_PATH` and `TXPPFAS_PATH`
- `Fingerprinter` class using RDKit; accepts XML, JSON, or YAML definitions
- Full support for all key CSRML atom features (pseudo-elements G/Z/Q/X, negated atomList,
  H-count ranges, combineAtomFeatures OR-of-AND trees, matchingQueryAtom folding,
  atomHeteroAttachedCount via recursive SMARTS)
- JSON and YAML input support in `Fingerprinter` (in addition to XML)
- 99.4 % concordance with ChemoTyper reference (Richard *et al.*, 2023); all 129 bits now ≥ 90 % accuracy
