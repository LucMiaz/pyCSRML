Changelog
=========

0.1.0 (2026-03-20)
-------------------

Initial public release.

**Features:**

* Parse CSRML XML (ToxPrint v2 and TxP_PFAS v1) to SMARTS
* Bundled fingerprint definitions: ToxPrint v2.0 (729 bits), TxP_PFAS v1.0.4 (129 bits)
* ``PFASFingerprinter`` and ``ToxPrintFingerprinter`` classes using RDKit
* ``Embedding`` and ``EmbeddingSet`` containers for multi-compound analysis
* Support for all key CSRML features:

  - Pseudo-elements G (metals/metalloids), Z / Q (heteroatoms), X (halogens)
  - ``atomList negate="true"``
  - ``attachedHydrogenCount`` range expressions (minInclusive / minExclusive)
  - ``combineAtomFeatures`` OR-of-AND trees (QRY atoms)
  - ``matchingQueryAtom`` exception folding → ``[!$(…)]`` recursive SMARTS
  - ``ringCountAtom``
  - ``atomHeteroAttachedCount`` → recursive SMARTS ``$([elem](~[!#6;!#1])…)``

* JSON and YAML input support in :class:`~pyCSRML.Fingerprinter`
* 99.44 % overall concordance with ChemoTyper reference (Richard et al., 2023);
  all 129 bits ≥ 90 % accuracy
