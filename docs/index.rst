pyCSRML documentation
======================

**pyCSRML** is an attempt to provide a pure-Python implementation of the Chemical Subgraph Representation Markup Language (CSRML). It parses CSRML XML files, converts subgraph patterns to SMARTS, and computes binary chemotype fingerprints for molecules using RDKit.

Two fingerprint sets are bundled and accessible via path constants:

* **ToxPrint v2.0** — 729 general toxicologically-relevant chemotypes (``TOXPRINT_PATH``)
* **TxP_PFAS v1.0** — 129 per- and polyfluoroalkyl substance (PFAS) chemotypes (``TXPPFAS_PATH``)

Concordance with the reference ChemoTyper tool:

* **99.99 %** overall accuracy across 14 710 PFAS compounds from Richard
  *et al.* (2023), using the TxP_PFAS v1.0 fingerprint.
* **98.17 %** overall accuracy on the full ToxCast library (9 014 compounds,
  711 / 729 bits ≥ 90 %) for ToxPrint v2.0.
* **99.98 %** overall accuracy on the ToxCast CF-containing subset (808 compounds,
  all 129 bits ≥ 90 %) for TxP_PFAS v1.0.
* **CLinventory** (pharmaceutical/industrial chemical set) — ToxPrint v2.0 and
  TxP_PFAS v1.0 results are generated automatically when running the test suite;
  see ``tests/concordance_report.md`` for the latest figures.

Run ``pytest tests/test_chemotyper_concordance.py -v -s`` to reproduce;
results (including per-bit MCC, balanced accuracy, and ROC-AUC) are written to
``tests/concordance_report.md``.  Known residual
discrepancies (G/Q pseudo-element approximations, RDKit vs ChemoTyper
aromaticity for fluorinated heterocycles) are documented in ``README.md``
under the *Known discrepancies* heading.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   quickstart
   benchmarks
   csrml_format
   json_yaml_format

.. toctree::
   :maxdepth: 2
   :caption: API reference

   api/fingerprinter
   api/csrml

.. toctree::
   :maxdepth: 1
   :caption: About

   changelog

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
