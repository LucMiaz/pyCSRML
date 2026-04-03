pyCSRML documentation
======================

**pyCSRML** is an attempt to provide a pure-Python implementation of the Chemical Subgraph Representation Markup Language (CSRML). It parses CSRML XML files, converts subgraph patterns to SMARTS, and computes binary chemotype fingerprints for molecules using RDKit.

Two fingerprint sets are bundled:

* **ToxPrint v2.0** — 729 general toxicologically-relevant chemotypes
* **TxP_PFAS v1.0** — 129 per- and polyfluoroalkyl substance (PFAS) chemotypes

Concordance with the reference ChemoTyper tool:

* **99.4 %** overall accuracy across 14 710 PFAS compounds from Richard
  *et al.* (2023), using the TxP_PFAS v1.0 fingerprint (bits ≥ 90 %:
  see ``tests/concordance_report.md`` after re-running the test suite).
* **98.17 %** overall accuracy on the full ToxCast library (9 014 compounds,
  711 / 729 bits ≥ 90 %) for ToxPrint v2.0.
* **99.98 %** overall accuracy on the ToxCast CF-containing subset (808 compounds,
  all 129 bits ≥ 90 %) for TxP_PFAS v1.0.

Run ``pytest tests/test_chemotyper_concordance.py -v -s`` to reproduce;
results are written to ``tests/concordance_report.md``.  Known residual
discrepancies (G/Q pseudo-element approximations, RDKit vs ChemoTyper
aromaticity for fluorinated heterocycles) are documented in ``README.md``
under the *Known discrepancies* heading.

.. note::

   See ``README.md`` for the performance table and a "Known discrepancies"
   section explaining systematic differences vs the reference tool.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   quickstart
   csrml_format
   json_yaml_format

.. toctree::
   :maxdepth: 2
   :caption: API reference

   api/fingerprinter
   api/embedding
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
