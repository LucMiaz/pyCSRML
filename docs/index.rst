pyCSRML documentation
======================

**pyCSRML** is a pure-Python implementation of the Chemical Subgraph Representation
Markup Language (CSRML). It parses CSRML XML files, converts subgraph patterns to
SMARTS, and computes binary chemotype fingerprints for molecules using RDKit.

Two fingerprint sets are bundled:

* **ToxPrint v2.0** — 729 general toxicologically-relevant chemotypes
* **TxP_PFAS v1.0** — 129 per- and polyfluoroalkyl substance (PFAS) chemotypes

Concordance with the reference ChemoTyper tool: **99.4 %** overall accuracy
across 14 710 compounds (Richard *et al.*, 2023).

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
