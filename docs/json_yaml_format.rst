Custom fingerprints: JSON and YAML format
==========================================

In addition to CSRML XML files, the :class:`~pyCSRML.Fingerprinter` accepts
fingerprint definitions written in **JSON** or **YAML**. This lets you:

* Distribute pre-built fingerprint specs without requiring the XML parser.
* Craft lightweight custom fingerprint sets by hand.
* Share fingerprint definitions that combine patterns from multiple sources.

The format is intentionally simple: a list of bit descriptors, each containing
an id, a human-readable label, a SMARTS pattern, and an optional list of
exception SMARTS.


File format
-----------

Both JSON and YAML use the same schema.  Two layouts are accepted:

**Wrapper dict (recommended)** — includes metadata:

.. code-block:: json

   {
     "id": "my-fingerprints",
     "title": "My custom chemotype fingerprints",
     "bits": [
       {
         "id": "fp-001",
         "label": "bond:CF_monofluoro",
         "smarts": "[#6]-[#9]",
         "exception_smarts": []
       },
       {
         "id": "fp-002",
         "label": "chain:CF2_gem-difluoro",
         "smarts": "[#6](-[#9])-[#9]",
         "exception_smarts": []
       }
     ]
   }

**Plain list** — bits only, no metadata:

.. code-block:: json

   [
     {
       "id": "fp-001",
       "label": "bond:CF_monofluoro",
       "smarts": "[#6]-[#9]",
       "exception_smarts": []
     }
   ]

Both layouts are also valid YAML. For example, the same spec in YAML:

.. code-block:: yaml

   id: my-fingerprints
   title: My custom chemotype fingerprints
   bits:
     - id: fp-001
       label: bond:CF_monofluoro
       smarts: "[#6]-[#9]"
       exception_smarts: []
     - id: fp-002
       label: chain:CF2_gem-difluoro
       smarts: "[#6](-[#9])-[#9]"
       exception_smarts: []


Field reference
---------------

.. list-table::
   :widths: 15 10 75
   :header-rows: 1

   * - Field
     - Required
     - Description
   * - ``id``
     - yes
     - Unique identifier for the bit (e.g. ``"fp-001"``).  Used internally
       but not exposed in the fingerprint output.
   * - ``label``
     - yes
     - Human-readable bit name returned by :attr:`~pyCSRML.Fingerprinter.bit_names`.
   * - ``smarts``
     - yes
     - RDKit-compatible SMARTS pattern for the primary substructure match.
       A bit is set if the molecule contains this substructure **and** none of
       the exception SMARTS match.
   * - ``exception_smarts``
     - no
     - List of SMARTS strings.  When any of these match the molecule, the bit
       is forced to 0 even if the main SMARTS matched.  Defaults to ``[]``.


Loading a custom file
---------------------

.. code-block:: python

   from pyCSRML import Fingerprinter
   from rdkit import Chem

   # JSON
   fp = Fingerprinter("my_fingerprints.json")

   # YAML (requires pyyaml: pip install pyyaml)
   fp = Fingerprinter("my_fingerprints.yaml")

   mol = Chem.MolFromSmiles("FC(F)(F)CCO")
   arr, names = fp.fingerprint(mol)
   print([n for n, b in zip(names, arr) if b])

The bundled pre-built JSON caches for ToxPrint and TxP_PFAS can also be loaded
directly, bypassing the XML parser entirely:

.. code-block:: python

   import importlib.resources
   from pyCSRML import Fingerprinter

   data_dir = importlib.resources.files("pyCSRML") / "data"
   fp = Fingerprinter(str(data_dir / "TxP_PFAS_v1.0.4.json"))

   from rdkit import Chem
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   arr, names = fp.fingerprint(mol)


Exporting the bundled specs to JSON
-------------------------------------

The :mod:`pyCSRML.convert_xml_to_json` script regenerates the bundled JSON
caches from the XML sources.  You can also export any XML file yourself:

.. code-block:: python

   from pyCSRML._csrml import parse_csrml_xml, ordered_bit_list
   import json

   data = parse_csrml_xml("my_fingerprints.xml")
   bits = []
   for bit_id in ordered_bit_list(data):
       sg = data["subgraph_index"][bit_id]
       bits.append({
           "id": sg["id"],
           "label": sg["label"],
           "smarts": sg["smarts"],
           "exception_smarts": sg["exception_smarts"],
       })

   spec = {"id": data["id"], "title": data["title"], "bits": bits}
   with open("my_fingerprints.json", "w") as f:
       json.dump(spec, f, indent=2)


YAML installation
-----------------

YAML support requires `PyYAML <https://pyyaml.org/>`_, which is not installed
by default:

.. code-block:: bash

   pip install pyyaml
   # or, with the yaml extra:
   pip install "pyCSRML[yaml]"
