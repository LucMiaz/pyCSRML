Installation
============

Requirements
------------

* Python 3.10 or later
* `RDKit <https://www.rdkit.org/>`_ ≥ 2022.3
* NumPy ≥ 1.23

Optional dependencies
---------------------

.. list-table::
   :widths: 25 30 45
   :header-rows: 1

   * - Extra
     - Packages installed
     - Enables
   * - ``analysis``
     - pandas ≥ 1.5, matplotlib ≥ 3.6
     - ``EmbeddingSet.plot()``, ``EmbeddingSet.to_dataframe()``
   * - ``umap``
     - umap-learn ≥ 0.5
     - ``EmbeddingSet.umap()``
   * - ``full``
     - all of the above
     - everything


From PyPI (pip)
---------------

.. code-block:: bash

   pip install pyCSRML

With optional extras:

.. code-block:: bash

   pip install "pyCSRML[analysis]"
   pip install "pyCSRML[full]"


Conda users
-----------

RDKit is best installed via conda-forge before installing pyCSRML with pip:

.. code-block:: bash

   conda install -c conda-forge rdkit
   pip install pyCSRML


Development install
-------------------

.. code-block:: bash

   git clone https://github.com/luc-miaz/pyCSRML
   cd pyCSRML
   pip install -e ".[dev]"
