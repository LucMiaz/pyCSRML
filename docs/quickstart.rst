Quick start
===========

Computing fingerprints for a single molecule
--------------------------------------------

.. code-block:: python

   from pyCSRML import Fingerprinter, TOXPRINT_PATH
   from rdkit import Chem

   fp = Fingerprinter(TOXPRINT_PATH)   # loads bundled ToxPrint v2.0 (729 bits)

   mol = Chem.MolFromSmiles("c1ccccc1")   # benzene
   arr, names = fp.fingerprint(mol)

   print(f"Bits set: {arr.sum()} / {fp.n_bits}")
   on_bits = [names[i] for i in range(len(arr)) if arr[i]]
   print(on_bits[:5])

``arr`` is a boolean :class:`numpy.ndarray` of length ``fp.n_bits``.
``names`` is the matching list of chemotype labels (e.g. ``"ring:aro_6_C"``).


Using the bundled TxP\_PFAS definition
---------------------------------------

.. code-block:: python

   from pyCSRML import Fingerprinter, TXPPFAS_PATH
   from rdkit import Chem

   fp = Fingerprinter(TXPPFAS_PATH)   # loads bundled TxP_PFAS v1.0.4 (129 bits)

   mol = Chem.MolFromSmiles(
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"  # PFOA
   )
   arr, names = fp.fingerprint(mol)
   print(f"Bits set: {arr.sum()} / {fp.n_bits}")


Batch processing
----------------

:meth:`~pyCSRML.Fingerprinter.fingerprint_batch` processes a list of molecules
and returns a 2-D boolean NumPy matrix of shape ``(n_molecules, n_bits)``.

.. code-block:: python

   from pyCSRML import Fingerprinter, TXPPFAS_PATH
   from rdkit import Chem

   smiles = [
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",   # PFOA
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
       "CCO",     # ethanol (negative control)
   ]
   mols = [Chem.MolFromSmiles(s) for s in smiles]

   fp = Fingerprinter(TXPPFAS_PATH)
   matrix = fp.fingerprint_batch(mols)   # shape (3, 129), dtype bool
   print(matrix.shape, matrix.dtype)


Using a custom CSRML XML file
------------------------------

.. code-block:: python

   from pyCSRML import Fingerprinter

   fp = Fingerprinter("path/to/my_fingerprints.xml")
   mol = Chem.MolFromSmiles("CCO")
   arr, names = fp.fingerprint(mol)
