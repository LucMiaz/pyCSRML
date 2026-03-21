Quick start
===========

Computing fingerprints for a single molecule
--------------------------------------------

.. code-block:: python

   from pyCSRML import PFASFingerprinter
   from rdkit import Chem

   fp = PFASFingerprinter()   # loads bundled TxP_PFAS v1.0 definitions

   mol = Chem.MolFromSmiles(
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"  # PFOA
   )
   arr, names = fp.fingerprint(mol)

   print(f"Bits set: {arr.sum()} / {fp.n_bits}")
   on_bits = [names[i] for i in range(len(arr)) if arr[i]]
   print(on_bits[:5])

``arr`` is a boolean :class:`numpy.ndarray` of length ``fp.n_bits``.
``names`` is the matching list of chemotype labels (e.g. ``"pfas_chain:perF-linear_C8_plus"``).


ToxPrint v2 fingerprints
-------------------------

.. code-block:: python

   from pyCSRML import ToxPrintFingerprinter
   from rdkit import Chem

   fp = ToxPrintFingerprinter()
   mol = Chem.MolFromSmiles("c1ccccc1")
   arr, names = fp.fingerprint(mol)
   print(f"Benzene: {arr.sum()} bits")


Using a custom CSRML XML file
------------------------------

.. code-block:: python

   from pyCSRML import Fingerprinter

   fp = Fingerprinter("path/to/my_fingerprints.xml")
   mol = Chem.MolFromSmiles("CCO")
   arr, names = fp.fingerprint(mol)


Batch processing with EmbeddingSet
------------------------------------

:class:`~pyCSRML.EmbeddingSet` collects fingerprints for a list of compounds and
provides analysis helpers.

.. code-block:: python

   from pyCSRML import PFASFingerprinter, from_fingerprinter

   fp = PFASFingerprinter()

   smiles = [
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",   # PFOA
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
       "FCCCF",   # short difluoro
       "CCO",     # ethanol (negative control)
   ]
   names = ["PFOA", "PFOS", "4F-propane", "EtOH"]

   eset = from_fingerprinter(fp, smiles_list=smiles, names=names)

   # Clustermap / heatmap (requires matplotlib)
   eset.plot(kind="heatmap")

   # Tanimoto similarity matrix
   sim = eset.similarity_matrix()
   print(sim)


Accessing individual Embeddings
---------------------------------

.. code-block:: python

   emb = eset[0]              # first compound
   print(emb.name)            # "PFOA"
   print(emb.array.sum())     # number of set bits
   print(emb.on_bits())       # list of set bit names


Pairwise Tanimoto similarity
-----------------------------

.. code-block:: python

   from pyCSRML import Embedding

   sim = emb1.tanimoto(emb2)
   print(f"Tanimoto(PFOA, PFOS) = {sim:.3f}")
