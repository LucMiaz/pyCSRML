CSRML format
============

The Chemical Subgraph Representation Markup Language (CSRML) is an XML-based
language for encoding molecular subgraph patterns.
pyCSRML parses the two fingerprint definition files shipped with the
`ChemoTyper <https://github.com/mn-am/chemotyper>`_ tool:
ToxPrint v2.0 (729 bits) and TxP_PFAS v1.0.4 (129 bits).


File structure
--------------

A CSRML file contains:

.. code-block:: xml

   <csrml id="..." version="..." csrmlVersion="2">
     <title>...</title>
     <description>...</description>
     <subgraphs>
       <subgraph id="txp-001" label="bond:CC_...">
         <subgraphMolecule feature="substructureMatch">
           <molecule feature="substructureMatch">
             <atomArray> ... </atomArray>
             <bondArray> ... </bondArray>
           </molecule>
           <molecule feature="substructureException"> ... </molecule>
         </subgraphMolecule>
       </subgraph>
       ...
     </subgraphs>
     <classes> ... </classes>
   </csrml>

Each ``<subgraph>`` defines one fingerprint bit. The ``label`` attribute becomes
the bit name.


Atom features (matchIf)
------------------------

Atoms within a CSRML pattern carry ``<matchIf>`` child elements that refine the
query beyond a simple element match:

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - CSRML feature
     - SMARTS mapping
   * - ``atomList`` (element list)
     - ``[#6,#8]``  (OR of atomic numbers)
   * - ``atomList negate="true"``
     - ``[!#9]``  (negated list → AND of ``!#n``)
   * - ``aliphaticAtom``
     - ``A``
   * - ``aromaticAtom``
     - ``a``
   * - ``attachedHydrogenCount value=1``
     - ``H1``
   * - ``attachedHydrogenCount`` (minInclusive / minExclusive)
     - ``!H0``, ``!H0;!H1``, …
   * - ``connectivity value=2``
     - ``X2``
   * - ``ringAtom``
     - ``R``
   * - ``chainAtom``
     - ``!R``
   * - ``ringCountAtom``
     - ``R``
   * - ``atomicFormalCharge value=-1``
     - ``-1``
   * - ``combineAtomFeatures combineBy="or"``
     - recursive OR-of-AND tree (see below)


Pseudo-elements
---------------

CSRML defines several symbolic element codes that are not in the periodic table:

.. list-table::
   :widths: 15 20 65
   :header-rows: 1

   * - Symbol
     - SMARTS
     - Meaning
   * - ``*``
     - ``*``
     - any atom (wildcard)
   * - ``Q``
     - ``[!#6;!#1]``
     - heteroatom (non-C, non-H)
   * - ``Z``
     - ``[!#6;!#1]``
     - heteroatom (synonym of Q in TxP_PFAS)
   * - ``X``
     - ``[#9,#17,#35,#53]``
     - any halogen (F, Cl, Br, I)
   * - ``G``
     - ``[#5,#14,#32,#33,#51,#52]``
     - metal / metalloid (B, Si, Ge, As, Sb, Te)
   * - ``QRY``
     - derived from ``combineAtomFeatures``
     - complex query atom (see below)


combineAtomFeatures (QRY atoms)
--------------------------------

Some atoms (labelled ``QRY``) use a recursive ``combineAtomFeatures`` structure
that encodes an OR-of-AND tree:

.. code-block:: xml

   <matchIf feature="combineAtomFeatures" combineBy="or">
     <matchIf feature="combineAtomFeatures" combineBy="and">
       <matchIf feature="atomList" value="C"/>
       <matchIf feature="aliphaticAtom"/>
       <matchIf feature="connectivity" value="3"/>
     </matchIf>
     <matchIf feature="combineAtomFeatures" combineBy="and">
       <matchIf feature="atomList" value="C"/>
       <matchIf feature="aliphaticAtom"/>
       <matchIf feature="attachedHydrogenCount" minInclusive="1"/>
     </matchIf>
     <!-- ... -->
   </matchIf>

This maps to the SMARTS atom primitive: ``[#6;A;X3,#6;A;!H0,...]``.


Exception handling
------------------

Fingerprint bits may have one or more ``substructureException`` molecules.
An exception molecule can either:

1. **Stand-alone** – the bit is set to 0 if the molecule also matches the exception
   as a global substructure query.

2. **Anchor-linked** (``matchingQueryAtom``) – one atom in the exception has a
   ``matchingQueryAtom`` feature referencing an atom in the main pattern.  In this
   case the exception is folded into the main SMARTS as a ``[!$(...)]`` recursive
   negation on the referenced atom, encoding the constraint directly rather than
   as a separate pass.

Example: the bit ``pfas_chain:perF-linear_cap_C1_excl`` has the main SMARTS
``[#6](-[#9])(-[#9])(-[#9])-[…]`` and an exception anchored at the terminal CF₃.
After folding the result is::

   [#6](-[#9])(-[#9])(-[#9])-[!$([#6;A](-[#9])-[#9])]


Low-level API
-------------

.. code-block:: python

   from pyCSRML._csrml import parse_csrml_xml, ordered_bit_list

   data = parse_csrml_xml("TxP_PFAS_v1.0.4.xml")
   bits = ordered_bit_list(data)

   for bit in bits[:3]:
       print(bit["id"], "→", bit["smarts"])
       for exc in bit["exception_smarts"]:
           print("  excl:", exc)


Known limitations
-----------------

The following CSRML features are recognised in the XML but cannot be fully
represented as SMARTS by pyCSRML.

``elSystems`` / ``piSystem`` / ``piElectronCount``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CSRML v1 uses ``<elSystems>`` elements to constrain atoms by membership in a
conjugated or aromatic pi-electron system (e.g. ``piElectronCount ≥ 4`` for
the sulfonyl O=S=O group).  These elements have no direct SMARTS equivalent
and are silently ignored by the parser.  In practice the explicit double-bond
topology already implied by the pattern graph produces the same constraint for
all tested compounds, so the impact on accuracy is negligible.

Fluorotelomer (``FT_n*``) bits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several fluorotelomer chain-length bits (``pfas_chain:FT_n1_*``,
``FT_n2_*``, ``FT_n3_*``) achieve 93–99 % concordance rather than 100 %.
The ChemoTyper reference implementation uses Java-side chain-length counting
that cannot be reproduced exactly with a single SMARTS query.  The patterns
generated by pyCSRML are therefore conservative approximations.

Concordance summary
~~~~~~~~~~~~~~~~~~~

After all implemented fixes the overall pyCSRML accuracy across 14 710 PFAS
compounds from Richard *et al.* (2023) is **99.44 %**, with **all 129 bits
reaching ≥ 90 % accuracy**.  The worst single bit is
``pfas_chain:FT_n2_hetero`` at 93.7 %.
