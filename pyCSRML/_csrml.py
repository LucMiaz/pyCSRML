"""
Parse CSRML (Chemical Subgraph Representation Markup Language) XML files
(ToxPrint v2 and TxP_PFAS v1) and convert subgraph patterns to SMARTS strings.

CSRML versions handled:
  - csrmlVersion="2" (ToxPrint v2.0): substructureMatch + substructureException molecules
  - csrmlVersion="1" (TxP_PFAS v1.0): single substructureMatch + mustMatch/mustNotMatch test cases

Key simplifications:
  - substructureException patterns are approximated or skipped (→ minor false positives)
  - `matchingQueryAtom` cross-references between exception and main molecules are ignored
  - Complex atom descriptors (atomDescriptorValue/Range, combineAtomFeatures, elementGroup)
    evaluate to wildcard `*`
  - `query` bonds (CSRML bondList OR-condition) map to `~` (any bond in SMARTS)
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import Optional

NS = {"csrml": "http://www.molecular-networks.com/schema/csrml"}

# ---------------------------------------------------------------------------
# Periodic table lookup (symbol → atomic number) using RDKit when available,
# else a static fallback table.
# ---------------------------------------------------------------------------
try:
    from rdkit.Chem import GetPeriodicTable as _gpt

    _PT = _gpt()

    def _atomic_num(symbol: str) -> Optional[int]:
        try:
            return _PT.GetAtomicNumber(symbol)
        except Exception:
            return None

except ImportError:
    _FALLBACK = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
        "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
        "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
        "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
        "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
        "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
        "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
        "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
    }

    def _atomic_num(symbol: str) -> Optional[int]:
        return _FALLBACK.get(symbol)


# ---------------------------------------------------------------------------
# XML helpers
# ---------------------------------------------------------------------------

def _tag(el) -> str:
    """Return local tag name without namespace prefix."""
    tag = el.tag
    return tag.split("}", 1)[1] if "}" in tag else tag


def _text(el, path: str, default: str = "") -> str:
    node = el.find(path, NS)
    return (node.text or "").strip() if node is not None else default


# ---------------------------------------------------------------------------
# Atom → SMARTS primitive
# ---------------------------------------------------------------------------

# Special CSRML element pseudo-symbols that do not correspond to real atoms
_SPECIAL_ELEMENTS = {
    # wildcard / any atom
    "*": "*",
    # heteroatom (non-C, non-H) – SMARTS: !#6;!#1
    "Q": "[!#6;!#1]",
    # any halogen
    "X": "[#9,#17,#35,#53]",  # F, Cl, Br, I
    # explicitly used as "any atom except"
    "G": "[!#6]",  # generic; context-dependent – used for metals/metalloids
}

_BOND_SMARTS = {
    "single": "-",
    "double": "=",
    "triple": "#",
    "aromatic": ":",
    "any": "~",
    "query": "~",  # CSRML bondList OR-condition → approximate with any
}


def _parse_atom_matchifs(atom_el) -> dict:
    """
    Return a dict of matchIf feature → data for a CSRML atom element.
    Data is a dict with keys: 'values' (list of str), and attribute keys
    (value, min, max, count, id) if present.
    """
    result = {}
    for mif in atom_el.findall("csrml:matchIf", NS):
        feat = mif.get("feature", "")
        if not feat:
            continue
        data = {}
        for attr in ("value", "min", "max", "count", "id", "combineBy"):
            v = mif.get(attr)
            if v is not None:
                data[attr] = v
        vals = [v.text.strip() for v in mif.findall("csrml:value", NS) if v.text]
        if vals:
            data["values"] = vals
        result[feat] = data
    return result


def _atom_to_smarts(element: Optional[str], matchifs: dict) -> str:
    """
    Convert a CSRML atom (element + matchIf dict) to a SMARTS atom string.
    Returns a bracket-form like `[#6;H1]` or a simple form like `C`.
    """
    element = element or ""
    parts_and: list[str] = []  # AND-conditions
    parts_or: list[str] = []   # OR-conditions (atomList)

    # --- Base element ---
    if element == "*" or not element:
        parts_and.append("*")
    elif element in _SPECIAL_ELEMENTS:
        # Return early for simple special elements (they're pre-formatted)
        special = _SPECIAL_ELEMENTS[element]
        # We'll still apply aromaticity / Hcount etc. if needed
        if special == "*":
            parts_and.append("*")
        else:
            # Already a bracket expression; if extra conditions, they're ignored
            # (these are rare and complex; return as-is)
            return special
    elif element == "QRY":
        # Defined entirely by matchIf; start with nothing (wildcard)
        pass
    else:
        n = _atomic_num(element)
        if n is not None:
            parts_and.append(f"#{n}")
        else:
            parts_and.append("*")  # unknown element → wildcard

    # --- Aromaticity ---
    if "aromaticAtom" in matchifs:
        parts_and.append("a")
    elif "aliphaticAtom" in matchifs:
        parts_and.append("A")

    # --- Hydrogen count ---
    if "attachedHydrogenCount" in matchifs:
        d = matchifs["attachedHydrogenCount"]
        hval = d.get("count") or d.get("value") or (d.get("values") or [None])[0]
        if hval is not None:
            parts_and.append(f"H{hval}")

    # --- Connectivity (total degree) ---
    if "connectivity" in matchifs:
        d = matchifs["connectivity"]
        cval = d.get("count") or d.get("value") or (d.get("values") or [None])[0]
        if cval is not None:
            parts_and.append(f"X{cval}")

    # --- Formal charge ---
    if "atomicFormalCharge" in matchifs:
        d = matchifs["atomicFormalCharge"]
        chval = d.get("value") or (d.get("values") or [None])[0]
        if chval is not None:
            c = int(chval)
            if c > 0:
                parts_and.append(f"+{c}")
            elif c < 0:
                parts_and.append(str(c))

    # --- Ring membership ---
    if "ringAtom" in matchifs:
        parts_and.append("R")
    elif "chainAtom" in matchifs:
        parts_and.append("!R")

    # --- Valency ---
    if "valency" in matchifs:
        d = matchifs["valency"]
        vval = d.get("value") or d.get("count") or (d.get("values") or [None])[0]
        if vval is not None:
            parts_and.append(f"v{vval}")

    # --- Heteroatom attached count ---
    if "atomHeteroAttachedCount" in matchifs:
        d = matchifs["atomHeteroAttachedCount"]
        hac = d.get("count") or d.get("value") or (d.get("values") or [None])[0]
        # In SMARTS no direct equivalent; skip

    # --- atomList (OR of elements) ---
    if "atomList" in matchifs:
        vals = matchifs["atomList"].get("values", [])
        for v in vals:
            if v.startswith("#"):
                parts_or.append(v)
            else:
                # Some values have SMARTS-like notation (e.g., '#8H1' = O-H)
                if any(c.isalpha() for c in v[1:]):
                    # Complex value; skip
                    continue
                n = _atomic_num(v)
                if n is not None:
                    parts_or.append(f"#{n}")

    # --- excludeAtomList (NOT of elements) ---
    if "excludeAtomList" in matchifs:
        vals = matchifs["excludeAtomList"].get("values", [])
        for v in vals:
            n = _atomic_num(v)
            if n is not None:
                parts_and.append(f"!#{n}")

    # --- Assemble ---
    # If only wildcard, return simple form
    if parts_and == ["*"] and not parts_or:
        return "*"

    # Build OR part (atomList) and AND part together
    cond_parts = parts_or  # goes first (OR-ed)
    if cond_parts:
        or_str = ",".join(cond_parts)
    else:
        or_str = ""

    # Filter out pure wildcard from and-parts
    and_parts = [p for p in parts_and if p != "*"]

    if not and_parts and not or_str:
        return "*"

    # If only one plain atomic number and nothing else, simplify
    if not or_str and len(and_parts) == 1 and and_parts[0].startswith("#"):
        # Can we use a simple symbol?
        num = int(and_parts[0][1:])
        # Use [#n] bracket form to be safe
        return f"[{and_parts[0]}]"

    # Combine: OR-list (comma-joined) + AND conditions (semicolon-joined)
    combined: list[str] = []
    if or_str:
        combined.append(or_str)
    combined.extend(and_parts)

    return "[" + ";".join(combined) + "]"


# ---------------------------------------------------------------------------
# Graph-to-SMARTS DFS builder
# ---------------------------------------------------------------------------

def _bond_char(order: str) -> str:
    """Return the SMARTS bond character for a CSRML bond order."""
    return _BOND_SMARTS.get(order, "~")


def _graph_to_smarts(atoms: dict, bonds: list[tuple]) -> Optional[str]:
    """
    Build a SMARTS string for a small query graph.

    Parameters
    ----------
    atoms : dict
        Mapping atom_id → {'element': str, 'matchifs': dict}
    bonds : list of (atom_id1, atom_id2, order_str)
    """
    if not atoms:
        return None

    # Build undirected adjacency list
    adj: dict[str, list[tuple]] = defaultdict(list)
    for a1, a2, order in bonds:
        adj[a1].append((a2, order))
        adj[a2].append((a1, order))

    atom_ids = list(atoms.keys())
    start = atom_ids[0]

    # ---- Two-pass ring-closure detection ----
    # Pass 1: DFS to find tree edges and back edges
    visited_p1: set = set()
    tree_edges: set = set()  # (parent, child)
    back_edges: list = []    # (desc, anc, order)

    def _find_back_edges(node: str, prev: Optional[str]):
        visited_p1.add(node)
        for nb, order in adj[node]:
            if nb not in visited_p1:
                tree_edges.add((node, nb))
                tree_edges.add((nb, node))  # bidirectional for lookup
                _find_back_edges(nb, node)
            elif nb != prev:
                back_edges.append((node, nb, order))

    _find_back_edges(start, None)

    # Assign ring closure numbers — process each (unordered) pair exactly once.
    # In the DFS, back edges can appear in BOTH directions, e.g. (a→b) and (b→a),
    # because the start node has no parent. We deduplicate by normalised key.
    rc_counter = [1]
    rc_map: dict[tuple, int] = {}

    # Build ring-closure annotation per atom.
    # Convention: the ANCESTOR atom (appears first in the SMARTS string) carries the
    # bond-type specifier; the DESCENDANT (appears second) carries just the number.
    atom_rc: dict[str, list[str]] = defaultdict(list)
    seen_rc_keys: set = set()

    for desc, anc, order in back_edges:
        key = (min(desc, anc), max(desc, anc))
        if key in seen_rc_keys:
            continue  # skip the reverse-direction duplicate
        seen_rc_keys.add(key)

        n = rc_counter[0]
        rc_map[key] = n
        rc_counter[0] += 1
        num_str = f"%{n:02d}" if n >= 10 else str(n)

        bond_c = _bond_char(order)
        # Omit explicit "-" for single bonds to produce cleaner SMARTS.
        if bond_c == "-":
            bond_c = ""

        # anc is visited first in DFS → appears first in SMARTS → takes bond spec.
        atom_rc[anc].append(f"{bond_c}{num_str}")
        # desc is visited second → just the closing number.
        atom_rc[desc].append(num_str)

    # ---- Pass 2: DFS to build SMARTS string ----
    visited_p2: set = set()

    def _dfs_smarts(node: str, prev: Optional[str]) -> str:
        visited_p2.add(node)
        atom_info = atoms[node]
        base = _atom_to_smarts(atom_info["element"], atom_info["matchifs"])
        rc_suffix = "".join(atom_rc.get(node, []))

        # Only follow tree edges downward: neighbor must be unvisited, not the
        # parent, AND the edge (node, nb) must be a tree edge from pass 1.
        children = [
            (nb, ord_)
            for (nb, ord_) in adj[node]
            if nb != prev and nb not in visited_p2 and (node, nb) in tree_edges
        ]

        if not children:
            return base + rc_suffix

        child_strs = []
        for nb, ord_ in children:
            bond_c = _bond_char(ord_)
            child_smarts = _dfs_smarts(nb, node)
            child_strs.append(f"{bond_c}{child_smarts}")

        # Last child is the "main chain"; all others are branches
        main_child = child_strs[-1]
        branch_children = child_strs[:-1]

        branches = "".join(f"({b})" for b in branch_children)
        return base + rc_suffix + branches + main_child

    try:
        return _dfs_smarts(start, None)
    except RecursionError:
        return None


# ---------------------------------------------------------------------------
# Parse a single CSRML <molecule> element
# ---------------------------------------------------------------------------

def _parse_molecule(mol_el) -> Optional[dict]:
    """
    Parse a CSRML <molecule> element into a dict with keys:
      'feature': 'substructureMatch' | 'substructureException' | None
      'atoms': {atom_id: {'element': str, 'matchifs': dict}}
      'bonds': [(a1_id, a2_id, order)]
    Returns None if the molecule has no atoms.
    """
    mif_el = mol_el.find("csrml:matchIf", NS)
    feature = mif_el.get("feature", "") if mif_el is not None else ""

    atoms: dict = {}
    for atom_el in mol_el.findall(".//csrml:atoms/csrml:atom", NS):
        aid = atom_el.get("id")
        if aid is None:
            continue
        elem = atom_el.get("element")
        matchifs = _parse_atom_matchifs(atom_el)
        atoms[aid] = {"element": elem, "matchifs": matchifs}

    if not atoms:
        return None

    bonds: list = []
    for bond_el in mol_el.findall(".//csrml:bonds/csrml:bond", NS):
        order = bond_el.get("order", "single")
        atom_refs = [a.get("id") for a in bond_el.findall("csrml:atom", NS)]
        if len(atom_refs) == 2 and all(ar in atoms for ar in atom_refs):
            bonds.append((atom_refs[0], atom_refs[1], order))

    return {"feature": feature, "atoms": atoms, "bonds": bonds}


# ---------------------------------------------------------------------------
# Parse a single <subgraph> element → metadata + SMARTS
# ---------------------------------------------------------------------------

def parse_subgraph(sg_el) -> dict:
    """
    Parse a CSRML <subgraph> element.

    Returns a dict:
      id, label, title, comment,
      smarts (str or None),
      exception_smarts (list of str)
    """
    sg_id = sg_el.get("id", "")
    label = _text(sg_el, "csrml:label")
    title = _text(sg_el, "csrml:title")
    comment = _text(sg_el, "csrml:comment")

    match_mol = None
    exception_mols: list[dict] = []

    for mol_el in sg_el.findall("csrml:molecule", NS):
        parsed = _parse_molecule(mol_el)
        if parsed is None:
            continue
        feat = parsed["feature"]
        if feat == "substructureMatch" and match_mol is None:
            match_mol = parsed
        elif feat == "substructureException":
            exception_mols.append(parsed)

    smarts: Optional[str] = None
    exception_smarts: list[str] = []

    if match_mol:
        smarts = _graph_to_smarts(match_mol["atoms"], match_mol["bonds"])

    for ex_mol in exception_mols:
        ex_smarts = _graph_to_smarts(ex_mol["atoms"], ex_mol["bonds"])
        if ex_smarts:
            exception_smarts.append(ex_smarts)

    return {
        "id": sg_id,
        "label": label,
        "title": title,
        "comment": comment,
        "smarts": smarts,
        "exception_smarts": exception_smarts,
    }


# ---------------------------------------------------------------------------
# Parse class hierarchy recursively
# ---------------------------------------------------------------------------

def _parse_classes(classes_el) -> list:
    """
    Recursively parse the <classes> / <class> hierarchy.
    Returns a list of dicts: {id, label, subgraph_refs, children}.
    """

    def _parse_class(el) -> dict:
        result = {
            "id": el.get("id", ""),
            "label": _text(el, "csrml:label"),
            "subgraph_refs": [c.get("subgraph") for c in el.findall("csrml:class", NS) if c.get("subgraph")],
            "children": [_parse_class(c) for c in el.findall("csrml:class", NS) if not c.get("subgraph")],
        }
        return result

    return [_parse_class(c) for c in classes_el.findall("csrml:class", NS)]


# ---------------------------------------------------------------------------
# Full XML → Python dict
# ---------------------------------------------------------------------------

def parse_csrml_xml(xml_path: str) -> dict:
    """
    Parse a CSRML XML file (ToxPrint or TxP_PFAS) and return a dict:

    {
      'id': str,
      'version': str,
      'csrml_version': str,
      'title': str,
      'description': str,
      'hierarchy': list,         # nested class hierarchy
      'subgraphs': list[dict],   # ordered list of parsed subgraphs
      'subgraph_index': dict,    # id → subgraph dict
    }
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    doc_id = root.get("id", "")
    csrml_version = root.get("csrmlVersion", "")
    title = _text(root, "csrml:title")
    description = _text(root, "csrml:description")

    # Class hierarchies
    hierarchy = []
    for classes_el in root.findall("csrml:classes", NS):
        classes_id = classes_el.get("id", "")
        hier = {
            "id": classes_id,
            "title": _text(classes_el, "csrml:title"),
            "classes": _parse_classes(classes_el),
        }
        hierarchy.append(hier)

    # Subgraphs
    subgraphs = []
    for sg_el in root.findall("csrml:subgraph", NS):
        sg = parse_subgraph(sg_el)
        subgraphs.append(sg)

    subgraph_index = {sg["id"]: sg for sg in subgraphs}

    return {
        "id": doc_id,
        "version": doc_id,
        "csrml_version": csrml_version,
        "title": title,
        "description": description,
        "hierarchy": hierarchy,
        "subgraphs": subgraphs,
        "subgraph_index": subgraph_index,
    }


# ---------------------------------------------------------------------------
# Build ordered fingerprint bit list from class hierarchy
# ---------------------------------------------------------------------------

def _collect_refs_from_class(cls: dict) -> list[str]:
    """Recursively collect all subgraph references in DFS order."""
    refs = list(cls["subgraph_refs"])
    for child in cls["children"]:
        refs.extend(_collect_refs_from_class(child))
    return refs


def ordered_bit_list(parsed: dict) -> list[str]:
    """
    Return the ordered list of subgraph IDs (fingerprint bit order)
    following the class hierarchy.

    Falls back to the subgraphs list order if no hierarchy is present.
    """
    refs: list[str] = []
    for tier in parsed["hierarchy"]:
        for cls in tier["classes"]:
            refs.extend(_collect_refs_from_class(cls))

    # Remove duplicates while preserving order
    seen: set = set()
    ordered: list[str] = []
    for r in refs:
        if r not in seen:
            seen.add(r)
            ordered.append(r)

    # Any subgraphs not referenced in hierarchy → append at end
    for sg in parsed["subgraphs"]:
        if sg["id"] not in seen:
            ordered.append(sg["id"])

    return ordered
