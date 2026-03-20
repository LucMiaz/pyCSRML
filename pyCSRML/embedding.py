"""
Embedding and EmbeddingSet: containers for ToxPrint / TxP_PFAS fingerprint results.

Usage::

    from pyToxPrint.fingerprinter import ToxPrintFingerprinter, PFASFingerprinter
    from pyToxPrint.embedding import Embedding, EmbeddingSet
    from rdkit import Chem

    fp = PFASFingerprinter()

    # Single embedding
    mol = Chem.MolFromSmiles("FCCCF")
    arr, names = fp.fingerprint(mol)
    emb = Embedding(
        array=arr,
        bit_names=names,
        smiles="FCCCF",
        dtxsid="DTXSID1234",
        inchi="InChI=1S/...",
        inchikey="XXXX-XX",
        formula="C3H4F2",
        name="Sample",
        fingerprint_type="TxP_PFAS",
    )
    emb.plot()

    # Set of embeddings
    eset = EmbeddingSet([emb1, emb2, emb3])
    eset.plot(kind="heatmap")
    eset.umap()
    eset.cluster(method="hierarchical", n_clusters=3)
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Iterable, Optional, Union

import numpy as np

# Optional heavy deps — imported lazily so the module loads without them.
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False


# ---------------------------------------------------------------------------
# Similarity metrics
# ---------------------------------------------------------------------------

def _tanimoto(a: np.ndarray, b: np.ndarray) -> float:
    """Binary Tanimoto / Jaccard coefficient."""
    a_and_b = int(np.dot(a, b))
    a_or_b = int(np.sum(a)) + int(np.sum(b)) - a_and_b
    return a_and_b / a_or_b if a_or_b > 0 else 0.0


def _dice(a: np.ndarray, b: np.ndarray) -> float:
    """Dice coefficient (Sorensen–Dice)."""
    a_and_b = int(np.dot(a, b))
    denom = int(np.sum(a)) + int(np.sum(b))
    return 2 * a_and_b / denom if denom > 0 else 0.0


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else 0.0


_METRICS = {"tanimoto": _tanimoto, "dice": _dice, "cosine": _cosine}


# ---------------------------------------------------------------------------
# Embedding
# ---------------------------------------------------------------------------

@dataclass
class Embedding:
    """
    A single chemotype fingerprint with associated compound metadata.

    Parameters
    ----------
    array : np.ndarray
        Binary fingerprint vector (dtype bool or int).
    bit_names : list[str]
        Ordered list of chemotype labels matching the array.
    smiles : str, optional
    dtxsid : str, optional
        DSSTox Substance Identifier.
    inchi : str, optional
    inchikey : str, optional
    formula : str, optional
    name : str, optional
        Human-readable compound name.
    fingerprint_type : str, optional
        E.g. "ToxPrint", "TxP_PFAS".
    metadata : dict, optional
        Any extra key-value pairs.
    """

    array: np.ndarray
    bit_names: list[str]
    smiles: str = ""
    dtxsid: str = ""
    inchi: str = ""
    inchikey: str = ""
    formula: str = ""
    name: str = ""
    fingerprint_type: str = ""
    metadata: dict = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def n_bits(self) -> int:
        return len(self.array)

    @property
    def on_bits(self) -> np.ndarray:
        """Indices of set (1) bits."""
        return np.where(self.array)[0]

    @property
    def on_names(self) -> list[str]:
        """Labels of set bits."""
        return [self.bit_names[i] for i in self.on_bits]

    @property
    def density(self) -> float:
        """Fraction of bits that are set."""
        return float(np.mean(self.array))

    def __repr__(self) -> str:
        label = self.name or self.smiles or self.dtxsid or "?"
        return (
            f"Embedding({label!r}, type={self.fingerprint_type!r}, "
            f"bits={self.n_bits}, on={int(np.sum(self.array))})"
        )

    # ------------------------------------------------------------------
    # Similarity
    # ------------------------------------------------------------------

    def similarity(
        self,
        other: "Embedding",
        metric: str = "tanimoto",
    ) -> float:
        """
        Compute similarity to another Embedding.

        Parameters
        ----------
        other : Embedding
        metric : {'tanimoto', 'dice', 'cosine'}
        """
        fn = _METRICS.get(metric)
        if fn is None:
            raise ValueError(f"Unknown metric {metric!r}. Choose from {list(_METRICS)}")
        a = self.array.astype(float)
        b = other.array.astype(float)
        if len(a) != len(b):
            raise ValueError(
                f"Array length mismatch: {len(a)} vs {len(b)}. "
                "Both embeddings must use the same fingerprint."
            )
        return fn(a, b)

    # ------------------------------------------------------------------
    # I/O helpers
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Convert to a serialisable dict (array as list)."""
        return {
            "smiles": self.smiles,
            "dtxsid": self.dtxsid,
            "inchi": self.inchi,
            "inchikey": self.inchikey,
            "formula": self.formula,
            "name": self.name,
            "fingerprint_type": self.fingerprint_type,
            "bit_names": self.bit_names,
            "array": self.array.tolist(),
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Embedding":
        arr = np.array(d["array"], dtype=bool)
        return cls(
            array=arr,
            bit_names=d.get("bit_names", []),
            smiles=d.get("smiles", ""),
            dtxsid=d.get("dtxsid", ""),
            inchi=d.get("inchi", ""),
            inchikey=d.get("inchikey", ""),
            formula=d.get("formula", ""),
            name=d.get("name", ""),
            fingerprint_type=d.get("fingerprint_type", ""),
            metadata=d.get("metadata", {}),
        )

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------

    def plot(
        self,
        kind: str = "bar",
        max_bits: Optional[int] = None,
        figsize: tuple = (10, 3),
        title: Optional[str] = None,
        ax=None,
    ):
        """
        Visualise the fingerprint.

        Parameters
        ----------
        kind : {'bar', 'bits'}
            'bar'  – horizontal bar chart of set bits (bit-on count / label).
            'bits' – bitmap-style row of coloured squares.
        max_bits : int, optional
            Limit the number of displayed bits (for 'bar': most-informative).
        figsize : tuple
        title : str, optional
        ax : matplotlib Axes, optional
        """
        if not _HAS_MPL:
            raise ImportError("matplotlib is required for plotting.")

        label = title or self.name or self.smiles or self.dtxsid or "Embedding"
        on = self.on_bits
        on_labels = [self.bit_names[i] for i in on]

        if kind == "bar":
            if not on_labels:
                on_labels = ["(no bits set)"]
                vals = [0]
            else:
                vals = [1] * len(on_labels)
            if max_bits and len(on_labels) > max_bits:
                on_labels = on_labels[:max_bits]
                vals = vals[:max_bits]

            create_ax = ax is None
            if create_ax:
                fig, ax = plt.subplots(figsize=figsize)
            ax.barh(on_labels, vals, color="steelblue")
            ax.set_xlabel("Chemotype active")
            ax.set_title(label)
            ax.invert_yaxis()
            if create_ax:
                plt.tight_layout()
                plt.show()

        elif kind == "bits":
            n = max_bits or self.n_bits
            arr = self.array[:n].astype(float).reshape(1, -1)
            create_ax = ax is None
            if create_ax:
                fig, ax = plt.subplots(figsize=figsize)
            ax.imshow(arr, aspect="auto", cmap="Blues", vmin=0, vmax=1)
            ax.set_yticks([])
            ax.set_xlabel("Bit index")
            ax.set_title(label)
            if create_ax:
                plt.tight_layout()
                plt.show()
        else:
            raise ValueError(f"Unknown kind {kind!r}. Choose 'bar' or 'bits'.")


# ---------------------------------------------------------------------------
# EmbeddingSet
# ---------------------------------------------------------------------------

class EmbeddingSet:
    """
    A collection of Embedding objects with batch analysis tools.

    Parameters
    ----------
    embeddings : iterable of Embedding
    """

    def __init__(self, embeddings: Iterable[Embedding]):
        self.embeddings: list[Embedding] = list(embeddings)
        if not self.embeddings:
            raise ValueError("EmbeddingSet must contain at least one Embedding.")

        # Validate consistent fingerprints
        n = self.embeddings[0].n_bits
        if not all(e.n_bits == n for e in self.embeddings):
            raise ValueError("All embeddings must have the same number of bits.")
        self._n_bits = n

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self.embeddings)

    def __iter__(self):
        return iter(self.embeddings)

    def __getitem__(self, idx):
        return self.embeddings[idx]

    @property
    def n_bits(self) -> int:
        return self._n_bits

    @property
    def bit_names(self) -> list[str]:
        return self.embeddings[0].bit_names

    @property
    def labels(self) -> list[str]:
        """Best available label for each embedding."""
        return [e.name or e.smiles or e.dtxsid or f"#{i}" for i, e in enumerate(self.embeddings)]

    # ------------------------------------------------------------------
    # Matrix helpers
    # ------------------------------------------------------------------

    def to_matrix(self) -> np.ndarray:
        """Return (n_compounds, n_bits) boolean matrix."""
        return np.vstack([e.array for e in self.embeddings]).astype(bool)

    def to_dataframe(self, include_metadata: bool = True):
        """
        Return a pandas DataFrame with one row per compound.

        Columns: smiles, dtxsid, inchi, inchikey, formula, name,
        fingerprint_type (if include_metadata), then one column per bit.
        """
        if not _HAS_PANDAS:
            raise ImportError("pandas is required for to_dataframe().")
        meta_cols = ["name", "smiles", "dtxsid", "inchi", "inchikey", "formula", "fingerprint_type"]
        meta = {col: [getattr(e, col) for e in self.embeddings] for col in meta_cols}
        bits = {name: self.to_matrix()[:, i].tolist() for i, name in enumerate(self.bit_names)}
        if include_metadata:
            data = {**meta, **bits}
        else:
            data = bits
        return pd.DataFrame(data)

    # ------------------------------------------------------------------
    # Similarity
    # ------------------------------------------------------------------

    def similarity_matrix(self, metric: str = "tanimoto") -> np.ndarray:
        """
        Compute the pairwise similarity matrix.

        Returns
        -------
        S : np.ndarray of shape (n, n), dtype float64
        """
        fn = _METRICS.get(metric)
        if fn is None:
            raise ValueError(f"Unknown metric {metric!r}. Choose from {list(_METRICS)}")
        n = len(self.embeddings)
        S = np.zeros((n, n), dtype=np.float64)
        mat = self.to_matrix().astype(float)
        for i in range(n):
            for j in range(i, n):
                s = fn(mat[i], mat[j])
                S[i, j] = s
                S[j, i] = s
        return S

    # ------------------------------------------------------------------
    # Bit statistics
    # ------------------------------------------------------------------

    def bit_frequency(self) -> np.ndarray:
        """Return fraction of compounds having each bit set."""
        return self.to_matrix().mean(axis=0)

    # ------------------------------------------------------------------
    # Dimensionality reduction
    # ------------------------------------------------------------------

    def pca(
        self,
        n_components: int = 2,
        return_model: bool = False,
    ) -> Union[np.ndarray, tuple]:
        """
        PCA on the binary fingerprint matrix.

        Returns
        -------
        coords : np.ndarray of shape (n, n_components)
        If return_model is True: (coords, sklearn PCA object)
        """
        try:
            from sklearn.decomposition import PCA
        except ImportError:
            raise ImportError("scikit-learn is required for PCA.")
        mat = self.to_matrix().astype(float)
        pca_model = PCA(n_components=n_components)
        coords = pca_model.fit_transform(mat)
        return (coords, pca_model) if return_model else coords

    def umap(
        self,
        n_components: int = 2,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        metric: str = "jaccard",
        random_state: int = 42,
        return_model: bool = False,
        **kwargs,
    ) -> Union[np.ndarray, tuple]:
        """
        UMAP dimensionality reduction on the binary fingerprint matrix.

        Parameters
        ----------
        n_components : int, default 2
        n_neighbors : int, default 15
        min_dist : float, default 0.1
        metric : str, default 'jaccard' (appropriate for binary fingerprints)
        random_state : int, default 42
        return_model : bool
            If True, return (coords, umap model).
        **kwargs
            Extra keyword arguments forwarded to umap.UMAP().

        Returns
        -------
        coords : np.ndarray of shape (n, n_components)
        """
        try:
            import umap as _umap
        except ImportError:
            raise ImportError(
                "The 'umap-learn' package is required for UMAP. "
                "Install it with: pip install umap-learn"
            )
        mat = self.to_matrix().astype(float)
        reducer = _umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=random_state,
            **kwargs,
        )
        coords = reducer.fit_transform(mat)
        return (coords, reducer) if return_model else coords

    # ------------------------------------------------------------------
    # Clustering
    # ------------------------------------------------------------------

    def cluster(
        self,
        method: str = "hierarchical",
        n_clusters: Optional[int] = None,
        similarity_metric: str = "tanimoto",
        linkage: str = "average",
        return_labels: bool = True,
    ) -> np.ndarray:
        """
        Cluster the embeddings.

        Parameters
        ----------
        method : {'hierarchical', 'kmeans', 'dbscan'}
        n_clusters : int, optional
            Required for 'hierarchical' and 'kmeans'.
        similarity_metric : str
            Used for hierarchical clustering (distance = 1 − similarity).
        linkage : str
            Linkage for hierarchical clustering: 'average', 'complete', 'ward', etc.
        return_labels : bool
            If True (default) return cluster label array.

        Returns
        -------
        labels : np.ndarray of shape (n,) with cluster indices (0-based)
        """
        try:
            from sklearn.cluster import AgglomerativeClustering, KMeans, DBSCAN
        except ImportError:
            raise ImportError("scikit-learn is required for clustering.")

        mat = self.to_matrix().astype(float)

        if method == "hierarchical":
            if n_clusters is None:
                raise ValueError("n_clusters is required for hierarchical clustering.")
            dist = 1 - self.similarity_matrix(similarity_metric)
            clust = AgglomerativeClustering(
                n_clusters=n_clusters,
                metric="precomputed",
                linkage=linkage,
            )
            labels = clust.fit_predict(dist)

        elif method == "kmeans":
            if n_clusters is None:
                raise ValueError("n_clusters is required for kmeans.")
            clust = KMeans(n_clusters=n_clusters, random_state=42, n_init="auto")
            labels = clust.fit_predict(mat)

        elif method == "dbscan":
            dist = 1 - self.similarity_matrix(similarity_metric)
            clust = DBSCAN(metric="precomputed")
            labels = clust.fit_predict(dist)

        else:
            raise ValueError(f"Unknown method {method!r}. Choose 'hierarchical', 'kmeans', or 'dbscan'.")

        return labels

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------

    def plot(
        self,
        kind: str = "heatmap",
        metric: str = "tanimoto",
        top_n_bits: Optional[int] = 50,
        labels: Optional[list[str]] = None,
        figsize: Optional[tuple] = None,
        title: Optional[str] = None,
        cmap: str = "Blues",
        cluster_rows: bool = False,
        ax=None,
        **kwargs,
    ):
        """
        Visualise the EmbeddingSet.

        Parameters
        ----------
        kind : {'heatmap', 'similarity', 'bits', 'umap', 'pca', 'distribution'}
            'heatmap'    – compound × bit presence/absence matrix (top_n_bits)
            'similarity' – compound × compound similarity heatmap
            'bits'       – stacked bit diagram
            'umap'       – 2-D UMAP scatter
            'pca'        – 2-D PCA scatter
            'distribution' – histogram of pairwise similarities
        metric : str
            Similarity metric for 'similarity' and 'distribution' plots.
        top_n_bits : int, optional
            For 'heatmap': display only the top-N most frequent bits.
        labels : list[str], optional
            Compound labels for axis tick labels. Defaults to self.labels.
        figsize : tuple, optional
        title : str, optional
        cmap : str
            Matplotlib colormap name.
        cluster_rows : bool
            For 'heatmap': reorder rows by hierarchical clustering.
        ax : matplotlib Axes, optional
            Only used for simple single-Axes plots.
        **kwargs
            Forwarded to the underlying plot / UMAP call.
        """
        if not _HAS_MPL:
            raise ImportError("matplotlib is required for plotting.")

        row_labels = labels or self.labels

        if kind == "heatmap":
            _plot_heatmap(self, row_labels, top_n_bits, figsize, title, cmap, cluster_rows, ax)
        elif kind == "similarity":
            _plot_similarity(self, row_labels, metric, figsize, title, cmap, ax)
        elif kind == "bits":
            _plot_bits(self, row_labels, top_n_bits, figsize, title, cmap, ax)
        elif kind == "umap":
            _plot_embedding(self, "umap", row_labels, metric, figsize, title, **kwargs)
        elif kind == "pca":
            _plot_embedding(self, "pca", row_labels, metric, figsize, title, **kwargs)
        elif kind == "distribution":
            _plot_sim_distribution(self, metric, figsize, title, ax)
        else:
            raise ValueError(
                f"Unknown kind {kind!r}. Choose from "
                "'heatmap', 'similarity', 'bits', 'umap', 'pca', 'distribution'."
            )


# ---------------------------------------------------------------------------
# Plot helpers (internal)
# ---------------------------------------------------------------------------

def _plot_heatmap(eset, labels, top_n_bits, figsize, title, cmap, cluster_rows, ax_in):
    mat = eset.to_matrix()
    freq = mat.mean(axis=0)

    # Select top-N bits by frequency
    if top_n_bits and top_n_bits < eset.n_bits:
        top_idx = np.argsort(freq)[::-1][:top_n_bits]
    else:
        top_idx = np.arange(eset.n_bits)

    sub = mat[:, top_idx].astype(float)
    bit_lbls = [eset.bit_names[i] for i in top_idx]

    if cluster_rows:
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            from scipy.spatial.distance import pdist
            dist_vecs = pdist(sub, metric="jaccard")
            dist_vecs = np.nan_to_num(dist_vecs, nan=1.0)
            Z = linkage(dist_vecs, method="average")
            order = leaves_list(Z)
            sub = sub[order]
            labels = [labels[i] for i in order]
        except ImportError:
            pass

    n, m = sub.shape
    default_fig = (max(8, m * 0.25), max(4, n * 0.3))
    fig, ax = plt.subplots(figsize=figsize or default_fig)
    im = ax.imshow(sub.T, aspect="auto", cmap=cmap, interpolation="nearest", vmin=0, vmax=1)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(m))
    ax.set_yticklabels(bit_lbls, fontsize=6)
    ax.set_xlabel("Compound")
    ax.set_ylabel("Chemotype bit")
    ax.set_title(title or "Chemotype Fingerprint Heatmap")
    plt.colorbar(im, ax=ax, fraction=0.02, label="Active")
    plt.tight_layout()
    plt.show()


def _plot_similarity(eset, labels, metric, figsize, title, cmap, ax_in):
    S = eset.similarity_matrix(metric)
    n = len(eset)
    default_fig = (max(6, n * 0.4), max(5, n * 0.35))
    fig, ax = plt.subplots(figsize=figsize or default_fig)
    im = ax.imshow(S, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title(title or f"Pairwise {metric.capitalize()} Similarity")
    plt.colorbar(im, ax=ax, fraction=0.02, label=metric)
    plt.tight_layout()
    plt.show()


def _plot_bits(eset, labels, top_n_bits, figsize, title, cmap, ax_in):
    mat = eset.to_matrix()
    if top_n_bits:
        freq = mat.mean(axis=0)
        top_idx = np.argsort(freq)[::-1][:top_n_bits]
        mat = mat[:, top_idx]
    n, m = mat.shape
    fig, ax = plt.subplots(figsize=figsize or (max(6, m * 0.2), max(3, n * 0.5)))
    ax.imshow(mat.astype(float), aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Bit index")
    ax.set_title(title or "Chemotype Bits")
    plt.tight_layout()
    plt.show()


def _plot_embedding(eset, method, labels, metric, figsize, title, color_by=None, **kwargs):
    color_labels = color_by  # optional list of cluster/category labels

    if method == "umap":
        try:
            coords = eset.umap(**kwargs)
        except ImportError as e:
            raise e
    elif method == "pca":
        coords = eset.pca(**kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")

    fig, ax = plt.subplots(figsize=figsize or (7, 6))

    if color_labels is not None:
        unique = sorted(set(color_labels))
        cmap = plt.get_cmap("tab20")
        color_map = {u: cmap(i / max(1, len(unique) - 1)) for i, u in enumerate(unique)}
        for lbl in unique:
            mask = [i for i, c in enumerate(color_labels) if c == lbl]
            ax.scatter(coords[mask, 0], coords[mask, 1], label=str(lbl), s=40)
        ax.legend(fontsize=8, title="Cluster")
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=40, alpha=0.7)

    for i, lbl in enumerate(labels):
        ax.annotate(lbl, (coords[i, 0], coords[i, 1]), fontsize=6, alpha=0.7)

    ax.set_xlabel(f"{method.upper()} 1")
    ax.set_ylabel(f"{method.upper()} 2")
    ax.set_title(title or f"{method.upper()} of Fingerprints")
    plt.tight_layout()
    plt.show()


def _plot_sim_distribution(eset, metric, figsize, title, ax_in):
    S = eset.similarity_matrix(metric)
    n = len(eset)
    upper = S[np.triu_indices(n, k=1)]
    create_ax = ax_in is None
    fig, ax = plt.subplots(figsize=figsize or (6, 4)) if create_ax else (None, ax_in)
    ax.hist(upper, bins=30, edgecolor="white", color="steelblue", alpha=0.85)
    ax.set_xlabel(f"{metric.capitalize()} similarity")
    ax.set_ylabel("Pair count")
    ax.set_title(title or f"Pairwise {metric.capitalize()} Distribution")
    if create_ax:
        plt.tight_layout()
        plt.show()


# ---------------------------------------------------------------------------
# Factory: build EmbeddingSet from a fingerprinter + list of compounds
# ---------------------------------------------------------------------------

def from_fingerprinter(
    fingerprinter,
    mols=None,
    smiles_list: Optional[list[str]] = None,
    dtxsids: Optional[list[str]] = None,
    inchis: Optional[list[str]] = None,
    inchikeys: Optional[list[str]] = None,
    formulas: Optional[list[str]] = None,
    names: Optional[list[str]] = None,
    **metadata,
) -> EmbeddingSet:
    """
    Build an EmbeddingSet by fingerprinting a collection of compounds.

    Parameters
    ----------
    fingerprinter : Fingerprinter
        A ToxPrintFingerprinter or PFASFingerprinter instance.
    mols : list of rdkit.Chem.Mol, optional
    smiles_list : list of str, optional
        If mols is None and smiles_list is provided, mols are created from SMILES.
    dtxsids : list of str, optional
    inchis : list of str, optional
    inchikeys : list of str, optional
    formulas : list of str, optional
    names : list of str, optional
    **metadata : dict
        Additional key-value pairs stored in each Embedding.metadata.

    Returns
    -------
    EmbeddingSet
    """
    try:
        from rdkit import Chem as _Chem
    except ImportError:
        raise ImportError("RDKit is required.")

    if mols is None and smiles_list is not None:
        mols = [_Chem.MolFromSmiles(s) for s in smiles_list]
    elif mols is None:
        raise ValueError("Either mols or smiles_list must be provided.")

    n = len(mols)
    bit_names = fingerprinter.bit_names
    fp_type = fingerprinter.title or type(fingerprinter).__name__

    def _get(lst, i):
        return lst[i] if lst is not None and i < len(lst) else ""

    embeddings = []
    for i, mol in enumerate(mols):
        arr, _ = fingerprinter.fingerprint(mol)
        smi = _get(smiles_list, i)
        if not smi and mol is not None:
            try:
                smi = _Chem.MolToSmiles(mol)
            except Exception:
                smi = ""
        emb = Embedding(
            array=arr,
            bit_names=bit_names,
            smiles=smi,
            dtxsid=_get(dtxsids, i),
            inchi=_get(inchis, i),
            inchikey=_get(inchikeys, i),
            formula=_get(formulas, i),
            name=_get(names, i),
            fingerprint_type=fp_type,
            metadata=dict(metadata),
        )
        embeddings.append(emb)

    return EmbeddingSet(embeddings)
