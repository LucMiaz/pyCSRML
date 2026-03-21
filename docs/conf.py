# Sphinx configuration for pyCSRML docs
import os
import sys

# Allow Sphinx to find the package source
sys.path.insert(0, os.path.abspath(".."))

# ---------------------------------------------------------------------------
# Project information
# ---------------------------------------------------------------------------
project = "pyCSRML"
copyright = "2026, Luc Miaz"
author = "Luc Miaz"
release = "0.1.0"

# ---------------------------------------------------------------------------
# Extensions
# ---------------------------------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",          # Google / NumPy-style docstrings
    "sphinx.ext.viewcode",          # [source] links
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "myst_parser",                  # Markdown support
]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "rdkit": ("https://www.rdkit.org/docs/", None),
}

# ---------------------------------------------------------------------------
# File handling
# ---------------------------------------------------------------------------
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 4,
    "titles_only": False,
}
html_static_path = ["_static"]
