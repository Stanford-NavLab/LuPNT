import os
import sys
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parent
REPO_ROOT = DOCS_DIR.parent
os.makedirs(DOCS_DIR / "_build", exist_ok=True)

sys.path.insert(0, str(REPO_ROOT / "python"))

project = "LuPNT"
author = "Stanford NAV Lab"
copyright = "2026, Stanford NAV Lab"
version = "0.1.0"
release = version

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "breathe",
    "exhale",
    "nbsphinx",
    "nbsphinx_link",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_title = "LuPNT"
html_static_path = ["_static"]
html_favicon = "_static/lupnt_mark.svg"
html_theme_options = {
    "light_logo": "lupnt_logo_horizontal.svg",
    "dark_logo": "lupnt_logo_horizontal_dark.svg",
}

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "inherited-members": True,
}
autosummary_generate = True
napoleon_google_docstring = True
napoleon_numpy_docstring = True
nbsphinx_execute = "never"

breathe_projects = {"lupnt": str(DOCS_DIR / "_build" / "doxygen" / "xml")}
breathe_default_project = "lupnt"

exhale_args = {
    "containmentFolder": "./cpp_api",
    "rootFileName": "cpp_library_root.rst",
    "rootFileTitle": "C++ API",
    "doxygenStripFromPath": str(REPO_ROOT / "cpp" / "lupnt"),
    "createTreeView": True,
    "exhaleExecutesDoxygen": True,
    "exhaleUseDoxyfile": False,
    "exhaleDoxygenStdin": f"""
PROJECT_NAME = LuPNT
INPUT = {REPO_ROOT / "cpp" / "lupnt"}
RECURSIVE = YES
# Document only public headers. Implementation files (*.cc) hold internal
# anonymous-namespace helpers that Doxygen would otherwise surface as
# `lupnt::@<number>` namespaces in the API; the public API lives in the headers.
FILE_PATTERNS = *.h
EXCLUDE_PATTERNS = */environment/plasma/fortran/*
GENERATE_XML = YES
GENERATE_HTML = NO
GENERATE_LATEX = NO
XML_OUTPUT = xml
EXTRACT_ALL = YES
EXTRACT_ANON_NSPACES = NO
HIDE_UNDOC_MEMBERS = NO
HIDE_UNDOC_CLASSES = NO
QUIET = YES
""",
}

primary_domain = "cpp"
highlight_language = "cpp"
