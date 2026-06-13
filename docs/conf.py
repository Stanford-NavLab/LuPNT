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
FILE_PATTERNS = *.h *.cc
EXCLUDE_PATTERNS = */environment/plasma/fortran/*
GENERATE_XML = YES
GENERATE_HTML = NO
GENERATE_LATEX = NO
XML_OUTPUT = xml
EXTRACT_ALL = YES
HIDE_UNDOC_MEMBERS = NO
HIDE_UNDOC_CLASSES = NO
QUIET = YES
""",
}

primary_domain = "cpp"
highlight_language = "cpp"
