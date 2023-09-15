# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "LuPNT"
copyright = "2023, Stanford NAV Lab"
author = "Keidai Iiyama, Guillem Casadesus Villa, Marta Cortinovis, Ashwin Kanhere"
release = '"0.1"'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx_sitemap",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
    "breathe",
    "exhale",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "restructuredtext",
    ".md": "markdown",
}

# Setup the breathe extension
breathe_projects = {"My Project": "./_build/xml"}
breathe_default_project = "My Project"

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "doxygenStripFromPath": "..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle": "Library API",
    # Suggested optional arguments
    "createTreeView": True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": """
         INPUT = ../../lupnt/
         EXCLUDE = "../../lupnt/3rdparty/"
     """,
}

# Tell sphinx what the primary language being documented is.
primary_domain = "cpp"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "cpp"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_logo = "img/nav_lab_logo.png"

html_theme_options = {
    "canonical_url": "",
    "analytics_id": "",
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    "style_nav_header_background": "#8C1515",
    "logo_only": False,
    "display_version": True,
    # Toc options
    "collapse_navigation": False,
    # "sticky_navigation" : False,
    # "navigation_depth" : 4,
    # "includehidden" : True,
    # "titles_only" : True,
}

# Default to the main branch, default to main and tags not existing
linkcode_revision = "development"
in_main = False
tagged = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []
