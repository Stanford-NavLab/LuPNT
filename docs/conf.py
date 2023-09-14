# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'LuPNT'
copyright = '2023, Stanford NAV Lab'
author = 'Keidai Iiyama, Guillem Casadesus Villa, Marta Cortinovis, Ashwin Kanhere'
release = '"0.1"'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# extensions = ["myst_parser"]
extensions = ['sphinx_rtd_theme',
            'sphinx.ext.autodoc',
            'sphinx.ext.napoleon']
            # 'sphinx.ext.linkcode']

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'restructuredtext',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
# html_static_path = []

html_logo = "img/nav_lab_logo.png"

html_theme_options = {
    "style_nav_header_background" : "#8C1515",
    # "display_version" : True,
    "collapse_navigation" : False,
    # "sticky_navigation" : False,
    # "navigation_depth" : 4,
    # "includehidden" : True,
    # "titles_only" : True,
    "logo_only" : False,
    "display_version" : True,
}

#Default to the main branch, default to main and tags not existing
linkcode_revision = "main"
in_main = False
tagged = False
