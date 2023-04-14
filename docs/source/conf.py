# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'bhpwave'
copyright = '2023, Zachary Nasipak'
author = 'Zachary Nasipak'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              "nbsphinx"]

templates_path = ['_templates']
exclude_patterns = []

autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__call__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_theme_options = {
    "repository_url": "https://github.com/znasipak/bhpwave",
    "use_repository_button": True,
    "navbar_end": ["navbar-icon-links"],
}
html_title = "bhpwave"
html_context = {
    "default_mode": "light"
}

