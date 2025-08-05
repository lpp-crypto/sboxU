# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# import os
# import sys
# sys.path.insert(0, os.path.abspath('../'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'sboxU'
copyright = '2025, Léo Perrin'
author = 'Léo Perrin'
release = '2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",  # needed to generate doc from docstrings
    "sphinx.ext.napoleon", # needed to use google-style docstrings
    "sphinx.ext.viewcode",
] 

# source: https://gist.github.com/blepfo/6c37b6c39efcdc74504d33418607bba0
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False


templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
