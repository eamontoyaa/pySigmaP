#!/usr/bin/env python
#
# pysigmap documentation build configuration file, created by
# sphinx-quickstart on Fri Jun  9 13:47:02 2017.
#
# If extensions (or modules to document with autodoc) are in another
# directory, add these directories to sys.path here. If the directory is
# relative to the documentation root, use os.path.abspath to make it
# absolute, like shown here.
#
# import pysigmap
import os
import sys

sys.path.insert(0, os.path.abspath("../"))
# sys.path.insert(0, os.path.abspath('../../'))
# sys.path.insert(0, os.path.abspath('../../pysigmap/'))
sys.path.insert(0, os.path.abspath("../pysigmap/"))


# -- General configuration ---------------------------------------------

# For not sorting alphabetically
autodoc_member_order = "bysource"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.napoleon",
    "matplotlib.sphinxext.plot_directive",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "pySigmaP"
copyright = "2020, E. A. Montoya-Araque, A. J. Aparicio-Ortube, D. G. Zapata-Medina, L. G. Arboleda-Monsalve and Universidad Nacional de Colombia"
author = "E. A. Montoya-Araque\\\A. J. Aparicio-Ortube\\\D. G. Zapata-Medina\\\L. G. Arboleda-Monsalve"
release = "0.1.9"

# -- Options for HTML output -------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Theme options are theme-specific and customize the look and feel of a
# theme further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_build/_static"]


# -- Options for LaTeX output ------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass
# [howto, manual, or own class]).
latex_documents = [
    (
        master_doc,
        "pysigmap.tex",
        "pySigmaP Documentation",
        "E. A. Montoya-Araque\\\A. J. Aparicio-Ortube\\\D. G. Zapata-Medina\\\L. G. Arboleda-Monsalve",
        "manual",
    ),
]
