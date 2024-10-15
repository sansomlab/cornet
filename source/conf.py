# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import sphinx_rtd_theme
from pathlib import Path
from ruffus import *
import shutil
import glob
import sqlite3
import numpy as np
import pandas as pd
import textwrap
import subprocess
from scipy.stats.mstats import gmean
import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools

# import local pipeline utility functions
#from ../pipelines/pipeline_utils import templates

sys.path.insert(0, os.path.abspath('..'))
#sys.path.insert(0, os.path.abspath('../pipelines')
sys.path.insert(1, os.path.abspath('../pipelines'))
print(sys.path)

#from pipeline_utils import templates

print(sys.executable)
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CORNET'
copyright = '2024, Tarrion Baird, Stephen Sansom'
author = 'Tarrion Baird, Stephen Sansom'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration



# Autodoc defaults
autodoc_default_options = {
#        'members': 'var1, var2',
        'member-order': 'bysource'
#        'special-members': '__init__',
#        'undoc-members': True,
#        'exclude-members': '__weakref__'
    }

extensions = ['sphinx.ext.autosectionlabel','sphinx.ext.autodoc','sphinx.ext.coverage', 'sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
