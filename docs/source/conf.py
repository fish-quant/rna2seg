# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'RNA2seg'
copyright = '2025, Thomas Defard, Alice Blondel'
author = 'Thomas Defard, Alice Blondel'
release = '0.0.2'

html_logo = '../../img/mini_logo.svg'
html_static_path = ['_static']

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    "nbsphinx",
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]
templates_path = ['_templates']
exclude_patterns = []
autodoc_mock_imports = [
    "scanpy", "cv2", "albumentations", "pyarrow", "anndata", "pandas", "tqdm",
    "sklearn", "numpy", "networkx", "leidenalg", "scikit-image", "seaborn", "tifffile",
    "ssam", "scipy", "matplotlib", "torch", "torchvision", "dask", "time", "spatialdata",
    "pathlib", "rasterio", "cellpose", "geopandas", "xarray", "sopa", "shapely", "skimage",
    "json", "logging", "scipy.ndimage", "torch.utils.data", "torch.nn", "torch.nn.functional",
    "torch.optim", "cellpose.transforms", "albumentations.core.transforms_interface.ImageOnlyTransform", 
    "numpy as np", "pandas as pd", "spatialdata as sd", "scipy.ndimage as ndi", 
    "torch.utils.data.Dataset", "cellpose.transforms as tf_cp", "rna_seg", "instanseg",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'
html_static_path = ['_static']

# -- Path setup --------------------------------------------------------------
# Add the source directory to the sys.path for module imports

import sys
import os
import nbsphinx

source_dirs = [
    '../../src/', '.', '..', '../..', '../../..', '../../../'
]
for source_dir in source_dirs:
    sys.path.insert(0, os.path.abspath(source_dir))