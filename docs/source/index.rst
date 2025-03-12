.. RNA2seg documentation master file, created by
   sphinx-quickstart on Tue Feb 11 16:27:05 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RNA2seg's documentation!
===================================

**RNA2seg** is a deep learning-based segmentation model designed to improve cell segmentation in **Imaging-based Spatial Transcriptomics (IST)**. Traditional IST methods rely on nuclear and membrane staining to define cell boundaries, but segmentation can be challenging due to the variable quality of membrane markers.  

RNA2seg addresses this issue by integrating an **arbitrary number of staining channels** along with **RNA spatial distributions** to enhance segmentation accuracy, particularly in regions with low-quality membrane staining. It is built on **SpatialData**, enabling seamless processing and analysis of spatial transcriptomics data.  

.. image:: ../../img/overview.png


Check out the :doc:`install` section for further information about how to install the package.

realesed version and code
------------
the package code is at  https://github.com/fish-quant/rna2seg

12/03/25 RNA2seg 0.0.7 :
 - fix RNA embbeding bug
 - add pretrained model for brain data

Contents
------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. toctree::

   install
   userguide
   pretrained_model
   module_rna2seg
   
Support
------------

If you have any questions relative to the package, please open an issue on  `GitHub <https://github.com/fish-quant/rna2seg>`_.

Citation
------------

If you use this library, please be sure to cite:

.. code-block:: text

   @article {Defard2025.03.03.641259,
	author = {Defard, Thomas and Blondel, Alice and Coleon, Anthony and Dias de Melo, Guilherme and Walter, Thomas and Mueller, Florian},
	title = {RNA2seg: a generalist model for cell segmentation in image-based spatial transcriptomics},
	elocation-id = {2025.03.03.641259},
	year = {2025},
	doi = {10.1101/2025.03.03.641259},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/03/11/2025.03.03.641259},
	eprint = {https://www.biorxiv.org/content/early/2025/03/11/2025.03.03.641259.full.pdf},
	journal = {bioRxiv}
}

Indices and tables
===================================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

