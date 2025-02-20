.. _installation:

Installation Guide
=====

It is recommended to create a virtual environment before installing RNA2seg to isolate dependencies:  

.. code-block:: console  

   $ conda create --name rna2seg-env python=3.10

Then, install RNA2seg and its dependencies:  

.. code-block:: console  

   (rna2seg-env) $ pip install instanseg-torch==0.0.5
   (rna2seg-env) $ pip install rna2seg

You might need to manually install the pytorch compatible with your system. ( https://pytorch.org/get-started )
RNA2seg is now installed and ready to use. 🚀
