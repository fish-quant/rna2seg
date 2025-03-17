# RNA2seg - A generalist model for cell segmentation in image-based spatial transcriptomics

<p align="center">
    <img src="./img/logo.png" width="200"/>
</p>

## Overview

**RNA2seg** is a deep learning-based segmentation model designed to improve cell segmentation in **Imaging-based Spatial Transcriptomics (IST)**. Traditional IST methods rely on nuclear and membrane staining to define cell boundaries, but segmentation can be challenging due to the variable quality of membrane markers.  

RNA2seg addresses this issue by integrating an **arbitrary number of staining channels** along with **RNA spatial distributions** to enhance segmentation accuracy, particularly in regions with low-quality membrane staining. It is built on **SpatialData**, enabling seamless processing and analysis of spatial transcriptomics data.  

![overview](./img/overview.png)

### **Key Features:**  
- **Multi-channel input**: Leverages nuclear, membrane, and RNA spatial data.  
- **Robust segmentation**: Outperforms state-of-the-art methods across multiple datasets.  
- **Zero-shot learning**: Generalizes to unseen datasets without requiring manual annotations.  
- **Built on SpatialData**: Ensures compatibility with modern spatial transcriptomics pipelines.  

RNA2seg is designed for researchers working with spatial transcriptomics data who need **accurate and adaptable segmentation** across different experimental conditions.


## Documentation

Check RNA2seg's [documentation](https://rna2seg.readthedocs.io/en/latest/) to get started. It contains installation explanations and tutorials.

## new released 

14/03/25 RNA2seg 0.1:
 - fix bug patch for parameter min_points_per_patch
 - remove dimater parameter in RNA2geg
 - re-order dependency

12/03/25 RNA2seg 0.0.7 :
 - add pretrained model for brain data
 - fix RNA embbeding bug


## Installation

It is recommended to create a virtual environment before installing RNA2seg to isolate dependencies:  

```
$ conda create --name rna2seg-env python=3.10
```
Then, install RNA2seg and its dependencies:  

```
(rna2seg-env) $ pip install rna2seg
```

RNA2seg is now installed and ready to use. 

## Support

If you have any bug or questions relative to the package, please open an issue 

## Citation

If you use this library, please cite:

```
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
```

## Datasets

All datasets, required to reproduce the results of publication "RN2Aseg: a generalist model for cell segmentation1
in image-based spatial transcriptomics2" are available at ...


## test


Toy dataset to run test and notebook can be downloaded at: https://cloud.minesparis.psl.eu/index.php/s/qw2HaDVxwwy1EOK
(data from Petukhov. et al. Nat Biotechnol 40, 345–354 (2022). https://doi.org/10.1038/s41587-021-01044-w )
