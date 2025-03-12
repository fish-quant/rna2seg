


import shutil
import sys

sys.path =['/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/pydev',
           '/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/third_party/thriftpy',
           '/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/pydev',
           '/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/pycharm_display', '/home/tom/anaconda3/envs/rna2seg-env/lib/python310.zip', '/home/tom/anaconda3/envs/rna2seg-env/lib/python3.10', '/home/tom/anaconda3/envs/rna2seg-env/lib/python3.10/lib-dynload', '/home/tom/anaconda3/envs/rna2seg-env/lib/python3.10/site-packages', '/home/tom/Bureau/phd/rna_seg2paper/rna_seg_pkg/src', '/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/pycharm_matplotlib_backend', '/home/tom/anaconda3/envs/rna2seg-env/lib/python3.10/site-packages/setuptools/_vendor', '/home/tom/Bureau/phd/rna_seg2paper', '/home/tom/Bureau/phd/rna_seg2paper/rna_seg_pkg/src']


from pathlib import Path


#from _constant_path import PathTest
import pytest
import spatialdata as sd
import torch
from tqdm import tqdm

from rna2seg.dataset_zarr import (RNA2segDataset, compute_consistent_cell,
                                  custom_collate_fn)

# ADD YOUR PATH
MERFISH_ZARR_PATH = Path("/media/tom/Transcend/open_merfish/test_spatial_data/test005/sub_mouse_ileum.zarr")


class VariableTest:
    image_key = "staining_z3"
    points_key = "transcripts"
    patch_width = 1200
    patch_overlap = 50
    min_transcripts_per_patch = 0
    merfish_zarr_path = MERFISH_ZARR_PATH
    folder_patch_rna2seg = Path(merfish_zarr_path) / f".rna2seg_{patch_width}_{patch_overlap}"

    channels_dapi = ["DAPI"]
    channels_cellbound = ["Cellbound1"]
    gene_column_name = "gene"

    key_cell_segmentation = "Cellbound1"
    key_nuclei_segmentation = "DAPI"
    # to name for future shape that will be created in the sdata
    key_cell_consistent = "Cellbound1_consistent"
    key_nucleus_consistent = "DAPI_consistent"


def clean():
    patch_width = VariableTest.patch_width
    patch_overlap = VariableTest.patch_overlap
    key_shape = f"sopa_patches_rna2seg_{patch_width}_{patch_overlap}"
    list_shape = [VariableTest.key_cell_consistent + "_with_nuclei",
                  VariableTest.key_nucleus_consistent + "_in_cell",
                  VariableTest.key_nucleus_consistent + "_not_in_cell",
                  VariableTest.key_cell_consistent + "_without_nuclei",
                  ]

    list_folder_to_remove = [VariableTest.merfish_zarr_path / f".{key_shape}",
                             VariableTest.folder_patch_rna2seg / f".sopa_cache",
                             VariableTest.merfish_zarr_path / f"shapes/{list_shape[0]}",
                             VariableTest.merfish_zarr_path / f"shapes/{list_shape[1]}",
                             VariableTest.merfish_zarr_path / f"shapes/{list_shape[2]}",
                             VariableTest.merfish_zarr_path / f"shapes/{list_shape[3]}",
                             VariableTest.merfish_zarr_path / f"shapes/sopa_patches_rna2seg_1200_150",
                             ]
    for folder in list_folder_to_remove:
        if folder.exists():
            print(f"remove {folder}")
            shutil.rmtree(folder)

@pytest.mark.run(order=1)
def test_clean_before():
    clean()


@pytest.mark.run(order=2)
def test_consistent():
    sdata = sd.read_zarr(VariableTest.merfish_zarr_path)

    sdata, _ = compute_consistent_cell(
        sdata=sdata,
        key_shape_nuclei_seg=VariableTest.key_nuclei_segmentation,
        key_shape_cell_seg=VariableTest.key_cell_segmentation,
        key_cell_consistent=VariableTest.key_cell_consistent,
        key_nuclei_consistent=VariableTest.key_nucleus_consistent,
        image_key=VariableTest.image_key,
        threshold_intersection_contain=0.95,
        threshold_intersection_intersect=0.05,
        accepted_nb_nuclei_per_cell=None,
        max_cell_nb_intersecting_nuclei=1,
    )

    assert len(sdata['Cellbound1_consistent_with_nuclei']) == len(sdata['DAPI_consistent_in_cell']) == 204
    assert len(sdata['DAPI_consistent_not_in_cell']) == 171
    assert len(sdata['Cellbound1_consistent_without_nuclei']) == 248


@pytest.mark.run(order=3)
def test_after_before():
    clean()