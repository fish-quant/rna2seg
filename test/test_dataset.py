import sys


sys.path = ['/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/pydev',
            '/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/third_party/thriftpy',
            '/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/pydev',
            '/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/pycharm_display',
            '/home/tom/anaconda3/envs/test_rna2seg4/lib/python310.zip',
            '/home/tom/anaconda3/envs/test_rna2seg4/lib/python3.10',
            '/home/tom/anaconda3/envs/test_rna2seg4/lib/python3.10/lib-dynload',
            '',
            '/home/tom/anaconda3/envs/test_rna2seg4/lib/python3.10/site-packages',
            '/home/tom/Bureau/phd/rna_seg2paper/rna_seg_pkg/src',
            '/home/tom/Bureau/spatialdata-main/src',
            '/home/tom/.local/share/JetBrains/Toolbox/apps/pycharm-professional/plugins/python/helpers/pycharm_matplotlib_backend',
            '/home/tom/anaconda3/envs/test_rna2seg4/lib/python3.10/site-packages/setuptools/_vendor',
            '/home/tom/Bureau/phd/rna_seg2paper',
            '/home/tom/Bureau/phd/rna_seg2paper/rna_seg_pkg/src']

import warnings

# from _constant_path import PathTest
import pytest

warnings.filterwarnings("ignore")

import warnings
from pathlib import Path

import cv2
import spatialdata as sd

### create patch rna2seg



warnings.filterwarnings("ignore")
import json
import shutil
from pathlib import Path

import cv2
import pandas as pd
import spatialdata as sd

from rna2seg.dataset_zarr.patches import create_patch_rna2seg





class VariableTest:
    image_key = "staining_z3"
    points_key = "transcripts"
    patch_width = 1200
    patch_overlap = 150
    min_transcripts_per_patch = 0
    merfish_zarr_path = Path("/media/tom/Transcend/open_merfish/test_spatial_data/test005/sub_mouse_ileum.zarr")
    folder_patch_rna2seg = Path(merfish_zarr_path) / f".rna2seg_{patch_width}_{patch_overlap}"

    channels_dapi = ["DAPI"]
    channels_cellbound = ["Cellbound1"]
    gene_column_name = "gene"


def clean():
    patch_width = VariableTest.patch_width
    patch_overlap = VariableTest.patch_overlap
    key_shape = f"rna2seg_{patch_width}_{patch_overlap}"

    list_folder_to_remove = [VariableTest.merfish_zarr_path / f".{key_shape}",
                             VariableTest.merfish_zarr_path / f".sopa_cache",
                             VariableTest.merfish_zarr_path / f"shapes/sopa_patches_rna2seg_{patch_width}_{patch_overlap}"]

    for folder in list_folder_to_remove:
        if folder.exists():
            print(f"remove {folder}")
            shutil.rmtree(folder)


@pytest.mark.run(order=1)
def test_clean_before_test():
    clean()


@pytest.mark.run(order=2)
def test_create_patch_rna2seg():
    # MODIFY WITH YOUR PATH

    # load sdata and set path parameters
    sdata = sd.read_zarr(VariableTest.merfish_zarr_path)

    #create patch in the sdata and precompute transcipt.csv for each patch with sopa
    create_patch_rna2seg(sdata=sdata,
                         image_key=VariableTest.image_key,
                         points_key=VariableTest.points_key,
                         patch_width=VariableTest.patch_width,
                         patch_overlap=VariableTest.patch_overlap,
                         min_transcripts_per_patch=0,
                         folder_patch_rna2seg=VariableTest.folder_patch_rna2seg,
                         overwrite=True)
    print(sdata)

    path_cache = sdata.path / ".rna2seg_1200_150/5"

    with open(path_cache / "bounds.json") as f:
        bounds = json.load(f)
    print(bounds)
    assert bounds['bounds'] == [1050.0, 1050.0, 2250.0, 2250.0]

    df = pd.read_csv(path_cache / "transcripts.csv")
    assert df.x[5] == 1457.131390465281

    key_shape = f"sopa_patches_rna2seg_{VariableTest.patch_width}_{VariableTest.patch_overlap}"
    assert key_shape in sdata.shapes.keys()


# check rna2seg dataset
@pytest.mark.run(order=3)
def test_create_patch_rna2seg():


    import albumentations as A

    from rna2seg.dataset_zarr import RNA2segDataset

    transform_resize = A.Compose([
        A.Resize(width=512, height=512, interpolation=cv2.INTER_NEAREST),
    ])

    sdata = sd.read_zarr(VariableTest.merfish_zarr_path)

    dataset = RNA2segDataset(
        sdata=sdata,
        channels_dapi=VariableTest.channels_dapi,
        channels_cellbound=VariableTest.channels_cellbound,
        patch_width=VariableTest.patch_width,
        patch_overlap=VariableTest.patch_overlap,
        gene_column=VariableTest.gene_column_name,
        transform_resize=transform_resize,
        patch_dir=VariableTest.folder_patch_rna2seg
    )

    assert len(dataset) == 9
    assert float(dataset[6]['rna_img'][0, 244, 244]) == 2.8158267014077865e-05
    assert len(dataset[6]) == 8

    from rna2seg.models import RNA2seg
    device = "cpu"

    rna2seg = RNA2seg(
        device,
        net='unet',
        flow_threshold=0.9,
        cellbound_flow_threshold=0.4,
        pretrained_model="default_pretrained"
    )

    input_dict = dataset[6]
    flow, cellprob, masks_pred, cells = rna2seg.run(
        path_temp_save=VariableTest.folder_patch_rna2seg,
        input_dict=input_dict
    )

    assert len(cells) == 87


# add test for save_shapes2zarr(dataset, segmentation_shape_name)


@pytest.mark.run(order=4)
def test_clean_after_test():
    clean()

###
