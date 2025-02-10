import sys
import tifffile
import rasterio
import numpy as np
from tqdm import tqdm 
import spatialdata as sd
from pathlib import Path
from rasterio.features import rasterize
from cellpose.dynamics import  labels_to_flows

import logging
log = logging.getLogger(__name__)

from rna2seg.dataset_zarr.utils.utils_preprocessing import labels_to_flows_omnipose
from rna2seg._constant import RNAsegFiles

def get_segmentation_crop(cell_segmentation, bounds, shape):
    try :
        polygons = list(cell_segmentation.geometry)
        print(len(polygons))
        x_trans = bounds[0]
        y_trans = bounds[1]
        transform = rasterio.Affine(a=1, b=0, c=x_trans, d=0, e=1, f=y_trans)
        cell_segmentation = rasterize(((polygons[i], i+1) for i in range(len(polygons))),
                                        out_shape=shape, transform=transform,
                                        fill=0, all_touched=True, dtype=np.uint16)
    except ValueError as e:
        assert len(polygons) == 0, "no polygon in the agreement_segmentation"
        return np.zeros(shape)
    return cell_segmentation

def precompute_flow(sdata,
                    image_key: str, 
                    shape_patch_key: str,
                    segmentation_key: str,
                    channel_dapi: int,
                    patch_dir_csv = str,
                    list_path_index: list[int] | None = None,
                    compute_cellpose: bool = True,
                    compute_omnipose: bool = False,
                    key_cell : str = None,
                    shape  = (1200, 1200)):
    
    assert compute_cellpose + compute_omnipose >= 1, "at least one of the two should be computed"
    if list_path_index is None:
        list_path_index = list(range(len(sdata[shape_patch_key].geometry)))

    for patch_index in tqdm(list_path_index, desc="precompute flow", file = sys.stdout):
        if patch_index >= len(sdata[shape_patch_key].geometry):
            log.warning(f"patch_index {patch_index} is out of bound")
            continue
        patch = sdata[shape_patch_key].geometry[patch_index]
        bounds = [int(x) for x in patch.bounds]

        if shape is None:
            
            image = sdata[image_key]["scale0"].image.sel(
                c=channel_dapi,
                x=slice(bounds[0], bounds[2]),
                y=slice(bounds[1], bounds[3]),
            ).values
            shape = image.shape

        cell_segmentation = sdata[segmentation_key].cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]
        full_segmentation = get_segmentation_crop(cell_segmentation, bounds = bounds, shape = shape,)

        path_save = Path(patch_dir_csv) / f'{patch_index}/{key_cell}'
        print(path_save)
        path_save.mkdir(exist_ok=True, parents=True)
        if shape is None:
            tifffile.imwrite(path_save / RNAsegFiles.IMAGE, image)

        if compute_cellpose:
            flow = labels_to_flows([full_segmentation], files=None,
                                        device=None, redo_flows=False)
            np.save(path_save / RNAsegFiles.LABEL_AND_MASK_FLOW, flow[0])

        if compute_omnipose:
            flow = labels_to_flows_omnipose(labels=[full_segmentation],
            links=None,
            files=None,
            use_gpu=False,
            device=None,
            omni=True,
            redo_flows=True,
            dim=2)
            np.save(path_save / RNAsegFiles.LABEL_AND_MASK_FLOW_OMNIPOSE, flow[0])
        ## save the nb of cell in a file
        nb_cell = len(np.unique(flow[0][0])) - 1
        with open(path_save / RNAsegFiles.NB_CELL_FILE, "w") as f:
            f.write(str(nb_cell))


if __name__ == "__main__":
    channel_dapi = ["DAPI"]
    list_path_index = None
    segmentation_key = "manual_annotation_for_finetuning"
    key_cell = "manual_annotation_for_finetuning"
    shape_patch_key = "sopa_patches_rna2seg_1200_50"
    image_key = "data_release_baysor_merfish_gut_RNAseg_DATASET_z3"
    patch_dir_csv = "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/mouse_ileum/RNAseg_DATASET.zarr/.sopa_cache/rna_seg_1200_50"
    sdata = sd.read_zarr("/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/mouse_ileum/RNAseg_DATASET.zarr")

    precompute_flow(
        sdata, image_key, shape_patch_key, segmentation_key, channel_dapi, 
        patch_dir_csv,  key_cell=key_cell, list_path_index=list_path_index, shape  = (1200, 1200),
    )