import numpy as np
from tqdm import tqdm
from scipy import ndimage as ndi
from matplotlib import colors as mcolors
from sopa._sdata import to_intrinsic

from rna2seg.dataset_zarr.RNA2segDataset import rna2img


### RNA to Image

def apply_transformation(row, transformation_matrix=None):
    x, y = row["x"], row["y"]
    if transformation_matrix is None:
        return x, y
    vector = np.array([x, y, 1])
    transformed = transformation_matrix @ vector
    return transformed[0], transformed[1]

def get_patch_df(bounds, df, transformation_matrix=None):

    # Get bounds in df coordinates
    if transformation_matrix is None:
        new_bounds = bounds
    else:
        x_min, y_min, x_max, y_max = bounds
        inverse_matrix = np.linalg.inv(transformation_matrix)
        corners = np.array([[x_min, y_min, 1],[x_max, y_min, 1], [x_min, y_max, 1],[x_max, y_max, 1]])
        transformed_corners = np.dot(inverse_matrix, corners.T).T
        x_coords = transformed_corners[:, 0]
        y_coords = transformed_corners[:, 1]
        new_bounds = x_coords.min(), y_coords.min(), x_coords.max(), y_coords.max()

    # Extract patch points
    filtered_df = df[
        (df["x"] >= new_bounds[0]) & (df["x"] <= new_bounds[2]) &
        (df["y"] >= new_bounds[1]) & (df["y"] <= new_bounds[3])
    ]

    # Transform coordinates
    filtered_df[["x", "y"]] = filtered_df.map_partitions(
        lambda df: df.apply(apply_transformation, axis=1, result_type="expand")
    )
    return filtered_df

def get_rna_img(dataset, bounds, key_transcripts="transcripts", transformation_matrix=None, dict_gene_value=None):
    print("Get RNA image ...")
    st_segmentation = dataset.st_segmentation
    df = st_segmentation.sdata[key_transcripts]
    img = st_segmentation.sdata[st_segmentation.image_key]
    df = to_intrinsic(st_segmentation.sdata, df, img)
    patch_df = get_patch_df(bounds, df, transformation_matrix)
    patch_df = patch_df.compute()
    if dict_gene_value is None:
        dict_gene_value = {gene: 1 for gene in patch_df["gene"].unique()}

    img = rna2img(
        patch_df, 
        dict_gene_value=dict_gene_value,  
        image_shape=(bounds[3]-bounds[1], bounds[2]-bounds[0], 3),
        offset_x=bounds[0], offset_y=bounds[1], gene_column="gene", 
        max_filter_size=5,
        addition_mode=True
    )
    return img

### Stainings

def get_staining_img(dataset, bounds):
    print("Get image ...")
    st_segmentation = dataset.st_segmentation
    xmin, ymin, xmax, ymax = bounds
    stainings = [*st_segmentation.channels_dapi, *st_segmentation.channels_cellbound]
    image = st_segmentation.image.sel(
        c=stainings,
        x=slice(xmin, xmax),
        y=slice(ymin, ymax),
    ).values
    return image


### Segmentation

def create_cell_contours(seg_mask, size_line = 5, min_size= 2500): 
    
    contour_mask = np.zeros(seg_mask.shape)
    cell_to_contour = np.unique(seg_mask)

    nb_cell_in = 0
    for cell in tqdm(cell_to_contour):
        mask = (seg_mask == cell).astype(int)
        if cell == 0 or mask.sum()<min_size:
            continue
        nb_cell_in += 1
        contour_mask += (ndi.maximum_filter(mask, size=size_line) - ndi.minimum_filter(mask, size=size_line))
    
    return contour_mask

def get_segmentation_img(dataset, image, bounds, color="red", size_line=5, key_cell="cell"):
    print("Get segmentation image ...")
    st_segmentation = dataset.st_segmentation
    image_with_mask=image/image.max()
    if image_with_mask.ndim == 2:
        image_with_mask = np.stack([image_with_mask] * 3, axis=-1)

    segmentation =  st_segmentation.get_segmentation_crop(
        bounds = bounds, shape = image_with_mask.shape[0:2], key_cell= key_cell)
    cell_contour_mask = create_cell_contours(segmentation, size_line = size_line)
    
    image_with_mask[cell_contour_mask > 0] = mcolors.hex2color(color)

    return image_with_mask