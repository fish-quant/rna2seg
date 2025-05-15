import shapely
import numpy as np
from tqdm import tqdm
import geopandas as gpd
import spatialdata as sd
from sopa.utils.utils import get_spatial_image
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from sopa.utils.utils import to_intrinsic

def compute_polygon_intersection(list_polygon_patch, list_polygon_annotation):
    tree = shapely.STRtree(list_polygon_patch)
    conflicts = tree.query(list_polygon_annotation, predicate="intersects")
    return conflicts


def compute_consistent_cell(
        sdata: sd.SpatialData,
        key_shape_cell_seg: str,
        key_shape_nuclei_seg: str,
        key_cell_consistent: str,
        key_nuclei_consistent: str,
        image_key: str,
        threshold_intersection_contain: float = 0.95,
        threshold_intersection_intersect: float = 0.05,
        accepted_nb_nuclei_per_cell: set = {1},
        max_cell_nb_intersecting_nuclei: int = 1,
        key_unconsistent_cell: str | None = None,
        shape_in_intrinsic: bool = True,

):
    """

    Compute consistent cell and nuclei from the initial cell and nuclei segmentation

    :param key_shape_nuclei_seg: key of the nuclei segmentation in sdata.shapes
    :type str
    :param key_shape_cell_seg: key of the cell segmentation in sdata.shapes
    :type str
    :param key_cell_consistent: key of the consistent cell in sdata.shapes
    :type str
    :param key_nuclei_consistent: key of the consistent nuclei in sdata.shapes
    :type str
    :param image_key: key of the image in sdata.images
    :type str
    :param threshold_intersection_contain: a cell is considered as containing a nuclei if
    cell.intersection(nucleus).area >= threshold_intersection * nucleus.area.
    It can be interpreted as the proportion of
    the nuclei need to be inside the cell for the (cell, nuclei) pair to be considered consistent
    :type float
    :param accepted_nb_nuclei_per_cell: set of accepted number of nuclei per cell. If None, only one nucleus per cell is
    accepted. if {0,1} then only cell with 0 or 1 nucleus are accepted
    :type set
    :param threshold_intersection_intersect: threshold of intersection apply as follow: a cell is considered as in the
    intersection with a nuclei if cell.intersection(nucleus).area >= threshold_intersection_intersect * nucleus.area
    :param shape_in_intrinsic: if True, the shape are in intrinsic coordinates, if False, not in intrinsic coordinates
     and are save in intrinsic coordinates
    :type float
    :return: None
    """

    if not shape_in_intrinsic:
        sdata[key_shape_cell_seg] = to_intrinsic(
            sdata, sdata[key_shape_cell_seg], sdata[image_key])
        sdata[key_shape_nuclei_seg] = to_intrinsic(
            sdata, sdata[key_shape_nuclei_seg], sdata[image_key])

    if threshold_intersection_intersect is None:
        threshold_intersection_intersect = 1 - threshold_intersection_contain
    if accepted_nb_nuclei_per_cell is None:
        accepted_nb_nuclei_per_cell = {1}
    list_cell = list(sdata[key_shape_cell_seg].geometry)  # can be optimise ssee read_patches_cells sopa function ?
    list_nuclei = list(sdata[key_shape_nuclei_seg].geometry)  # can be optimise ssee read_patches_cells sopa function ?

    (cell_consistent_with_nuclei, nuclei_consistent_with_cell, cell_consistent_without_nuclei,
     nuclei_consistent_not_in_cell, list_polygon_cell_unconsistent, dict_polygon_index) = _calculate_consistent_cells(
        list_polygon_cell=list_cell,
        list_polygon_nuclei=list_nuclei,
        threshold_intersection_contain=threshold_intersection_contain,
        accepted_nb_nuclei_per_cell=accepted_nb_nuclei_per_cell,
        threshold_intersection_intersect=threshold_intersection_intersect,
        max_cell_nb_intersecting_nuclei=max_cell_nb_intersecting_nuclei,
    )
    dict_polygon_index["index_all_polygon"] = list(range(len(list_cell)))
    dict_polygon_index["list_cell_polygon"] = list_cell

    if len(cell_consistent_with_nuclei) == 0:
        raise ValueError("No consistent cell found")

    # add the consistante cells and nuclei to the spatialdata

    def _save_consistent_segmentation(sdata, list_polygon, image_key, key_consistent):
        image = get_spatial_image(sdata, image_key)
        geo_df = gpd.GeoDataFrame({"geometry": list_polygon})
        geo_df.index = image_key + geo_df.index.astype(str)
        geo_df = ShapesModel.parse(
            geo_df, transformations=get_transformation(image, get_all=True).copy()
        )
        sdata.shapes[key_consistent] = geo_df

    key_cell_consistent_with_nuclei = f"{key_cell_consistent}_with_nuclei"
    key_nuclei_consistent_with_cell = f"{key_nuclei_consistent}_in_cell"
    key_cell_consistent_without_nuclei = f"{key_cell_consistent}_without_nuclei"
    key_nuclei_consistent_not_in_cell = f"{key_nuclei_consistent}_not_in_cell"

    _save_consistent_segmentation(
        sdata=sdata,
        list_polygon=nuclei_consistent_with_cell,
        image_key=image_key,
        key_consistent=key_cell_consistent_with_nuclei)
    sdata.write_element(key_cell_consistent_with_nuclei, overwrite=True)

    _save_consistent_segmentation(
        sdata=sdata,
        list_polygon=cell_consistent_with_nuclei,
        image_key=image_key,
        key_consistent=key_nuclei_consistent_with_cell)
    sdata.write_element(key_nuclei_consistent_with_cell, overwrite=True)

    _save_consistent_segmentation(
        sdata=sdata,
        list_polygon=cell_consistent_without_nuclei,
        image_key=image_key,
        key_consistent=key_cell_consistent_without_nuclei)
    sdata.write_element(key_cell_consistent_without_nuclei, overwrite=True)

    _save_consistent_segmentation(sdata=sdata,
                                  list_polygon=nuclei_consistent_not_in_cell,
                                  image_key=image_key,
                                  key_consistent=key_nuclei_consistent_not_in_cell)
    sdata.write_element(key_nuclei_consistent_not_in_cell, overwrite=True)

    if key_unconsistent_cell is not None:
        _save_consistent_segmentation(
            sdata=sdata,
            list_polygon=list_polygon_cell_unconsistent,
            image_key=image_key,
            key_consistent=key_unconsistent_cell
        )
        sdata.write_element(key_nuclei_consistent_not_in_cell, overwrite=True)

    return sdata, dict_polygon_index


def _calculate_consistent_cells(list_polygon_cell,
                               list_polygon_nuclei,
                               threshold_intersection_contain=0.95,
                               accepted_nb_nuclei_per_cell=None,
                               threshold_intersection_intersect=None,
                               max_cell_nb_intersecting_nuclei=1,
                               ):

    assert 0 not in accepted_nb_nuclei_per_cell, "0 is not an accepted number of nuclei per cell"
    if accepted_nb_nuclei_per_cell is None:
        accepted_nb_nuclei_per_cell = {1}
    if threshold_intersection_intersect is None:
        threshold_intersection_intersect = 1 - threshold_intersection_contain

    n_cells = len(list_polygon_cell)
    n_nuclei = len(list_polygon_nuclei)
    list_all_cell = list_polygon_cell + list_polygon_nuclei
    tree = shapely.STRtree(list_all_cell)
    conflicts = tree.query(list_all_cell, predicate="intersects")
    conflicts = conflicts[:, conflicts[0] != conflicts[1]].T

    cell_indice_to_keep = {}  # key cell , value set of nuclei indice inside the cell
    cell_indice_to_remove = {}  # key cell, value set of nuclei indice inside the cell
    nuclei_intersecting_cells = {}  # key nuclei, value set of cell indice intersecting the nuclei

    # initialize the cell_indice_to_keep
    for i in range(n_cells):
        cell_indice_to_keep[i] = set()
        cell_indice_to_remove[i] = set()

    for i in range(n_nuclei):
        nuclei_intersecting_cells[i] = set()

    for i, (i1, i2) in tqdm(enumerate(conflicts), desc="Resolving conflicts"):
        if i1 == i2:  # same polygon
            continue
        if i1 < n_cells and i2 < n_cells:  # both polygons are cells
            continue
        if i1 >= n_cells and i2 >= n_cells:  # both polygons are nuclei
            continue
        if i1 < i2:
            cell_index_all, nucleus_index_all = i1, i2
        else:
            nucleus_index_all, cell_index_all = i1, i2
        cell, nucleus = list_all_cell[cell_index_all], list_all_cell[nucleus_index_all]

        if not cell.is_valid:
            cell = cell.buffer(0)
            assert cell.is_valid
        if not nucleus.is_valid:
            nucleus = nucleus.buffer(0)
            assert nucleus.is_valid

        if cell.contains(nucleus):
            cell_indice_to_keep[cell_index_all].add(nucleus_index_all - n_cells)
        elif threshold_intersection_contain > 0:
            cell_nuc_inter = cell.intersection(nucleus).area / nucleus.area
            if cell_nuc_inter >= threshold_intersection_contain:
                cell_indice_to_keep[cell_index_all].add(nucleus_index_all - n_cells)
            else:
                if cell_nuc_inter > threshold_intersection_intersect:
                    cell_indice_to_remove[cell_index_all].add(nucleus_index_all - n_cells)
                    nuclei_intersecting_cells[nucleus_index_all - n_cells].add(cell_index_all)

    # compute consistent nuclei
    # nb_cell_per_nuclei = [len(cell_indice_to_keep[cell]) for cell in range(n_cells, len(list_all_cell))]

    # remove cell that have many nuclei inside
    cell_consistent_with_nuclei = [cell for cell in cell_indice_to_keep if
                                   len(cell_indice_to_keep[cell]) in accepted_nb_nuclei_per_cell]
    # remove cell that have "unconsistent" nuclei
    cell_set2remove = set([cell for cell in cell_indice_to_remove if len(cell_indice_to_remove[cell]) > 0])
    cell_consistent_with_nuclei = list(set(cell_consistent_with_nuclei) - cell_set2remove)

    # get cell without nuclei

    # keep cell that does not contain any nuclei
    cell_consistent_without_nuclei = [cell for cell in cell_indice_to_keep if len(cell_indice_to_keep[cell]) in [0]]

    # keep cell that have no nuclei intersecting
    cell_consistent_without_nuclei = list(set(cell_consistent_without_nuclei) - cell_set2remove)

    set_all_cell_consistent = set(cell_consistent_with_nuclei).union(set(cell_consistent_without_nuclei))
    assert len(set(cell_consistent_with_nuclei).intersection(set(cell_consistent_without_nuclei))) == 0, \
        "cell_consistent_without_nuclei and unique_cell_consistent have common elements, this should not happen"

    nuclei_consistent_with_cell = [list(cell_indice_to_keep[cell]) for cell in cell_consistent_with_nuclei]
    nuclei_consistent_with_cell = np.concatenate(nuclei_consistent_with_cell).astype(int)

    nuclei_consistent_not_in_cell = []
    for nucleus in range(n_nuclei):
        if len(nuclei_intersecting_cells[nucleus]) <= max_cell_nb_intersecting_nuclei:
            if nucleus not in nuclei_consistent_with_cell:
                nuclei_consistent_not_in_cell.append(nucleus)
                assert len(nuclei_intersecting_cells[nucleus].intersection(set_all_cell_consistent)) == 0

    nuclei_consistent_not_in_cell = list(set(nuclei_consistent_not_in_cell) - set(nuclei_consistent_with_cell))

    list_polygon_cell_consistent = [list_polygon_cell[i] for i in cell_consistent_with_nuclei]
    list_polygon_nuclei_consistent = [list_polygon_nuclei[i] for i in nuclei_consistent_with_cell]
    list_polygon_cell_consistent_without_nuclei = [list_polygon_cell[i] for i in cell_consistent_without_nuclei]
    list_polygon_nuclei_consistent_not_in_cell = [list_polygon_nuclei[i] for i in nuclei_consistent_not_in_cell]

    # compute unconcistent cell
    set_all_cell = set(list(range(n_cells)))
    set_all_cell_consistent = set(cell_consistent_with_nuclei).union(set(cell_consistent_without_nuclei))
    set_all_cell_unconsistent = set_all_cell - set_all_cell_consistent
    list_polygon_cell_unconsistent = [list_polygon_cell[i] for i in set_all_cell_unconsistent]

    dict_polygon_index = {
        "cell_consistent_with_nuclei": cell_consistent_with_nuclei,
        "nuclei_consistent_with_cell": nuclei_consistent_with_cell,
        "cell_consistent_without_nuclei": cell_consistent_without_nuclei,
        "nuclei_consistent_not_in_cell": nuclei_consistent_not_in_cell,
        "cell_unconsistent": list(set_all_cell_unconsistent),
        "cell_indice_to_keep": cell_indice_to_keep

    }
    return (list_polygon_cell_consistent, list_polygon_nuclei_consistent,
            list_polygon_cell_consistent_without_nuclei, list_polygon_nuclei_consistent_not_in_cell,
            list_polygon_cell_unconsistent,
            dict_polygon_index)
