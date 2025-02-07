







from spatialdata._io import write_shapes
from spatialdata import SpatialData

# fonction from sopa 1.0.14
def save_shapes(
        sdata: SpatialData,
        name: str,
        overwrite: bool = False,
) -> None:
    if not sdata.is_backed():
        return

    elem_group = sdata._init_add_element(name=name, element_type="shapes", overwrite=overwrite)

    write_shapes(
        shapes=sdata.shapes[name],
        group=elem_group,
        name=name,
    )





