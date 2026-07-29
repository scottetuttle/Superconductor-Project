from shs.boundaries import load_boundary_set
from shs.boundaries import (
    BoundarySide,
    BoundaryType,
)


def test_load_boundary_json():

    boundaries = load_boundary_set(
        "configs/simulations/nbn_hotspot_test.json"
    )

    left = boundaries.get(
        BoundarySide.LEFT
    )

    assert left.type == BoundaryType.FIXED_TEMPERATURE

    assert left.temperature == 4.2


def test_insulating_boundary():

    boundaries = load_boundary_set(
        "configs/simulations/nbn_hotspot_test.json"
    )

    top = boundaries.get(
        BoundarySide.TOP
    )

    assert top.type == BoundaryType.INSULATING