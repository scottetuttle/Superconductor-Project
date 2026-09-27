from shs.boundaries import load_boundary_set
from shs.boundaries import (
    BoundarySide,
    BoundaryType,
)
import json


def test_load_boundary_json():

    boundaries = load_boundary_set(
        "configs/simulations/nbn_hotspot_test.json"
    )

    left = boundaries.get(
        BoundarySide.LEFT
    )

    assert left.type == BoundaryType.FIXED_TEMPERATURE

    with open("configs/simulations/nbn_hotspot_test.json", encoding="utf-8") as file:
        expected = json.load(file)["boundaries"]["left"]["temperature"]
    assert left.temperature == expected


def test_insulating_boundary():

    boundaries = load_boundary_set(
        "configs/simulations/nbn_hotspot_test.json"
    )

    top = boundaries.get(
        BoundarySide.TOP
    )

    assert top.type == BoundaryType.INSULATING
