import numpy as np

from shs.geometry.database import load_geometry
from shs.geometry.mesh import create_mesh

from shs.mapping import (
    build_region_map,
    build_material_map,
)

from shs.materials.database import get_material


def test_material_map_shape():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    region_map = build_region_map(
        geometry,
        mesh
    )

    material = get_material("NbN")

    material_map = build_material_map(
        region_map,
        material
    )

    assert material_map.material_ids.shape == (
        mesh.ny,
        mesh.nx
    )


def test_material_dictionary():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    region_map = build_region_map(
        geometry,
        mesh
    )

    material = get_material("NbN")

    material_map = build_material_map(
        region_map,
        material
    )

    assert 0 in material_map.materials

    assert (
        material_map.materials[0].name
        == "NbN"
    )


def test_property_arrays():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    region_map = build_region_map(
        geometry,
        mesh
    )

    material = get_material("NbN")

    material_map = build_material_map(
        region_map,
        material
    )

    assert np.all(
        material_map.thermal_conductivity
        == material.thermal_conductivity
    )

    assert np.all(
        material_map.heat_capacity
        == material.heat_capacity
    )

    assert np.all(
        material_map.normal_resistivity
        == material.normal_resistivity
    )

    assert np.all(
        material_map.Tc
        == material.Tc
    )


def test_all_cells_are_material_zero():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    region_map = build_region_map(
        geometry,
        mesh
    )

    material = get_material("NbN")

    material_map = build_material_map(
        region_map,
        material
    )

    assert np.all(
        material_map.material_ids == 0
    )