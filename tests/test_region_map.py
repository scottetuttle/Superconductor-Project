import numpy as np

from shs.geometry.database import load_geometry
from shs.geometry.mesh import create_mesh

from shs.mapping.region_map import build_region_map


def test_region_map_shape():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    region_map = build_region_map(
        geometry,
        mesh
    )

    assert region_map.region_ids.shape == (
        mesh.ny,
        mesh.nx
    )


def test_region_map_all_film():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    region_map = build_region_map(
        geometry,
        mesh
    )

    assert np.all(
        region_map.region_ids == 0
    )


def test_region_dictionary():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    region_map = build_region_map(
        geometry,
        mesh
    )

    assert 0 in region_map.region_names

    assert (
        region_map.region_names[0]
        == "film"
    )