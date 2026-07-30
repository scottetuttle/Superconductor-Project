import numpy as np

from shs.geometry.database import load_geometry
from shs.geometry.mesh import create_mesh
from shs.mapping import build_contact_map


def test_contact_map_exists():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    contact_map = build_contact_map(
        geometry,
        mesh,
    )

    assert len(contact_map.contact_masks) > 0


def test_contact_masks_are_boolean():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    contact_map = build_contact_map(
        geometry,
        mesh,
    )

    for mask in contact_map.contact_masks.values():

        assert mask.dtype == bool


def test_contact_mask_shape():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    contact_map = build_contact_map(
        geometry,
        mesh,
    )

    for mask in contact_map.contact_masks.values():

        assert mask.shape == (
            mesh.ny,
            mesh.nx,
        )


def test_contact_contains_cells():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    contact_map = build_contact_map(
        geometry,
        mesh,
    )

    for mask in contact_map.contact_masks.values():

        assert np.any(mask)