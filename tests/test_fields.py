from shs.geometry import load_geometry, create_mesh
from shs.physics import Fields

import numpy as np


def test_create_fields():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    fields = Fields.create(
        mesh,
        initial_temperature=3.0
    )


    assert fields.temperature.shape == (
        100,
        100
    )

    assert fields.voltage.shape == (
        100,
        100
    )

    assert fields.temperature[0,0] == 3.0


    # Superconducting order parameter

    assert fields.psi.shape == (
        mesh.ny,
        mesh.nx
    )

    assert np.iscomplexobj(
        fields.psi
    )


    # Initial superconducting state:
    #
    # psi = 1 + 0i
    #
    # amplitude = 1
    # phase = 0

    assert np.allclose(
        np.abs(fields.psi),
        1.0
    )

    assert np.allclose(
        np.angle(fields.psi),
        0.0
    )



def test_electromagnetic_fields():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    fields = Fields.create(
        mesh,
        initial_temperature=3.0
    )


    assert fields.electric_field_x.shape == (
        mesh.ny,
        mesh.nx
    )

    assert fields.current_density_x.shape == (
        mesh.ny,
        mesh.nx
    )

    assert fields.magnetic_field_x.shape == (
        mesh.ny,
        mesh.nx
    )