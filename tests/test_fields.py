from shs.geometry import load_geometry, create_mesh
from shs.physics import Fields


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