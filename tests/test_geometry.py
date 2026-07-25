from shs.geometry import load_geometry


def test_load_geometry():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    assert geometry.name == "NbN_test_film"
    assert geometry.width == 5e-6
    assert geometry.nx == 100