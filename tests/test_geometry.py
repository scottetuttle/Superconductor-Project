from shs.geometry import load_geometry, create_mesh


def test_geometry_region():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    assert geometry.regions[0].name == "NbN film"
    assert geometry.regions[0].material == "NbN"