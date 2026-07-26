from shs.geometry import load_geometry, create_mesh


def test_geometry_contacts():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    assert geometry.contacts[0].name == "left_current"
    assert geometry.contacts[0].contact_type == "current"

def test_geometry_container():
    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )
    assert geometry.film is not None
    assert geometry.regions is not None
    assert geometry.contacts is not None
    assert geometry.film.width is not None
    assert geometry.film.height is not None