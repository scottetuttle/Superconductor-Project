from shs.geometry import load_geometry, create_mesh


def test_geometry_contacts():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    assert geometry.contacts[0].name == "left_current"
    assert geometry.contacts[0].contact_type == "current"