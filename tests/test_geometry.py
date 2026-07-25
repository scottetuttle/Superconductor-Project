from shs.geometry import load_geometry, create_mesh


def test_load_geometry():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    assert geometry.name == "NbN_test_film"
    assert geometry.nx == 100
    assert geometry.width > 0



def test_mesh_generation():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    assert mesh.nx == 100
    assert mesh.ny == 100
    assert mesh.dx > 0
    assert mesh.dy > 0