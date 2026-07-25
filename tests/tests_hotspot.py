from shs.geometry import load_geometry, create_mesh
from shs.optics import GaussianHotspot


def test_gaussian_hotspot():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)


    hotspot = GaussianHotspot(
        x0=2.5e-6,
        y0=2.5e-6,
        amplitude=1e8,
        sigma=0.5e-6
    )


    Q = hotspot.generate(mesh)


    assert Q.shape == (
        100,
        100
    )

    assert Q.max() > 0