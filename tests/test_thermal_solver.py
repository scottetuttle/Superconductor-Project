import numpy as np

from shs.geometry import load_geometry, create_mesh
from shs.physics import ThermalModel
from shs.solvers import thermal_step



def test_heat_diffusion():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)


    temperature = np.ones(
        (
            mesh.nx,
            mesh.ny
        )
    ) * 3.0


    # Create a hot spot
    temperature[50,50] = 10.0


    thermal = ThermalModel(
        thermal_conductivity=10,
        heat_capacity=2e6,
        bath_temperature=3
    )


    updated = thermal_step(
        temperature,
        mesh,
        thermal,
        dt=1e-6
    )


    assert updated.shape == temperature.shape

    # Hot spot should cool
    assert updated[50,50] < temperature[50,50]