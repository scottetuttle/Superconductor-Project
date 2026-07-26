import numpy as np

from shs.geometry import load_geometry, create_mesh
from shs.physics import ThermalModel
from shs.solvers import thermal_step
from shs.physics import Fields



def test_heat_diffusion():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)


    fields = Fields.create(
        mesh,
        initial_temperature=3.0
)

    # Create a hot spot
    fields.temperature[50,50] = 10.0

    thermal = ThermalModel(
        thermal_conductivity=10,
        heat_capacity=2e6,
        bath_temperature=3,
        thermal_relaxation_rate=0.0
    )

    updated_fields = thermal_step(
        fields,
        mesh,
        thermal,
        dt=1e-9
    )


    assert (
        updated_fields.temperature.shape
        ==
        fields.temperature.shape
)

    # Hot spot should cool
    assert (
        updated_fields.temperature[50,50]
        <=
        fields.temperature[50,50]
)


    
def test_boundary_stability():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)


    fields = Fields.create(
        mesh,
        initial_temperature=3.0
)

    fields.heat_source[50,50] = 10.0


    thermal = ThermalModel(
        thermal_conductivity=10,
        heat_capacity=2e6,
        bath_temperature=3,
        thermal_relaxation_rate=0.0
    )


    updated = thermal_step(
        fields,
        mesh,
        thermal,
        dt=1e-6
    )


    # Edges should remain close to bath temperature
    assert updated.temperature[0,0] == 3.0

def test_hotspot_heating():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)


    fields = Fields.create(
        mesh,
        initial_temperature=3.0
    )


    fields.heat_source[50,50] = 1e12


    thermal = ThermalModel(
        thermal_conductivity=10,
        heat_capacity=2e6,
        bath_temperature=3,
        thermal_relaxation_rate=0.0
    )


    updated_fields = thermal_step(
        fields,
        mesh,
        thermal,
        dt=1e-9
    )


    assert updated_fields.temperature[50,50] > 3.0

def test_uniform_temperature_remains_constant():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    fields = Fields.create(
        mesh,
        initial_temperature=3.0
    )

    thermal = ThermalModel(
        thermal_conductivity=10,
        heat_capacity=2e6,
        bath_temperature=3,
        thermal_relaxation_rate=0.0
    )

    updated_fields = thermal_step(
        fields,
        mesh,
        thermal,
        dt=1e-6
    )

    assert np.allclose(
        updated_fields.temperature,
        fields.temperature
    )

def test_bath_cooling():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    fields = Fields.create(
        mesh,
        initial_temperature=10.0
    )

    thermal = ThermalModel(
        thermal_conductivity=10,
        heat_capacity=2e6,
        bath_temperature=3.0,
        thermal_relaxation_rate=100.0
    )

    updated_fields = thermal_step(
        fields,
        mesh,
        thermal,
        dt=1e-6
    )

    assert (
        updated_fields.temperature[50, 50]
        <
        10.0
    )

def test_bath_equilibrium():

    geometry = load_geometry(
        "configs/geometry/NbN_film.json"
    )

    mesh = create_mesh(geometry)

    fields = Fields.create(
        mesh,
        initial_temperature=3.0
    )

    thermal = ThermalModel(
        thermal_conductivity=10,
        heat_capacity=2e6,
        bath_temperature=3.0,
        thermal_relaxation_rate=100.0
    )

    updated_fields = thermal_step(
        fields,
        mesh,
        thermal,
        dt=1e-6
    )

    assert np.allclose(
        updated_fields.temperature,
        3.0
    )