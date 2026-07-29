import numpy as np

from shs.config.builder import build_simulation
from shs.physics import ThermalModel
from shs.solvers import thermal_step


CONFIG = "configs/simulations/nbn_hotspot_test.json"


def create_test_simulation():

    return build_simulation(CONFIG)


def test_heat_diffusion():

    simulation = create_test_simulation()

    fields = simulation.fields

    fields.temperature[50,50] = 10.0


    thermal = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=0.0
    )


    updated_fields = thermal_step(
        simulation,
        dt=1e-9,
        thermal_model=thermal
    )


    assert updated_fields.temperature.shape == (
        simulation.mesh.ny,
        simulation.mesh.nx
    )


    assert (
        updated_fields.temperature[50,50]
        <=
        10.0
    )


def test_boundary_stability():

    simulation = create_test_simulation()

    fields = simulation.fields

    fields.heat_source[50,50] = 10.0


    thermal = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=0.0
    )


    updated = thermal_step(
        simulation,
        dt=1e-6,
        thermal_model=thermal
    )


    assert updated.temperature[0,0] == 3.0



def test_hotspot_heating():

    simulation = create_test_simulation()

    fields = simulation.fields

    fields.heat_source[50,50] = 1e12


    thermal = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=0.0
    )


    updated_fields = thermal_step(
        simulation,
        dt=1e-9,
        thermal_model=thermal
    )


    assert updated_fields.temperature[50,50] > 3.0



def test_uniform_temperature_remains_constant():

    simulation = create_test_simulation()

    fields = simulation.fields


    thermal = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=0.0
    )


    updated_fields = thermal_step(
        simulation,
        dt=1e-6,
        thermal_model=thermal
    )


    assert np.allclose(
        updated_fields.temperature,
        fields.temperature
    )



def test_bath_cooling():

    simulation = create_test_simulation()

    fields = simulation.fields

    fields.temperature[:] = 10.0


    thermal = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=100.0
    )


    updated_fields = thermal_step(
        simulation,
        dt=1e-6,
        thermal_model=thermal
    )


    assert (
        updated_fields.temperature[50,50]
        <
        10.0
    )



def test_bath_equilibrium():

    simulation = create_test_simulation()

    fields = simulation.fields


    thermal = ThermalModel(
        bath_temperature=3.0,
        thermal_relaxation_rate=100.0
    )


    updated_fields = thermal_step(
        simulation,
        dt=1e-6,
        thermal_model=thermal
    )


    assert np.allclose(
        updated_fields.temperature,
        3.0
    )