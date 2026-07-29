from shs.config.builder import build_simulation
import pytest


def test_invalid_material_map():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.material_map.Tc = simulation.material_map.Tc[:-1]

    with pytest.raises(ValueError):

        simulation.validate()

def test_simulation_validation():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    assert simulation.validate()


def test_invalid_field_shape():

    simulation = build_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    simulation.fields.temperature = (
        simulation.fields.temperature[:-1]
    )

    with pytest.raises(ValueError):

        simulation.validate()