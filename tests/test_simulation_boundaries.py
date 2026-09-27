from shs.config.simulation import load_simulation
import json


def test_simulation_loads_boundaries():

    simulation = load_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    assert "left" in simulation.boundaries

    assert (
        simulation.boundaries["left"]["type"]
        == "fixed_temperature"
    )


def test_boundary_temperature():

    simulation = load_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    with open("configs/simulations/nbn_hotspot_test.json", encoding="utf-8") as file:
        expected = json.load(file)["boundaries"]["left"]["temperature"]

    assert simulation.boundaries["left"]["temperature"] == expected
