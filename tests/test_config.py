from shs.config.simulation import load_simulation


def test_simulation_config():

    config = load_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    assert config.temperature == 3.0
    assert config.current == 0.001
    assert config.dt == 1e-9