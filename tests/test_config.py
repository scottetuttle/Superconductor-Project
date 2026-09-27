from shs.config.simulation import load_simulation


def test_simulation_config():

    config = load_simulation(
        "configs/simulations/nbn_hotspot_test.json"
    )

    assert config.temperature == 9.0
    assert config.current == 0.01
    assert config.dt == 1e-9


def test_reference_transport_uses_consistent_tdgl_models():
    config = load_simulation("configs/simulations/nbn_tdgl_transport.json")

    assert config.tdgl.normalization == "pytdgl"
    assert config.tdgl.temperature_model == "tc_over_t_minus_one"
    assert config.electrical.normal_conductivity_model == "constant"
    assert config.electrical.drive_mode == "voltage"
    assert config.electrical.source_contact == "left_current"
    assert config.electrical.sink_contact == "right_current"
    assert not config.electromagnetic.include_self_field
    assert config.electromagnetic.screening_tolerance > 0.0
