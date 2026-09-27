import json
from copy import deepcopy

import numpy as np
import pytest

from shs.config.builder import build_simulation
from shs.config.simulation import LaserWaypointConfig, PinningSiteConfig
from shs.optics.moving_laser import gaussian_laser_heat_source
from shs.physics.pinning import pinning_suppression
from shs.solvers.coupled_solver import build_tdgl_model, build_thermal_model


CONFIG = "configs/simulations/nbn_tdgl_laser_protocol_demo.json"


def test_simulation_json_controls_constructed_geometry_material_and_state():
    simulation = build_simulation(CONFIG)
    with open(CONFIG, encoding="utf-8") as handle:
        source = json.load(handle)
    assert simulation.geometry.film.name == "NbN_ring"
    assert simulation.mesh.nx == 101
    assert simulation.mesh.ny == 101
    assert simulation.material_map.materials[0].name == "NbN_tdgl_reference"
    assert np.all(simulation.fields.temperature == source["temperature"])
    assert simulation.config.current == source["current"]["value"]


def test_solver_models_receive_thermal_and_tdgl_configuration_values():
    simulation = build_simulation(CONFIG)
    thermal = build_thermal_model(simulation)
    tdgl = build_tdgl_model(simulation)
    assert thermal.bath_temperature == simulation.config.thermal.bath_temperature
    assert thermal.thermal_relaxation_rate == simulation.config.thermal.thermal_relaxation_rate
    assert tdgl.parameters.u == simulation.config.tdgl.u
    assert tdgl.parameters.gamma == simulation.config.tdgl.gamma
    assert tdgl.parameters.include_scalar_potential == simulation.config.tdgl.include_scalar_potential


def test_laser_configuration_changes_source_position_power_and_width():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    laser = simulation.config.laser
    laser.enabled = True
    laser.absorbed_power_W = 2e-7
    laser.sigma_m = 2e-9
    laser.waypoints = [LaserWaypointConfig(0.0, 25e-9, 25e-9, 0.5)]
    narrow, _ = gaussian_laser_heat_source(simulation, 0.0)
    deposited = np.sum(narrow * simulation.material_map.thickness) * simulation.mesh.dx * simulation.mesh.dy
    assert deposited == pytest.approx(1e-7, rel=2e-6)
    laser.sigma_m = 5e-9
    broad, _ = gaussian_laser_heat_source(simulation, 0.0)
    assert broad.max() < narrow.max()


def test_pinning_configuration_changes_tdgl_coefficient_landscape():
    simulation = build_simulation(CONFIG)
    simulation.config.pinning.enabled = False
    assert np.max(pinning_suppression(simulation)) == 0
    simulation.config.pinning.enabled = True
    simulation.config.pinning.sites = [PinningSiteConfig(40e-9, 25e-9, 2e-9, 0.3)]
    assert np.max(pinning_suppression(simulation)) == pytest.approx(0.3)


def test_source_config_is_independent_from_mutable_loaded_config():
    first = build_simulation(CONFIG)
    second = build_simulation(CONFIG)
    original = deepcopy(second.config.laser.waypoints)
    first.config.laser.waypoints.clear()
    assert second.config.laser.waypoints == original
