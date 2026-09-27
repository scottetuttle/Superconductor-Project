import numpy as np
import pytest

from shs.config.builder import build_simulation
from shs.config.simulation import LaserWaypointConfig
from shs.optics.moving_laser import (
    apply_laser_heat_source, gaussian_laser_heat_source, laser_position,
    laser_power_fraction,
)
from shs.physics import ThermalModel
from shs.solvers.thermal_solver import thermal_step
from tools.tdgl_diagnostics import track_vortices


def _enable_laser(simulation, power=2e-7, sigma=3e-9):
    config = simulation.config.laser
    config.enabled = True
    config.absorbed_power_W = power
    config.sigma_m = sigma
    config.waypoints = [
        LaserWaypointConfig(0.0, 1e-8, 2e-8),
        LaserWaypointConfig(2e-12, 4e-8, 3e-8),
    ]
    config.validate()
    return config


def test_laser_path_interpolates_physical_position():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    config = _enable_laser(simulation)
    midpoint = laser_position(config, 1e-12)
    assert midpoint.x_m == pytest.approx(2.5e-8)
    assert midpoint.y_m == pytest.approx(2.5e-8)
    assert laser_position(config, 3e-12).x_m == pytest.approx(4e-8)


def test_laser_protocol_interpolates_absorbed_power_envelope():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    config = _enable_laser(simulation, power=2e-7, sigma=2e-9)
    config.waypoints = [
        LaserWaypointConfig(0.0, 2.5e-8, 2.5e-8, 0.0),
        LaserWaypointConfig(2e-12, 2.5e-8, 2.5e-8, 1.0),
    ]
    assert laser_power_fraction(config, 1e-12) == pytest.approx(0.5)
    heat, _ = gaussian_laser_heat_source(simulation, 1e-12)
    deposited = np.sum(
        heat * simulation.material_map.thickness
    ) * simulation.mesh.dx * simulation.mesh.dy
    assert deposited == pytest.approx(0.5 * config.absorbed_power_W, rel=2e-6)


def test_laser_rejects_invalid_power_fraction():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    config = _enable_laser(simulation)
    config.waypoints = [LaserWaypointConfig(0.0, 1e-8, 1e-8, 1.1)]
    with pytest.raises(ValueError, match="power_fraction"):
        config.validate()


def test_centered_gaussian_deposits_configured_power_on_full_film():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    config = _enable_laser(simulation, sigma=2e-9)
    config.waypoints = [LaserWaypointConfig(0.0, 2.5e-8, 2.5e-8)]
    heat, _ = gaussian_laser_heat_source(simulation, 0.0)
    deposited = np.sum(
        heat * simulation.material_map.thickness
    ) * simulation.mesh.dx * simulation.mesh.dy
    assert deposited == pytest.approx(config.absorbed_power_W, rel=2e-6)


def test_power_over_ring_hole_is_not_renormalized_into_film():
    simulation = build_simulation("configs/simulations/nbn_ring_validation.json")
    config = _enable_laser(simulation, sigma=3e-9)
    config.waypoints = [LaserWaypointConfig(0.0, 2.5e-8, 2.5e-8)]
    heat, _ = gaussian_laser_heat_source(simulation, 0.0)
    deposited = np.sum(
        heat * simulation.material_map.thickness
    ) * simulation.mesh.dx * simulation.mesh.dy
    assert deposited < 0.1 * config.absorbed_power_W
    assert np.all(heat[~simulation.region_map.active_mask] == 0.0)


def test_laser_energy_enters_active_thermal_domain():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    config = _enable_laser(simulation, sigma=2e-9)
    config.waypoints = [LaserWaypointConfig(0.0, 2.5e-8, 2.5e-8)]
    simulation.boundaries = None
    simulation.material_map.thermal_conductivity.fill(0.0)
    simulation.fields.temperature.fill(14.0)
    apply_laser_heat_source(simulation, 0.0)
    deposited = np.sum(
        simulation.fields.laser_heat_source * simulation.material_map.thickness
    ) * simulation.mesh.dx * simulation.mesh.dy
    dt = 1e-13
    initial_energy = np.sum(
        simulation.fields.temperature * simulation.material_map.heat_capacity
        * simulation.material_map.thickness
    ) * simulation.mesh.dx * simulation.mesh.dy
    thermal_step(simulation, dt, ThermalModel(14.0, 0.0, dt))
    final_energy = np.sum(
        simulation.fields.temperature * simulation.material_map.heat_capacity
        * simulation.material_map.thickness
    ) * simulation.mesh.dx * simulation.mesh.dy
    assert final_energy - initial_energy == pytest.approx(deposited * dt, rel=2e-11)


def test_vortex_tracker_preserves_sign_and_identity():
    frames, rows = [], []
    for index in range(3):
        winding = np.zeros((10, 10), dtype=int)
        winding[2 + index, 2 + index] = 1
        winding[7 - index, 7] = -1
        frames.append({"winding": winding, "mesh_dx_m": 1e-9, "mesh_dy_m": 1e-9})
        rows.append({"step": index, "time_s": index * 1e-14})
    records = track_vortices(frames, rows, max_displacement_cells=2.0)
    assert len(records) == 6
    assert len({row["track_id"] for row in records}) == 2
    for track_id in {row["track_id"] for row in records}:
        track = [row for row in records if row["track_id"] == track_id]
        assert len(track) == 3
        assert len({row["charge"] for row in track}) == 1
