"""Conservation checks for electrical and thermal flow around a hole."""

import numpy as np
import pytest

from shs.config.builder import build_simulation
from shs.physics import ThermalModel
from shs.solvers.electrical_solver import electrical_simulation_step
from shs.solvers.thermal_solver import thermal_step


@pytest.fixture
def ring():
    return build_simulation("configs/simulations/nbn_ring_validation.json")


def _section_currents(simulation):
    weights = np.ones(simulation.mesh.ny)
    weights[[0, -1]] = 0.5
    return np.sum(
        simulation.fields.current_density_x[:, :-1]
        * simulation.material_map.thickness[:, :-1]
        * weights[:, None]
        * simulation.mesh.dy,
        axis=0,
    )


def test_voltage_driven_normal_current_flows_around_hole_symmetrically(ring):
    ring.config.electrical.drive_mode = "voltage"
    ring.config.electrical.voltage_left = 1e-6
    ring.config.electrical.voltage_right = 0.0
    ring.config.electrical.solver.backend = "sparse_direct"
    ring.config.electrical.solver.tolerance = 1e-11
    electrical_simulation_step(ring, include_superconductivity=False)

    active = ring.region_map.active_mask
    assert np.all(ring.fields.normal_current_density_x[~active] == 0.0)
    assert np.all(ring.fields.normal_current_density_y[~active] == 0.0)
    assert np.allclose(ring.fields.voltage, ring.fields.voltage[::-1, :], atol=2e-16)
    # The legacy Dirichlet solve conserves its nodal (rectangle-rule) flux;
    # the current-controlled finite-volume test below validates physical
    # trapezoidal cross-section integration.
    nodal_sections = np.sum(
        ring.fields.current_density_x[:, :-1], axis=0
    ) * ring.geometry.film.thickness * ring.mesh.dy
    assert np.ptp(nodal_sections) / np.mean(np.abs(nodal_sections)) < 2e-12
    center = ring.mesh.nx // 2
    upper = np.sum(np.abs(ring.fields.current_density_y[:center, center]))
    lower = np.sum(np.abs(ring.fields.current_density_y[center + 1:, center]))
    assert upper == pytest.approx(lower, rel=2e-12)
    assert upper > 0.0


def test_current_driven_perforated_domain_recovers_requested_current(ring):
    requested = 1e-8
    ring.config.current = requested
    ring.config.electrical.drive_mode = "current"
    ring.config.electrical.solver.backend = "sparse_direct"
    ring.config.electrical.solver.tolerance = 1e-11
    electrical_simulation_step(ring, include_superconductivity=False)

    sections = _section_currents(ring)
    assert sections == pytest.approx(requested, rel=3e-12)
    assert ring.fields.electrical_solver_residual < 1e-11


def test_insulating_hole_blocks_thermal_leakage(ring):
    active = ring.region_map.active_mask
    ring.boundaries = None
    ring.fields.temperature.fill(4.0)
    ring.fields.temperature[~active] = 1e6
    ring.fields.heat_source.fill(0.0)
    thermal_step(
        ring,
        1e-13,
        ThermalModel(bath_temperature=4.0, thermal_relaxation_rate=0.0,
                     max_substep=1e-13),
    )
    assert np.all(ring.fields.temperature[active] == 4.0)
    assert np.all(ring.fields.temperature[~active] == 1e6)


def test_perforated_thermal_diffusion_conserves_energy_and_symmetry(ring):
    active = ring.region_map.active_mask
    ring.boundaries = None
    yy, xx = np.indices(active.shape)
    initial = 4.0 + 2.0 * np.exp(-((xx - 24.0) ** 2 + (yy - 50.0) ** 2) / 80.0)
    ring.fields.temperature = np.where(active, initial, 4.0)
    ring.fields.heat_source.fill(0.0)
    wx = np.ones(ring.mesh.nx)
    wy = np.ones(ring.mesh.ny)
    wx[[0, -1]] = wy[[0, -1]] = 0.5
    weights = active * wy[:, None] * wx[None, :]
    capacity = ring.material_map.heat_capacity
    energy_before = np.sum(weights * capacity * ring.fields.temperature)

    thermal_step(
        ring,
        2e-12,
        ThermalModel(bath_temperature=4.0, thermal_relaxation_rate=0.0,
                     max_substep=1e-13),
    )
    energy_after = np.sum(weights * capacity * ring.fields.temperature)
    assert energy_after == pytest.approx(energy_before, rel=3e-15)
    assert np.allclose(ring.fields.temperature, ring.fields.temperature[::-1], atol=2e-14)
    # Heat reaches both paths around the insulating obstacle equally.
    center = ring.mesh.nx // 2
    assert ring.fields.temperature[35, center] == pytest.approx(
        ring.fields.temperature[65, center], abs=2e-14
    )
