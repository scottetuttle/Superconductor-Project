"""Correctness checks for the experimental adaptive coupling benchmark."""

from copy import deepcopy
import importlib.util
import sys

import numpy as np
import pytest

from shs.physics.thermal import ThermalModel
from shs.config.builder import build_simulation
from shs.geometry.mesh import create_mesh
from shs.mapping import build_region_map, build_material_map
from shs.mapping.contact_map import ContactMap
from shs.physics.fields import Fields
from shs.boundaries import BoundarySet
from shs.physics.diagnostics import evaluate_physics_diagnostics


@pytest.fixture
def small():
    sim = build_simulation("configs/simulations/nbn_hotspot_test.json")
    sim.geometry.film.nx, sim.geometry.film.ny = 9, 7
    sim.mesh = create_mesh(sim.geometry)
    sim.region_map = build_region_map(sim.geometry, sim.mesh)
    sim.material_map = build_material_map(sim.region_map, sim.material_map.materials[0])
    left = np.zeros((7, 9), dtype=bool)
    right = np.zeros((7, 9), dtype=bool)
    left[:, 0] = True
    right[:, -1] = True
    sim.contact_map = ContactMap(
        {"left_current": left, "right_current": right},
        {"left_current": "current", "right_current": "current"},
    )
    sim.fields = Fields.create(sim.mesh, 3.0, np.ones((7, 9), dtype=complex))
    sim.boundaries = BoundarySet()
    sim.config.electrical.voltage_left = 0.0
    sim.config.thermal.bath_temperature = 3.0
    sim.config.thermal.thermal_relaxation_rate = 0.0
    sim.config.thermal.max_substep = 1e-10
    return sim


@pytest.fixture(scope="module")
def benchmark():
    name = "benchmark_adaptive_coupled_updated_test"
    spec = importlib.util.spec_from_file_location(
        name, "tools/benchmark_adaptive_coupled_updated.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_checkpoint_restores_derived_coupling_fields(small, benchmark):
    small.fields.supercurrent_density_x.fill(2.0)
    small.fields.joule_heat_source = np.full_like(small.fields.temperature, 3.0)
    small.fields.external_heat_source = np.full_like(small.fields.temperature, 5.0)
    small.fields.heat_source.fill(8.0)
    checkpoint = benchmark.copy_trial_state(small)

    small.fields.supercurrent_density_x.fill(20.0)
    small.fields.joule_heat_source.fill(30.0)
    small.fields.heat_source.fill(35.0)
    benchmark.apply_trial_state(small, checkpoint)

    assert np.all(small.fields.supercurrent_density_x == 2.0)
    assert np.all(small.fields.joule_heat_source == 3.0)
    assert np.all(small.fields.heat_source == 8.0)


def test_benchmark_thermal_step_uses_accepted_time_level(small, benchmark):
    accepted_temperature = np.full_like(small.fields.temperature, 3.0)
    small.fields.temperature.fill(10.0)
    benchmark.run_thermal(
        small,
        1e-12,
        ThermalModel(3.0, max_substep=1e-12),
        initial_temperature=accepted_temperature,
    )
    assert np.all(small.fields.temperature == 3.0)


def test_adapted_tdgl_model_preserves_physics_parameters(small, benchmark):
    small.config.tdgl.u = 7.0
    small.config.tdgl.gamma = 2.0
    small.config.tdgl.kappa = 4.0
    small.config.tdgl.stability_safety_factor = 0.6
    parameters = benchmark.NumericalParameters(1e-14, 0.5, 1e-7, 0.02, 1e-13)
    model = benchmark.build_tdgl_model(parameters, small)
    assert model.parameters.u == 7.0
    assert model.parameters.gamma == 2.0
    assert model.parameters.kappa == 4.0
    assert model.parameters.stability_safety_factor == 0.6
    assert model.parameters.max_normalized_timestep == 0.02


def test_benchmark_electrical_adapter_honors_iteration_budget(small, benchmark):
    small.fields.psi.fill(0.0)
    small.config.electrical.solver.backend = "red_black_sor"
    small.config.electrical.solver.max_iterations = 1
    small.config.electrical.solver.omega = 1.0
    with pytest.raises(RuntimeError, match="after 1 iterations"):
        benchmark.run_electrical(small, 1.0, 0.0, 1e-14)


def test_residual_is_mesh_independent_and_physically_scaled(benchmark):
    old_small = np.zeros((2, 2))
    new_small = np.full((2, 2), 0.5)
    old_large = np.zeros((20, 20))
    new_large = np.full((20, 20), 0.5)
    assert benchmark.relative_residual(old_small, new_small, scale=2.0) == pytest.approx(0.25)
    assert benchmark.relative_residual(old_large, new_large, scale=2.0) == pytest.approx(0.25)


def test_benchmark_heat_source_is_physically_centered(small, benchmark):
    heat = benchmark.make_gaussian_heat_source(small, 10.0, 2e-9)
    assert heat.shape == (7, 9)
    assert heat[3, 4] == pytest.approx(10.0)
    assert heat[3, 3] == pytest.approx(heat[3, 5])
    assert heat[2, 4] == pytest.approx(heat[4, 4])


def test_benchmark_heat_source_rejects_nonphysical_radius(small, benchmark):
    with pytest.raises(ValueError, match="positive finite meters"):
        benchmark.make_gaussian_heat_source(small, 10.0, 0.0)


def test_aitken_relaxation_accelerates_linear_fixed_point(benchmark):
    previous = np.array([1.0])
    current = np.array([0.5])
    sigma = benchmark.aitken_relaxation(previous, current, 0.5, 0.1, 1.0)
    assert sigma == pytest.approx(1.0)


def test_physics_diagnostics_measure_conservation_and_integrated_power(small):
    small.fields.current_density_x[:, :-1] = 4.0
    small.fields.current_density_y.fill(0.0)
    small.fields.joule_heat_source = np.full_like(small.fields.temperature, 2.0)
    result = evaluate_physics_diagnostics(small)
    expected_volume = (
        small.mesh.x[-1]
        * small.mesh.y[-1]
        * small.material_map.materials[0].thickness
    )
    assert result.current_continuity_rms == 0.0
    assert result.current_continuity_relative == 0.0
    assert result.total_joule_power == pytest.approx(2.0 * expected_volume)
