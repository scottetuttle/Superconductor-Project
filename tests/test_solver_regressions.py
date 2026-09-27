"""Regression tests for physical time, units, conservation, boundaries and failure isolation."""
from copy import deepcopy
import numpy as np
import pytest
from shs.config.builder import build_simulation
from shs.geometry.mesh import create_mesh
from shs.mapping import build_region_map, build_material_map, build_contact_map
from shs.mapping.contact_map import ContactMap
from shs.physics.fields import Fields
from shs.physics.thermal import ThermalModel
from shs.boundaries import BoundarySet, BoundaryCondition, BoundarySide, BoundaryType
from shs.solvers.coupled_solver import coupled_step, run_coupled_simulation, build_tdgl_model
from shs.solvers.electrical_solver import electrical_step, link_divergence
from shs.solvers.thermal_solver import thermal_step
from shs.solvers.tdgl_solver import tdgl_step, tdgl_scales
from shs.numerics.iterative import red_black_sor, sparse_direct
from shs.optics.hotspot import GaussianHotspot
from shs.physics.diagnostics import finite_volume_current_divergence


@pytest.fixture
def small():
    sim = build_simulation('configs/simulations/nbn_hotspot_test.json')
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
    sim.config.dt = 1e-14
    return sim


def test_coupling_iterations_do_not_advance_time(small):
    model = build_tdgl_model(small)
    initial = deepcopy(small)
    expected = deepcopy(small)
    dt = model.characteristic_time(15.5)*.01
    tdgl_step(expected, dt, model)
    for minimum in (1, 5, 12):
        sim = deepcopy(initial)
        sim.config.coupling.minimum_iterations = minimum
        _, result = coupled_step(sim, dt=dt, relaxation=.7)
        assert result.converged
        assert result.iterations >= minimum
        assert np.allclose(sim.fields.psi, expected.fields.psi, atol=1e-12)
        assert np.allclose(sim.fields.temperature, 3)


def test_coupled_thermal_endpoint_matches_analytical_cooling(small):
    small.fields.temperature.fill(10)
    small.fields.psi.fill(0)
    model = ThermalModel(3, thermal_relaxation_rate=1e10, max_substep=1e-15)
    dt = 1e-12
    _, result = coupled_step(small, dt=dt, thermal_model=model)
    assert result.converged
    assert np.allclose(small.fields.temperature, 3+7*np.exp(-1e10*dt), atol=1e-6)


def test_failed_coupling_preserves_state(small):
    small.config.electrical.voltage_left = .001
    small.config.coupling.max_iterations = 1
    before = deepcopy(vars(small.fields))
    identity = small.fields
    _, result = coupled_step(small)
    assert not result.converged
    assert small.fields is identity
    for name, value in before.items():
        assert np.array_equal(getattr(small.fields, name), value)


def test_exception_preserves_state(small, monkeypatch):
    import shs.solvers.coupled_solver as module
    before = deepcopy(small.fields)
    def fail(simulation, *args, **kwargs):
        simulation.fields.temperature.fill(999)
        raise RuntimeError('injected failure')
    monkeypatch.setattr(module, 'thermal_step', fail)
    with pytest.raises(RuntimeError, match='injected'):
        coupled_step(small)
    assert np.array_equal(small.fields.temperature, before.temperature)
    assert np.array_equal(small.fields.psi, before.psi)


def test_duration_uses_final_partial_step(small):
    small.fields.psi.fill(0)
    small.config.duration = 2.5e-14
    result = run_coupled_simulation(small)
    assert result.converged and result.steps == 3
    assert result.elapsed_time == pytest.approx(2.5e-14)


def test_physical_supercurrent_scale(small):
    model = build_tdgl_model(small)
    scales = tdgl_scales(small, model)
    a = .002
    small.fields.vector_potential_x.fill(a*scales.vector_potential_scale)
    tdgl_step(small, scales.time_scale*1e-8, model)
    # Uniform interior: Im(conj(psi)*(exp(-i A' dx')-1)*psi)/dx'.
    dx = small.mesh.dx/scales.xi
    expected = scales.supercurrent_density_to_physical(-np.sin(a*dx)/dx)
    assert np.allclose(small.fields.supercurrent_density_x[2:-2, 2:-2], expected, rtol=1e-7)
    assert np.allclose(small.fields.vector_potential_x, a*scales.vector_potential_scale)


def test_normal_current_source_sign_and_conservation(small):
    fields, mesh = small.fields, small.mesh
    # Manufactured potential and Js = sigma grad(V), so total current is zero.
    x, y = np.meshgrid(mesh.x, mesh.y)
    exact = 1e-3*(x/mesh.x[-1])**2
    jsx, jsy = np.zeros_like(x), np.zeros_like(x)
    sigma = small.material_map.electrical_conductivity
    jsx[:, :-1] = sigma[:, :-1]*np.diff(exact, axis=1)/mesh.dx
    left, right = np.zeros_like(x, dtype=bool), np.zeros_like(x, dtype=bool)
    left[:, 0], right[:, -1] = True, True
    small.contact_map.contact_masks.update(left_current=left, right_current=right)
    electrical_step(fields, mesh, small.material_map, small.contact_map,
                    voltage_left=0, voltage_right=.001, superconducting_current_x=jsx,
                    superconducting_current_y=jsy, solver_tolerance=1e-14)
    assert np.allclose(fields.voltage, exact, atol=1e-11)
    divergence = link_divergence(fields.current_density_x, fields.current_density_y, mesh.dx, mesh.dy)
    assert np.max(np.abs(divergence[:, 1:-1])) < 1e-6*np.max(np.abs(jsx))/mesh.dx


def test_direct_current_flux_drive_matches_requested_cross_section_current(small):
    fields, mesh = small.fields, small.mesh
    fields.psi.fill(0.0)
    requested = 2.5e-8

    electrical_step(
        fields, mesh, small.material_map, small.contact_map,
        drive_mode="current", applied_current=requested,
        solver_backend="sparse_direct", solver_tolerance=1e-12,
    )

    thickness = small.material_map.thickness[:, mesh.nx // 2]
    measured = np.trapezoid(
        fields.current_density_x[:, mesh.nx // 2] * thickness,
        mesh.y,
    )
    assert measured == pytest.approx(requested, rel=1e-11)
    assert np.ptp(fields.voltage[:, mesh.nx // 2]) < 1e-12
    divergence = link_divergence(
        fields.current_density_x, fields.current_density_y, mesh.dx, mesh.dy
    )
    assert np.max(np.abs(divergence[:, 1:-1])) < 1e-10 * requested / (
        np.mean(thickness) * mesh.dx * mesh.dy
    )
    finite_volume_divergence = finite_volume_current_divergence(
        fields.current_density_x, fields.current_density_y,
        mesh, small.material_map.thickness,
    )
    assert np.max(np.abs(finite_volume_divergence[:, 1:-1])) < (
        1e-10 * requested / (np.mean(thickness) * mesh.dx * mesh.dy)
    )


def test_direct_current_flux_drive_supports_matching_partial_contacts(small):
    left = np.zeros_like(small.fields.voltage, dtype=bool)
    right = np.zeros_like(left)
    left[2:5, 0] = True
    right[2:5, -1] = True
    small.contact_map.contact_masks.update(left_current=left, right_current=right)

    electrical_step(
        small.fields, small.mesh, small.material_map, small.contact_map,
        drive_mode="current", applied_current=1e-8,
        solver_backend="sparse_direct", solver_tolerance=1e-11,
    )

    center = small.mesh.nx // 2
    measured = np.trapezoid(
        small.fields.current_density_x[:, center]
        * small.material_map.thickness[:, center],
        small.mesh.y,
    )
    assert measured == pytest.approx(1e-8, rel=2e-3)


def test_direct_current_flux_drive_preserves_reversed_current_sign(small):
    requested = -1.25e-8
    electrical_step(
        small.fields, small.mesh, small.material_map, small.contact_map,
        drive_mode="current", applied_current=requested,
        solver_backend="sparse_direct", solver_tolerance=1e-12,
    )
    column = small.mesh.nx // 2
    measured = np.trapezoid(
        small.fields.current_density_x[:, column]
        * small.material_map.thickness[:, column],
        small.mesh.y,
    )
    assert measured == pytest.approx(requested, rel=1e-11)


def test_direct_current_flux_drive_supports_bottom_to_top_terminals(small):
    bottom = np.zeros_like(small.fields.voltage, dtype=bool)
    top = np.zeros_like(bottom)
    bottom[0, :] = True
    top[-1, :] = True
    small.contact_map = ContactMap(
        {"source": bottom, "sink": top},
        {"source": "current", "sink": "current"},
    )
    requested = 7.5e-9
    electrical_step(
        small.fields, small.mesh, small.material_map, small.contact_map,
        drive_mode="current", applied_current=requested,
        source_contact="source", sink_contact="sink",
        solver_backend="sparse_direct", solver_tolerance=1e-12,
    )
    row = small.mesh.ny // 2
    measured = np.trapezoid(
        small.fields.current_density_y[row, :]
        * small.material_map.thickness[row, :],
        small.mesh.x,
    )
    assert measured == pytest.approx(requested, rel=1e-11)


def test_coupled_step_converges_with_self_consistent_screening(small):
    small.config.electrical.drive_mode = "current"
    small.config.electrical.solver.backend = "sparse_direct"
    small.config.current = 1e-8
    magnetic = small.config.electromagnetic
    magnetic.include_self_field = True
    magnetic.screening_tolerance = 1e-3
    magnetic.screening_max_iterations = 60
    magnetic.screening_step_size = 0.5

    _, result = coupled_step(small, dt=1e-14)

    assert result.converged
    assert small.fields.magnetic_screening_residual <= magnetic.screening_tolerance
    assert np.any(np.abs(small.fields.induced_vector_potential_x) > 0.0)
    assert np.allclose(
        small.fields.vector_potential_x,
        small.fields.applied_vector_potential_x
        + small.fields.induced_vector_potential_x,
    )


def test_poisson_insulating_edges_and_tiny_omega():
    shape = (7, 11)
    fixed, values = np.zeros(shape, bool), np.zeros(shape)
    fixed[:, 0], fixed[:, -1], values[:, -1] = True, True, 1
    args = (np.zeros(shape), np.ones(shape), np.zeros(shape), fixed, values, 1., 1.)
    result = red_black_sor(*args, tolerance=1e-12)
    assert result.converged
    assert np.allclose(result.field, np.linspace(0, 1, 11)[None, :], atol=1e-10)
    result = red_black_sor(*args, omega=1e-12, max_iterations=1, tolerance=1e-8)
    assert not result.converged  # tiny iterate changes cannot masquerade as a solution


def test_sparse_direct_matches_sor_equation_and_boundary_conditions():
    shape = (9, 13)
    coefficient = np.linspace(1.0, 3.0, shape[1])[None, :] * np.ones(shape)
    fixed, values = np.zeros(shape, bool), np.zeros(shape)
    fixed[:, 0], fixed[:, -1] = True, True
    values[:, -1] = 0.7
    source = np.zeros(shape)
    args = (np.zeros(shape), coefficient, source, fixed, values, 0.3, 0.4)
    sor = red_black_sor(*args, tolerance=1e-11, max_iterations=100000)
    direct = sparse_direct(*args, tolerance=1e-11)
    assert sor.converged and direct.converged
    assert direct.residual <= 1e-11
    assert np.allclose(direct.field, sor.field, atol=2e-10)


def test_electrical_failure_does_not_commit(small):
    before = small.fields.voltage.copy()
    with pytest.raises(RuntimeError, match='Electrical solve failed'):
        electrical_step(small.fields, small.mesh, small.material_map, small.contact_map,
                        solver_max_iterations=1)
    assert np.array_equal(small.fields.voltage, before)


def test_external_heat_is_not_overwritten_or_accumulated(small):
    small.fields.external_heat_source = np.full_like(small.fields.temperature, 1e12)
    args = (small.fields, small.mesh, small.material_map, small.contact_map)
    electrical_step(*args, voltage_left=0, voltage_right=0)
    electrical_step(*args, voltage_left=0, voltage_right=0)
    assert np.all(small.fields.heat_source == 1e12)
    thermal_step(small, 1e-12, ThermalModel(3, max_substep=1e-12))
    assert np.allclose(small.fields.temperature, 3+1e-12*1e12/2e6, atol=1e-12)


def test_fixed_temperature_boundary(small):
    small.boundaries.add(BoundaryCondition(BoundarySide.LEFT, BoundaryType.FIXED_TEMPERATURE, temperature=5))
    thermal_step(small, 1e-10, ThermalModel(3, max_substep=1e-10))
    assert np.all(small.fields.temperature[:, 0] == 5)
    assert np.all(small.fields.temperature[:, 1] > 3)


def test_insulated_thermal_energy_conservation(small):
    small.fields.temperature[3, 4] = 12
    small.material_map.thermal_conductivity[:, 4:] *= 3
    weights = np.ones_like(small.fields.temperature)
    weights[:, [0, -1]] *= .5
    weights[[0, -1]] *= .5
    initial = np.sum(weights*small.fields.temperature*small.material_map.heat_capacity)
    thermal_step(small, 1e-9, ThermalModel(3, max_substep=1e-9))
    final = np.sum(weights*small.fields.temperature*small.material_map.heat_capacity)
    assert final == pytest.approx(initial, rel=1e-13)
    assert small.fields.temperature.min() >= 3
    assert small.fields.temperature.max() <= 12


@pytest.mark.parametrize('kind', [BoundaryType.FIXED_HEAT_FLUX, BoundaryType.CONVECTION])
def test_boundary_flux_energy_balance(small, kind):
    flux = 4.0
    boundary = BoundaryCondition(BoundarySide.LEFT, kind, heat_flux=flux,
                                  temperature=5., heat_transfer=2.)
    small.boundaries.add(boundary)
    initial = small.fields.temperature.copy()
    dt = 1e-12
    thermal_step(small, dt, ThermalModel(3, max_substep=dt))
    weights = np.ones_like(initial)
    weights[:, [0, -1]] *= .5
    weights[[0, -1]] *= .5
    energy = np.sum(weights*(small.fields.temperature-initial)*small.material_map.heat_capacity)*small.mesh.dx*small.mesh.dy
    assert energy == pytest.approx(flux*small.mesh.y[-1]*dt, rel=1e-4)


def test_gaussian_hotspot_rectangular_mesh(small):
    heat = GaussianHotspot(small.mesh.x[3], small.mesh.y[2], 10, 1e-7).generate(small.mesh)
    assert heat.shape == (7, 9)
    assert np.unravel_index(np.argmax(heat), heat.shape) == (2, 3)


def test_tdgl_seconds_guard(small):
    with pytest.raises(ValueError, match='seconds'):
        tdgl_step(small, .001, build_tdgl_model(small))

@pytest.mark.parametrize('filename', ['coupled_convergence_benchmark.py',
    'benchmark_adaptive_coupled_updated.py'])
def test_benchmark_external_heat_adapter(small, filename):
    import importlib.util
    import sys
    name = filename.removesuffix('.py')
    spec = importlib.util.spec_from_file_location(name, 'tools/'+filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    source = np.full_like(small.fields.temperature, 1e12)
    model = ThermalModel(3, max_substep=1e-12)
    module.run_thermal(small, 1e-12, model, source)
    module.run_thermal(small, 1e-12, model, source)
    assert np.allclose(small.fields.temperature, 3+2e-12*1e12/2e6, atol=1e-12)
    assert np.all(small.fields.heat_source == source)


def test_shipped_configuration_couples_with_consistent_tolerances():
    sim = build_simulation('configs/simulations/nbn_hotspot_test.json')
    sim.config.coupling.max_iterations = 10
    result = run_coupled_simulation(sim, steps=2, dt=1e-14)
    assert result.converged and result.elapsed_time == pytest.approx(2e-14)
    assert result.total_coupling_iterations <= 10
    assert np.all(
        sim.fields.temperature[:, 0]
        == sim.config.boundaries["left"]["temperature"]
    )


def test_thermal_diffusion_cosine_benchmark(small):
    length = small.mesh.x[-1]
    initial_mode = .5*np.cos(np.pi*small.mesh.x/length)
    small.fields.temperature[:] = 3+initial_mode
    dt = 1e-7
    thermal_step(small, dt, ThermalModel(3, max_substep=1e-9))
    expected = 3+initial_mode*np.exp(-(10/2e6)*(np.pi/length)**2*dt)
    assert np.max(np.abs(small.fields.temperature-expected)) < .002
