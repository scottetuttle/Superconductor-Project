import numpy as np
import pytest

from shs.config.builder import build_simulation
from shs.geometry.mesh import Mesh
from shs.physics.electromagnetics import (
    induced_vector_potential_from_sheet_current,
    perpendicular_magnetic_field,
    uniform_perpendicular_vector_potential,
)
from shs.physics.fluxoid import rectangular_fluxoid
from shs.solvers.coupled_solver import build_tdgl_model
from shs.solvers.magnetic_solver import magnetic_screening_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM


def _set_supercurrent(simulation):
    model = build_tdgl_model(simulation)
    scales = tdgl_scales(simulation, model)
    ax = scales.vector_potential_to_dimensionless(
        simulation.fields.vector_potential_x
    )
    ay = scales.vector_potential_to_dimensionless(
        simulation.fields.vector_potential_y
    )
    jx, jy = model.supercurrent_density(
        simulation.fields.psi, ax, ay,
        simulation.mesh.dx / scales.xi,
        simulation.mesh.dy / scales.xi,
    )
    simulation.fields.supercurrent_density_x = (
        scales.supercurrent_density_to_physical(jx)
    )
    simulation.fields.supercurrent_density_y = (
        scales.supercurrent_density_to_physical(jy)
    )
    simulation.fields.current_density_x = (
        simulation.fields.supercurrent_density_x.copy()
    )
    simulation.fields.current_density_y = (
        simulation.fields.supercurrent_density_y.copy()
    )
    return scales


def test_single_vortex_london_fluxoid_matches_phase_quantum():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    yy, xx = np.indices(simulation.fields.psi.shape)
    simulation.fields.psi = np.exp(1j * np.arctan2(yy - 49.5, xx - 49.5))
    _set_supercurrent(simulation)

    result = rectangular_fluxoid(simulation, (20, 79, 20, 79))

    assert result.winding == 1
    assert result.topological_fluxoid_Wb == pytest.approx(
        SUPERCONDUCTING_FLUX_QUANTUM
    )
    assert result.london_fluxoid_quanta == pytest.approx(1.0, rel=3e-4)


def test_fluxoid_and_current_are_invariant_under_discrete_gauge_transform():
    first = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    second = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    ax, ay = uniform_perpendicular_vector_potential(first.mesh, 0.05, "symmetric")
    for simulation in (first, second):
        simulation.fields.vector_potential_x = ax.copy()
        simulation.fields.vector_potential_y = ay.copy()
    yy, xx = np.indices(first.fields.psi.shape)
    chi = 0.2 * (xx / (first.mesh.nx - 1)) ** 2 + 0.1 * yy / (first.mesh.ny - 1)
    scales = tdgl_scales(second, build_tdgl_model(second))
    dx_dimensionless = second.mesh.dx / scales.xi
    dy_dimensionless = second.mesh.dy / scales.xi
    second.fields.psi *= np.exp(1j * chi)
    second.fields.vector_potential_x[:, :-1] += (
        np.diff(chi, axis=1) / dx_dimensionless * scales.vector_potential_scale
    )
    second.fields.vector_potential_y[:-1, :] += (
        np.diff(chi, axis=0) / dy_dimensionless * scales.vector_potential_scale
    )
    _set_supercurrent(first)
    _set_supercurrent(second)

    fluxoid_first = rectangular_fluxoid(first, (20, 79, 20, 79))
    fluxoid_second = rectangular_fluxoid(second, (20, 79, 20, 79))

    assert np.allclose(
        first.fields.supercurrent_density_x[:, :-1],
        second.fields.supercurrent_density_x[:, :-1], atol=1e-5,
    )
    assert np.allclose(
        first.fields.supercurrent_density_y[:-1, :],
        second.fields.supercurrent_density_y[:-1, :], atol=1e-5,
    )
    assert fluxoid_second.london_fluxoid_Wb == pytest.approx(
        fluxoid_first.london_fluxoid_Wb, abs=1e-20
    )
    for simulation, applied_x, applied_y in (
        (first, ax, ay),
        (second, second.fields.vector_potential_x.copy(),
         second.fields.vector_potential_y.copy()),
    ):
        simulation.fields.applied_vector_potential_x = applied_x.copy()
        simulation.fields.applied_vector_potential_y = applied_y.copy()
        simulation.config.electromagnetic.include_self_field = True
        simulation.config.electromagnetic.screening_step_size = 1.0
        magnetic_screening_step(simulation)
    assert np.allclose(
        first.fields.induced_vector_potential_x,
        second.fields.induced_vector_potential_x,
    )
    assert np.allclose(
        first.fields.induced_vector_potential_y,
        second.fields.induced_vector_potential_y,
    )
    assert np.allclose(first.fields.magnetic_field_z, second.fields.magnetic_field_z)


def test_first_screening_update_opposes_positive_applied_field_in_meissner_state():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    ax, ay = uniform_perpendicular_vector_potential(simulation.mesh, 0.01, "symmetric")
    fields = simulation.fields
    fields.applied_vector_potential_x = ax.copy()
    fields.applied_vector_potential_y = ay.copy()
    fields.vector_potential_x = ax.copy()
    fields.vector_potential_y = ay.copy()
    _set_supercurrent(simulation)
    simulation.config.electromagnetic.include_self_field = True
    simulation.config.electromagnetic.screening_step_size = 1.0

    magnetic_screening_step(simulation)

    induced_bz = perpendicular_magnetic_field(
        fields.induced_vector_potential_x,
        fields.induced_vector_potential_y,
        simulation.mesh.dx, simulation.mesh.dy,
    )
    center = (simulation.mesh.ny // 2, simulation.mesh.nx // 2)
    assert induced_bz[center] < 0.0
    assert fields.magnetic_field_z[center] < 0.01


def _uniform_sheet_current_center_potential(points):
    width = height = 50e-9
    x = np.linspace(0.0, width, points)
    y = np.linspace(0.0, height, points)
    mesh = Mesh(x, y, points, points, width / (points - 1), height / (points - 1))
    jx = np.full((points, points), 1e7)
    ax, _ = induced_vector_potential_from_sheet_current(
        mesh, jx, np.zeros_like(jx), 2e-9
    )
    return ax[points // 2, points // 2]


def test_induced_vector_potential_same_cell_regularization_refines():
    medium = _uniform_sheet_current_center_potential(51)
    fine = _uniform_sheet_current_center_potential(101)
    assert medium == pytest.approx(fine, rel=0.03)
