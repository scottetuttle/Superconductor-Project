import runpy

import numpy as np
import pytest

from shs.config.builder import build_simulation
from shs.solvers.coupled_solver import build_tdgl_model
from shs.physics.electromagnetics import (
    induced_vector_potential_from_sheet_current,
    magnetic_field_from_sheet_current,
    perpendicular_magnetic_field,
    uniform_perpendicular_vector_potential,
)
from shs.solvers.magnetic_solver import magnetic_screening_step


DIAGNOSTICS = runpy.run_path("tools/tdgl_diagnostics.py", run_name="em_diagnostics_test")


@pytest.mark.parametrize("gauge", ["symmetric", "landau_x", "landau_y"])
def test_uniform_vector_potential_has_requested_curl(gauge):
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    ax, ay = uniform_perpendicular_vector_potential(simulation.mesh, 0.75, gauge)

    bz = perpendicular_magnetic_field(
        ax, ay, simulation.mesh.dx, simulation.mesh.dy
    )

    assert np.allclose(bz, 0.75, atol=1e-12)


def test_cross_section_current_integrates_current_density_and_thickness():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    simulation.fields.current_density_x.fill(2.0e6)
    expected = 2.0e6 * simulation.geometry.film.height * simulation.geometry.film.thickness

    assert DIAGNOSTICS["cross_section_current"](simulation) == pytest.approx(expected)
    assert DIAGNOSTICS["current_uniformity"](simulation) == pytest.approx(0.0)


def test_positive_x_sheet_current_produces_negative_y_field_above_film():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    jx = np.ones_like(simulation.fields.temperature)
    jy = np.zeros_like(jx)

    bx, by, bz = magnetic_field_from_sheet_current(
        simulation.mesh, jx, jy, simulation.geometry.film.thickness,
        observation_height=5e-9, source_stride=10,
    )

    center = (simulation.mesh.ny // 2, simulation.mesh.nx // 2)
    assert by[center] < 0.0
    assert abs(bx[center]) < abs(by[center])
    assert abs(bz[center]) < abs(by[center])


def test_induced_vector_potential_is_linear_and_parallel_to_sheet_current():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    jx = np.ones_like(simulation.fields.temperature) * 2.0e6
    jy = np.zeros_like(jx)
    ax, ay = induced_vector_potential_from_sheet_current(
        simulation.mesh, jx, jy, simulation.material_map.thickness,
        source_stride=1,
    )
    doubled_x, doubled_y = induced_vector_potential_from_sheet_current(
        simulation.mesh, 2.0 * jx, jy, simulation.material_map.thickness,
        source_stride=1,
    )
    assert np.all(ax > 0.0)
    assert np.allclose(ay, 0.0)
    assert np.allclose(doubled_x, 2.0 * ax)
    assert np.allclose(doubled_y, 0.0)


def test_screening_update_separates_applied_induced_and_total_potential():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    fields = simulation.fields
    applied_x, applied_y = uniform_perpendicular_vector_potential(
        simulation.mesh, 0.2, "symmetric"
    )
    fields.applied_vector_potential_x = applied_x
    fields.applied_vector_potential_y = applied_y
    fields.vector_potential_x = applied_x.copy()
    fields.vector_potential_y = applied_y.copy()
    fields.current_density_x.fill(1.0e7)
    fields.current_density_y.fill(0.0)
    simulation.config.electromagnetic.include_self_field = True
    simulation.config.electromagnetic.screening_step_size = 0.25

    magnetic_screening_step(simulation)

    assert np.any(fields.induced_vector_potential_x > 0.0)
    assert np.allclose(fields.induced_vector_potential_y, 0.0)
    assert np.allclose(
        fields.vector_potential_x,
        fields.applied_vector_potential_x + fields.induced_vector_potential_x,
    )
    assert np.allclose(
        fields.vector_potential_y,
        fields.applied_vector_potential_y + fields.induced_vector_potential_y,
    )
    assert fields.magnetic_screening_residual > 0.0


def test_current_controller_reaches_target_without_advancing_extra_time_levels():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    target = 1e-8

    result, trials, measured, _ = DIAGNOSTICS["current_controlled_step"](
        simulation, 1e-14, build_tdgl_model(simulation), target,
        relative_tolerance=1e-3, max_iterations=8, voltage_limit=0.01,
    )

    assert result.converged
    assert trials <= 8
    assert measured == pytest.approx(target, rel=1e-3)
    assert DIAGNOSTICS["cross_section_current"](simulation) == pytest.approx(
        target, rel=1e-3
    )
