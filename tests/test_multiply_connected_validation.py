"""Reference checks for multiply connected superconducting domains."""

import numpy as np
import pytest

from shs.config.builder import build_simulation
from shs.physics.electromagnetics import uniform_perpendicular_vector_potential
from shs.physics.fluxoid import contour_around_hole, hole_fluxoid, rectangular_fluxoid
from shs.solvers.coupled_solver import build_tdgl_model
from shs.solvers.tdgl_solver import tdgl_scales, tdgl_step
from shs.tdgl.operators import covariant_laplacian


@pytest.fixture
def ring():
    return build_simulation("configs/simulations/nbn_ring_validation.json")


def test_circular_hole_area_converges_to_analytic_area(ring):
    inactive_area = np.count_nonzero(~ring.region_map.active_mask) * ring.mesh.dx * ring.mesh.dy
    expected = np.pi * ring.geometry.holes[0].radius**2
    assert inactive_area == pytest.approx(expected, rel=0.03)


def test_masked_constant_state_has_zero_neumann_laplacian(ring):
    active = ring.region_map.active_mask
    psi = np.ones(active.shape, dtype=complex)
    zeros = np.zeros(active.shape)
    result = covariant_laplacian(psi, zeros, zeros, 0.2, 0.2, active)
    assert np.max(np.abs(result[active])) == pytest.approx(0.0, abs=1e-14)


def test_inactive_values_cannot_couple_across_hole(ring):
    active = ring.region_map.active_mask
    zeros = np.zeros(active.shape)
    reference = np.ones(active.shape, dtype=complex)
    perturbed = reference.copy()
    perturbed[~active] = 1e6 * (1.0 + 1j)
    first = covariant_laplacian(reference, zeros, zeros, 0.2, 0.2, active)
    second = covariant_laplacian(perturbed, zeros, zeros, 0.2, 0.2, active)
    assert np.array_equal(first[active], second[active])


def test_masked_operator_is_gauge_covariant(ring):
    active = ring.region_map.active_mask
    yy, xx = np.indices(active.shape)
    chi = 0.013 * xx + 0.009 * yy
    psi = np.ones(active.shape, dtype=complex)
    zeros = np.zeros(active.shape)
    original = covariant_laplacian(psi, zeros, zeros, 1.0, 1.0, active)
    transformed = covariant_laplacian(
        psi * np.exp(1j * chi),
        np.full(active.shape, 0.013),
        np.full(active.shape, 0.009),
        1.0,
        1.0,
        active,
    )
    assert np.allclose(transformed[active], np.exp(1j * chi[active]) * original[active], atol=1e-12)


def test_integer_ring_winding_matches_one_flux_quantum(ring):
    yy, xx = np.indices(ring.fields.psi.shape)
    phase = np.arctan2(yy - 50, xx - 50)
    ring.fields.psi = np.exp(1j * phase)
    ring.fields.psi[~ring.region_map.active_mask] = 0.0
    model = build_tdgl_model(ring)
    scales = tdgl_scales(ring, model)
    jx, jy = model.supercurrent_density(
        ring.fields.psi,
        np.zeros_like(phase),
        np.zeros_like(phase),
        ring.mesh.dx / scales.xi,
        ring.mesh.dy / scales.xi,
        ring.region_map.active_mask,
    )
    ring.fields.supercurrent_density_x = scales.supercurrent_density_to_physical(jx)
    ring.fields.supercurrent_density_y = scales.supercurrent_density_to_physical(jy)
    result = rectangular_fluxoid(ring, (20, 80, 20, 80))
    assert result.winding == 1
    assert result.london_fluxoid_quanta == pytest.approx(1.0, rel=3e-4)


def test_named_hole_contour_uniform_field_flux_matches_b_times_area(ring):
    field = 0.037
    ax, ay = uniform_perpendicular_vector_potential(ring.mesh, field, "symmetric")
    ring.fields.vector_potential_x = ax
    ring.fields.vector_potential_y = ay
    bounds = contour_around_hole(ring, "central_hole", 3e-9)
    y0, y1, x0, x1 = bounds
    result = hole_fluxoid(ring, "central_hole", 3e-9)
    expected_flux = field * (x1 - x0) * ring.mesh.dx * (y1 - y0) * ring.mesh.dy
    assert result.flux_part_Wb == pytest.approx(expected_flux, rel=2e-14)


def test_tdgl_evolution_keeps_hole_inactive(ring):
    tdgl_step(ring, ring.config.dt, build_tdgl_model(ring))
    assert np.all(ring.fields.psi[~ring.region_map.active_mask] == 0.0)
    assert np.all(ring.fields.supercurrent_density_x[~ring.region_map.active_mask] == 0.0)
    assert np.all(ring.fields.supercurrent_density_y[~ring.region_map.active_mask] == 0.0)
