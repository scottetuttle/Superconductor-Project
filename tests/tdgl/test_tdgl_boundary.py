from shs.tdgl import(
    TDGLBoundaryCondition,
    TDGLBoundarySide,
    TDGLBoundaryType,
    TDGLBoundarySet
)


def test_tdgl_boundary_creation():

    boundary = TDGLBoundaryCondition(
        side=TDGLBoundarySide.LEFT,
        type=TDGLBoundaryType.INSULATING,
    )

    assert boundary.side == TDGLBoundarySide.LEFT

    assert boundary.type == TDGLBoundaryType.INSULATING


def test_tdgl_boundary_set():

    boundaries = TDGLBoundarySet()

    boundaries.add(
        TDGLBoundaryCondition(
            side=TDGLBoundarySide.LEFT,
            type=TDGLBoundaryType.INSULATING,
        )
    )

    assert TDGLBoundarySide.LEFT in boundaries

import numpy as np

from shs.tdgl.boundary import (
    apply_insulating_boundary,
    apply_normal_contact_boundary,
    apply_normal_contact_mask,
)


def test_insulating_boundary_zero_vector_potential():

    psi = np.ones(
        (20, 20),
        dtype=complex,
    )

    Ax = np.zeros(
        (20, 20),
        dtype=float,
    )

    Ay = np.zeros(
        (20, 20),
        dtype=float,
    )

    result = apply_insulating_boundary(
        psi,
        Ax,
        Ay,
        1.0,
        1.0,
    )

    assert np.allclose(
        result,
        psi,
    )


def test_insulating_boundary_preserves_shape():

    psi = np.random.random(
        (20, 20)
    ) + 1j * np.random.random(
        (20, 20)
    )

    Ax = np.zeros(
        (20, 20)
    )

    Ay = np.zeros(
        (20, 20)
    )

    result = apply_insulating_boundary(
        psi,
        Ax,
        Ay,
        1.0,
        1.0,
    )

    assert result.shape == psi.shape


def test_insulating_boundary_enforces_covariant_normal_derivative():

    psi = np.random.random(
        (20, 20)
    ) + 1j * np.random.random(
        (20, 20)
    )

    Ax = np.full(
        (20, 20),
        0.1,
    )

    Ay = np.full(
        (20, 20),
        0.2,
    )

    dx = 1.0
    dy = 1.0

    result = apply_insulating_boundary(
        psi,
        Ax,
        Ay,
        dx,
        dy,
    )

    # Left

        # Left

    Dx_left = (
        (np.exp(-1j*Ax[:, 0]*dx)*result[:, 1] - result[:, 0]) / dx
    )

    assert np.allclose(
        Dx_left,
        0.0,
    )

    # Right

    Dx_right = (
        (np.exp(-1j*Ax[:, -2]*dx)*result[:, -1] - result[:, -2]) / dx
    )

    assert np.allclose(
        Dx_right,
        0.0,
    )

    # Bottom

    Dy_bottom = (
        (np.exp(-1j*Ay[0, :]*dy)*result[1, :] - result[0, :]) / dy
    )

    assert np.allclose(
        Dy_bottom,
        0.0,
    )

    # Top

    Dy_top = (
        (np.exp(-1j*Ay[-2, :]*dy)*result[-1, :] - result[-2, :]) / dy
    )

    assert np.allclose(
        Dy_top,
        0.0,
    )


def test_normal_contact_boundary_sets_only_selected_edges_to_zero():
    psi = np.ones((5, 6), dtype=complex)

    result = apply_normal_contact_boundary(
        psi,
        [TDGLBoundarySide.LEFT, TDGLBoundarySide.TOP],
    )

    assert np.all(result[:, 0] == 0.0)
    assert np.all(result[-1, :] == 0.0)
    assert np.all(result[:-1, 1:] == 1.0)


def test_normal_contact_mask_supports_partial_edge_terminal():
    psi = np.ones((7, 9), dtype=complex)
    mask = np.zeros_like(psi, dtype=bool)
    mask[2:5, 0] = True

    result = apply_normal_contact_mask(psi, mask)

    assert np.all(result[2:5, 0] == 0.0)
    assert np.all(result[[0, 1, 5, 6], 0] == 1.0)
    assert np.all(result[:, 1:] == 1.0)


def test_normal_contact_mask_rejects_interior_terminal():
    psi = np.ones((7, 9), dtype=complex)
    mask = np.zeros_like(psi, dtype=bool)
    mask[3, 4] = True

    import pytest
    with pytest.raises(ValueError, match="outer mesh boundary"):
        apply_normal_contact_mask(psi, mask)
