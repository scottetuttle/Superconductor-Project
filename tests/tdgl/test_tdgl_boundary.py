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
        (result[:, 1] - result[:, 0]) / dx
        - 1j * Ax[:, 0] * result[:, 0]
    )

    assert np.allclose(
        Dx_left,
        0.0,
    )

    # Right

    Dx_right = (
        (result[:, -1] - result[:, -2]) / dx
        - 1j * Ax[:, -1] * result[:, -1]
    )

    assert np.allclose(
        Dx_right,
        0.0,
    )

    # Bottom

    Dy_bottom = (
        (result[1, :] - result[0, :]) / dy
        - 1j * Ay[0, :] * result[0, :]
    )

    assert np.allclose(
        Dy_bottom,
        0.0,
    )

    # Top

    Dy_top = (
        (result[-1, :] - result[-2, :]) / dy
        - 1j * Ay[-1, :] * result[-1, :]
    )

    assert np.allclose(
        Dy_top,
        0.0,
    )