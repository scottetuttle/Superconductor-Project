"""
SHS Electromagnetics Physics Module
===================================

Purpose:
--------
Models electromagnetic quantities associated with superconducting
devices.

This module describes the relationship between:

- Current density
- Electric potential
- Magnetic fields
- Vector potential
- Superconducting currents


Primary equations:
------------------

Current conservation:

        ∇ · J = 0


Electric field:

        E = -∇V - ∂A/∂t


where:

        V = electric potential
        A = magnetic vector potential


Magnetic response:

        ∇ × B = μ₀J


Superconducting coupling:
-------------------------

The electromagnetic module interacts with TDGL through the gauge
potential:

        (-i∇ - qA/ħ)ψ


This allows:

- Magnetic field effects
- Vortex formation
- Current distribution changes
- Self-field calculations


Planned capabilities:
---------------------

Phase 1:

- Solve electric potential
- Calculate current distribution
- Couple resistive regions


Phase 2:

- Magnetic field calculations
- Vector potential
- Self-field effects


Phase 3:

- Full Maxwell coupling
- Microwave response
- AC transport


Dependencies:
-------------

Requires:

- Geometry mesh
- Material electrical properties
- Superconductivity module


Current implementation status:
------------------------------

Architecture only.

Electromagnetic solver not implemented.
"""

from dataclasses import dataclass
import numpy as np
from scipy.signal import fftconvolve
from shs.utils.constants import VACUUM_PERMEABILITY
from shs.utils.defaults import default_section


_DEFAULTS = default_section("electromagnetic")


@dataclass
class ElectromagneticModel:
    """
    Global electromagnetic solver parameters.
    """

    reference_voltage: float = _DEFAULTS["reference_voltage"]

    ground_voltage: float = _DEFAULTS["ground_voltage"]


    magnetic_permeability: float = VACUUM_PERMEABILITY


    include_self_field: bool = _DEFAULTS["include_self_field"]

    include_displacement_current: bool = _DEFAULTS["include_displacement_current"]


def uniform_perpendicular_vector_potential(mesh, magnetic_field_z, gauge="symmetric"):
    """Return SI vector potential for a uniform perpendicular field.

    ``symmetric`` uses A=(-B(y-yc)/2, B(x-xc)/2). ``landau_x`` uses
    A=(-B(y-yc), 0), and ``landau_y`` uses A=(0, B(x-xc)). All have
    curl(A)_z = B and are useful for explicit gauge-invariance checks.
    """
    if not np.isfinite(magnetic_field_z):
        raise ValueError("Applied magnetic field must be finite.")
    x = mesh.x[None, :] - 0.5 * (mesh.x[0] + mesh.x[-1])
    y = mesh.y[:, None] - 0.5 * (mesh.y[0] + mesh.y[-1])
    shape = (mesh.ny, mesh.nx)
    if gauge == "symmetric":
        ax = np.broadcast_to(-0.5 * magnetic_field_z * y, shape).copy()
        ay = np.broadcast_to(0.5 * magnetic_field_z * x, shape).copy()
    elif gauge == "landau_x":
        ax = np.broadcast_to(-magnetic_field_z * y, shape).copy()
        ay = np.zeros(shape)
    elif gauge == "landau_y":
        ax = np.zeros(shape)
        ay = np.broadcast_to(magnetic_field_z * x, shape).copy()
    else:
        raise ValueError("Gauge must be symmetric, landau_x, or landau_y.")
    return ax, ay


def perpendicular_magnetic_field(vector_potential_x, vector_potential_y, dx, dy):
    """Calculate Bz=dAy/dx-dAx/dy at nodes from SI vector potential."""
    ax = np.asarray(vector_potential_x, dtype=float)
    ay = np.asarray(vector_potential_y, dtype=float)
    if ax.shape != ay.shape or ax.ndim != 2:
        raise ValueError("Vector-potential components must be matching 2D arrays.")
    if min(dx, dy) <= 0 or not np.isfinite(dx + dy):
        raise ValueError("Mesh spacing must be positive and finite.")
    return np.gradient(ay, dx, axis=1, edge_order=2) - np.gradient(
        ax, dy, axis=0, edge_order=2
    )


def magnetic_field_from_sheet_current(mesh, current_density_x, current_density_y,
                                      thickness, observation_height,
                                      source_stride=1):
    """Evaluate B from a thin-film current using the Biot-Savart law.

    The result is sampled on the mesh at a nonzero height above the film.
    This O(N^2) reference calculation is intended for diagnostics on small
    meshes, not for self-consistent magnetic screening iterations.
    """
    if thickness <= 0 or observation_height <= 0:
        raise ValueError("Thickness and observation height must be positive.")
    jx = np.asarray(current_density_x, dtype=float)
    jy = np.asarray(current_density_y, dtype=float)
    shape = (mesh.ny, mesh.nx)
    if jx.shape != shape or jy.shape != shape:
        raise ValueError("Current fields must match the mesh.")
    if not np.any(jx) and not np.any(jy):
        zeros = np.zeros(shape)
        return zeros.copy(), zeros.copy(), zeros.copy()
    if not isinstance(source_stride, int) or source_stride < 1:
        raise ValueError("source_stride must be a positive integer.")
    source_slice = np.s_[::source_stride, ::source_stride]
    kx = np.zeros(shape)
    ky = np.zeros(shape)
    kx[source_slice] = jx[source_slice] * thickness * source_stride**2
    ky[source_slice] = jy[source_slice] * thickness * source_stride**2
    # Convolution kernels depend only on observer-source displacement on this
    # uniform rectangular grid, so FFT convolution evaluates the same discrete
    # Biot-Savart sum without an O(N^2) Python loop.
    offset_x = np.arange(-(mesh.nx - 1), mesh.nx) * mesh.dx
    offset_y = np.arange(-(mesh.ny - 1), mesh.ny) * mesh.dy
    rx, ry = np.meshgrid(offset_x, offset_y)
    inverse_r3 = (rx**2 + ry**2 + observation_height**2) ** -1.5
    prefactor = VACUUM_PERMEABILITY * mesh.dx * mesh.dy / (4.0 * np.pi)
    bx = prefactor * fftconvolve(ky, observation_height * inverse_r3, mode="same")
    by = -prefactor * fftconvolve(kx, observation_height * inverse_r3, mode="same")
    bz = prefactor * (
        fftconvolve(kx, ry * inverse_r3, mode="same")
        - fftconvolve(ky, rx * inverse_r3, mode="same")
    )
    return bx, by, bz


def induced_vector_potential_from_sheet_current(
    mesh, current_density_x, current_density_y, thickness, source_stride=1,
):
    """Return the in-plane SI vector potential induced by a thin-film current.

    This discretizes ``A(r)=mu0/(4*pi) integral K(r')/|r-r'| dA'`` on the
    uniform grid. The singular same-cell contribution is replaced by the
    area-average for an equal-area circular cell, ``2/a`` with
    ``a=sqrt(dx*dy/pi)``.
    """
    jx = np.asarray(current_density_x, dtype=float)
    jy = np.asarray(current_density_y, dtype=float)
    shape = (mesh.ny, mesh.nx)
    if jx.shape != shape or jy.shape != shape:
        raise ValueError("Current fields must match the mesh.")
    thickness = np.broadcast_to(np.asarray(thickness, dtype=float), shape)
    if np.any(thickness <= 0) or not np.all(np.isfinite(thickness)):
        raise ValueError("Film thickness must be positive and finite.")
    if not isinstance(source_stride, int) or source_stride < 1:
        raise ValueError("source_stride must be a positive integer.")
    if not np.all(np.isfinite(jx)) or not np.all(np.isfinite(jy)):
        raise ValueError("Current fields must be finite.")
    if not np.any(jx) and not np.any(jy):
        zeros = np.zeros(shape)
        return zeros.copy(), zeros.copy()

    source_slice = np.s_[::source_stride, ::source_stride]
    kx = np.zeros(shape)
    ky = np.zeros(shape)
    scale = source_stride**2
    kx[source_slice] = jx[source_slice] * thickness[source_slice] * scale
    ky[source_slice] = jy[source_slice] * thickness[source_slice] * scale
    offset_x = np.arange(-(mesh.nx - 1), mesh.nx) * mesh.dx
    offset_y = np.arange(-(mesh.ny - 1), mesh.ny) * mesh.dy
    rx, ry = np.meshgrid(offset_x, offset_y)
    radius = np.hypot(rx, ry)
    kernel = np.zeros_like(radius)
    np.divide(1.0, radius, out=kernel, where=radius > 0)
    kernel[mesh.ny - 1, mesh.nx - 1] = 2.0 * np.sqrt(
        np.pi / (mesh.dx * mesh.dy)
    )
    prefactor = VACUUM_PERMEABILITY * mesh.dx * mesh.dy / (4.0 * np.pi)
    return (
        prefactor * fftconvolve(kx, kernel, mode="same"),
        prefactor * fftconvolve(ky, kernel, mode="same"),
    )
