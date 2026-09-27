"""Gauge-invariant rectangular-contour fluxoid diagnostics."""

from dataclasses import dataclass

import numpy as np

from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM, VACUUM_PERMEABILITY


@dataclass(frozen=True)
class FluxoidResult:
    flux_part_Wb: float
    supercurrent_part_Wb: float
    london_fluxoid_Wb: float
    topological_fluxoid_Wb: float
    winding: int

    @property
    def london_fluxoid_quanta(self):
        return self.london_fluxoid_Wb / SUPERCONDUCTING_FLUX_QUANTUM

    @property
    def topological_fluxoid_quanta(self):
        return self.topological_fluxoid_Wb / SUPERCONDUCTING_FLUX_QUANTUM


def contour_around_hole(simulation, hole_name, clearance):
    """Return grid bounds for a superconducting contour around a named hole.

    ``clearance`` is the physical distance between the hole's bounding box and
    the requested contour. The selected nodes are the nearest nodes outside
    that box, and the completed contour is checked by ``rectangular_fluxoid``.
    """
    if not np.isfinite(clearance) or clearance < 0:
        raise ValueError("Fluxoid contour clearance must be finite and nonnegative.")
    matches = [hole for hole in simulation.geometry.holes or [] if hole.name == hole_name]
    if len(matches) != 1:
        raise ValueError(f"Expected one hole named {hole_name!r}.")
    hole = matches[0]
    if hole.shape == "circle":
        half_width = half_height = hole.radius
    else:
        half_width, half_height = 0.5 * hole.width, 0.5 * hole.height
    x0 = int(np.searchsorted(simulation.mesh.x, hole.x - half_width - clearance, side="right") - 1)
    x1 = int(np.searchsorted(simulation.mesh.x, hole.x + half_width + clearance, side="left"))
    y0 = int(np.searchsorted(simulation.mesh.y, hole.y - half_height - clearance, side="right") - 1)
    y1 = int(np.searchsorted(simulation.mesh.y, hole.y + half_height + clearance, side="left"))
    bounds = (y0, y1, x0, x1)
    if not (0 <= y0 < y1 < simulation.mesh.ny and 0 <= x0 < x1 < simulation.mesh.nx):
        raise ValueError("Requested hole contour does not fit inside the mesh.")
    return bounds


def hole_fluxoid(simulation, hole_name, clearance, amplitude_floor=1e-8):
    """Evaluate a fluxoid on a generated contour around a named hole."""
    return rectangular_fluxoid(
        simulation,
        contour_around_hole(simulation, hole_name, clearance),
        amplitude_floor,
    )


def _closed_phase_winding(phase, y0, y1, x0, x1):
    wrapped = lambda value: np.angle(np.exp(1j * value))
    total = 0.0
    total += np.sum(wrapped(phase[y0, x0 + 1:x1 + 1] - phase[y0, x0:x1]))
    total += np.sum(wrapped(phase[y0 + 1:y1 + 1, x1] - phase[y0:y1, x1]))
    total -= np.sum(wrapped(phase[y1, x0 + 1:x1 + 1] - phase[y1, x0:x1]))
    total -= np.sum(wrapped(phase[y0 + 1:y1 + 1, x0] - phase[y0:y1, x0]))
    return int(np.rint(total / (2.0 * np.pi)))


def rectangular_fluxoid(simulation, bounds, amplitude_floor=1e-8):
    """Evaluate London and phase-winding fluxoids on a grid-aligned contour.

    ``bounds`` is ``(y0, y1, x0, x1)`` with inclusive corner nodes and a
    counterclockwise contour. The entire contour must remain superconducting.
    """
    y0, y1, x0, x1 = map(int, bounds)
    fields, mesh, material = simulation.fields, simulation.mesh, simulation.material_map
    if not (0 <= y0 < y1 < mesh.ny and 0 <= x0 < x1 < mesh.nx):
        raise ValueError("Fluxoid bounds must define a nonempty contour inside the mesh.")
    psi = fields.psi
    contour_amplitude = np.concatenate((
        np.abs(psi[y0, x0:x1 + 1]), np.abs(psi[y1, x0:x1 + 1]),
        np.abs(psi[y0 + 1:y1, x0]), np.abs(psi[y0 + 1:y1, x1]),
    ))
    active = simulation.region_map.active_mask
    contour_active = np.concatenate((
        active[y0, x0:x1 + 1], active[y1, x0:x1 + 1],
        active[y0 + 1:y1, x0], active[y0 + 1:y1, x1],
    ))
    if not np.all(contour_active):
        raise ValueError("Fluxoid contour crosses an inactive part of the device.")
    if np.any(contour_amplitude <= amplitude_floor):
        raise ValueError("Fluxoid contour crosses a point where phase is undefined.")

    ax, ay = fields.vector_potential_x, fields.vector_potential_y
    flux_part = (
        np.sum(ax[y0, x0:x1]) * mesh.dx
        + np.sum(ay[y0:y1, x1]) * mesh.dy
        - np.sum(ax[y1, x0:x1]) * mesh.dx
        - np.sum(ay[y0:y1, x0]) * mesh.dy
    )

    ns_x_bottom = 0.5 * (np.abs(psi[y0, x0:x1])**2 + np.abs(psi[y0, x0 + 1:x1 + 1])**2)
    ns_x_top = 0.5 * (np.abs(psi[y1, x0:x1])**2 + np.abs(psi[y1, x0 + 1:x1 + 1])**2)
    ns_y_left = 0.5 * (np.abs(psi[y0:y1, x0])**2 + np.abs(psi[y0 + 1:y1 + 1, x0])**2)
    ns_y_right = 0.5 * (np.abs(psi[y0:y1, x1])**2 + np.abs(psi[y0 + 1:y1 + 1, x1])**2)
    lambda2 = material.penetration_depth**2
    lambda_x_bottom = 0.5 * (lambda2[y0, x0:x1] + lambda2[y0, x0 + 1:x1 + 1])
    lambda_x_top = 0.5 * (lambda2[y1, x0:x1] + lambda2[y1, x0 + 1:x1 + 1])
    lambda_y_left = 0.5 * (lambda2[y0:y1, x0] + lambda2[y0 + 1:y1 + 1, x0])
    lambda_y_right = 0.5 * (lambda2[y0:y1, x1] + lambda2[y0 + 1:y1 + 1, x1])
    js_x, js_y = fields.supercurrent_density_x, fields.supercurrent_density_y
    circulation = (
        np.sum(lambda_x_bottom * js_x[y0, x0:x1] / ns_x_bottom) * mesh.dx
        + np.sum(lambda_y_right * js_y[y0:y1, x1] / ns_y_right) * mesh.dy
        - np.sum(lambda_x_top * js_x[y1, x0:x1] / ns_x_top) * mesh.dx
        - np.sum(lambda_y_left * js_y[y0:y1, x0] / ns_y_left) * mesh.dy
    )
    supercurrent_part = VACUUM_PERMEABILITY * circulation
    winding = _closed_phase_winding(np.angle(psi), y0, y1, x0, x1)
    return FluxoidResult(
        flux_part_Wb=float(flux_part),
        supercurrent_part_Wb=float(supercurrent_part),
        london_fluxoid_Wb=float(flux_part + supercurrent_part),
        topological_fluxoid_Wb=float(winding * SUPERCONDUCTING_FLUX_QUANTUM),
        winding=winding,
    )
