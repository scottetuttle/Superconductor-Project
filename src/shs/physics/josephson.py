"""Spatially resolved TDGL weak-link helpers."""

import numpy as np
from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM


def weak_link_suppression(simulation):
    """Return local GL-coefficient suppression for a vertical weak link."""
    config = getattr(getattr(simulation, "config", None), "josephson", None)
    shape = simulation.fields.psi.shape
    if config is None or not config.enabled or config.coefficient_suppression == 0:
        return np.zeros(shape)
    x = simulation.mesh.x
    mask_x = np.abs(x - config.center_x_m) <= 0.5 * config.width_m
    result = np.zeros(shape)
    result[:, mask_x] = config.coefficient_suppression
    result *= simulation.region_map.active_mask
    return result


def junction_observables(simulation):
    """Measure terminal voltage, bank amplitudes, and gauge-invariant phase drop."""
    config = simulation.config.josephson
    x = simulation.mesh.x
    left_index = max(0, int(np.searchsorted(x, config.center_x_m - config.width_m / 2)) - 1)
    right_index = min(simulation.mesh.nx - 1,
                      int(np.searchsorted(x, config.center_x_m + config.width_m / 2)))
    psi = simulation.fields.psi
    active = simulation.region_map.active_mask
    valid = active[:, left_index] & active[:, right_index]
    phase = np.angle(psi)
    ax = simulation.fields.vector_potential_x
    separation = (right_index - left_index) * simulation.mesh.dx
    delta = np.angle(np.exp(1j * (
        phase[:, right_index] - phase[:, left_index]
        - 0.5 * (ax[:, right_index] + ax[:, left_index]) * separation
        * (2.0 * np.pi / SUPERCONDUCTING_FLUX_QUANTUM)
    )))
    weights = np.abs(psi[:, left_index] * psi[:, right_index])
    selected = valid & (weights > 1e-12)
    phase_drop = float(np.angle(np.sum(weights[selected] * np.exp(1j * delta[selected])))) if np.any(selected) else np.nan
    source = simulation.contact_map.contact_masks.get(simulation.config.electrical.source_contact)
    sink = simulation.contact_map.contact_masks.get(simulation.config.electrical.sink_contact)
    voltage = float(np.mean(simulation.fields.voltage[source]) - np.mean(simulation.fields.voltage[sink]))
    link = np.abs(x - config.center_x_m) <= 0.5 * config.width_m
    return {
        "junction_voltage_V": voltage,
        "junction_phase_difference_rad": phase_drop,
        "junction_mean_amplitude": float(np.mean(np.abs(psi[:, link][active[:, link]]))),
        "left_bank_mean_amplitude": float(np.mean(np.abs(psi[:, left_index][valid]))),
        "right_bank_mean_amplitude": float(np.mean(np.abs(psi[:, right_index][valid]))),
    }
