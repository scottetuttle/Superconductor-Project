"""Physics-level diagnostics used to qualify adaptive solver decisions."""

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class PhysicsDiagnostics:
    current_continuity_rms: float
    current_continuity_max: float
    current_continuity_relative: float
    total_joule_power: float
    thermal_energy: float
    minimum_joule_density: float
    minimum_temperature: float
    maximum_order_parameter: float
    phonon_thermal_energy: float = 0.0
    electron_thermal_energy: float = 0.0

    def to_dict(self):
        return asdict(self)


def finite_volume_current_divergence(jx, jy, mesh, thickness):
    """Return charge divergence from face currents and nodal control volumes."""
    jx, jy = np.asarray(jx, dtype=float), np.asarray(jy, dtype=float)
    thickness = np.asarray(thickness, dtype=float)
    if jx.shape != jy.shape or jx.shape != thickness.shape:
        raise ValueError("Current and thickness arrays must have identical shapes.")
    wx, wy = np.ones(mesh.nx), np.ones(mesh.ny)
    wx[[0, -1]], wy[[0, -1]] = 0.5, 0.5
    tx = 0.5 * (thickness[:, :-1] + thickness[:, 1:])
    ty = 0.5 * (thickness[:-1] + thickness[1:])
    flux_x = jx[:, :-1] * tx * wy[:, None] * mesh.dy
    flux_y = jy[:-1, :] * ty * wx[None, :] * mesh.dx
    balance = np.zeros_like(jx, dtype=float)
    balance[:, :-1] += flux_x
    balance[:, 1:] -= flux_x
    balance[:-1] += flux_y
    balance[1:] -= flux_y
    volumes = thickness * wy[:, None] * wx[None, :] * mesh.dx * mesh.dy
    return balance / volumes


def _control_volume_weights(simulation):
    mesh = simulation.mesh
    wx = np.ones(mesh.nx)
    wy = np.ones(mesh.ny)
    wx[[0, -1]] = 0.5
    wy[[0, -1]] = 0.5
    return (
        wy[:, None]
        * wx[None, :]
        * mesh.dx
        * mesh.dy
        * simulation.material_map.thickness
    )


def evaluate_physics_diagnostics(simulation):
    """Measure conservation defects and basic physical bounds for one state."""
    fields = simulation.fields
    mesh = simulation.mesh
    divergence = finite_volume_current_divergence(
        fields.current_density_x,
        fields.current_density_y,
        mesh,
        simulation.material_map.thickness,
    )
    contact_mask = np.zeros_like(divergence, dtype=bool)
    for mask in simulation.contact_map.contact_masks.values():
        contact_mask |= mask
    free = ~contact_mask
    free_divergence = divergence[free]
    continuity_rms = float(
        np.sqrt(np.mean(free_divergence**2)) if free_divergence.size else 0.0
    )
    continuity_max = float(
        np.max(np.abs(free_divergence), initial=0.0)
    )
    current_scale = max(
        float(np.sqrt(np.mean(fields.current_density_x**2 + fields.current_density_y**2))),
        np.finfo(float).tiny,
    )
    continuity_relative = continuity_rms / (
        current_scale / min(mesh.dx, mesh.dy)
    )
    volumes = _control_volume_weights(simulation)
    thermal = getattr(getattr(simulation, "config", None), "thermal", None)
    two_temperature = getattr(thermal, "model", "single_temperature") == "two_temperature"
    phonon_temperature = fields.phonon_temperature if fields.phonon_temperature is not None else fields.temperature
    fraction = thermal.electron_heat_capacity_fraction if two_temperature else 1.0
    electron_energy = float(np.sum(fields.temperature * simulation.material_map.heat_capacity
                                   * fraction * volumes))
    phonon_energy = (float(np.sum(phonon_temperature * simulation.material_map.heat_capacity
                                 * (1-fraction) * volumes)) if two_temperature else 0.0)
    joule = (
        fields.joule_heat_source
        if fields.joule_heat_source is not None
        else np.zeros_like(fields.temperature)
    )
    return PhysicsDiagnostics(
        current_continuity_rms=continuity_rms,
        current_continuity_max=continuity_max,
        current_continuity_relative=float(continuity_relative),
        total_joule_power=float(np.sum(joule * volumes)),
        thermal_energy=electron_energy + phonon_energy,
        electron_thermal_energy=electron_energy,
        phonon_thermal_energy=phonon_energy,
        minimum_joule_density=float(np.min(joule)),
        minimum_temperature=float(np.min(fields.temperature)),
        maximum_order_parameter=float(np.max(np.abs(fields.psi))),
    )
