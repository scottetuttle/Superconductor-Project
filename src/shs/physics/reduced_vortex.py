"""Overdamped reduced dynamics for optically manipulated Abrikosov vortices.

This model is a screening tool for parameter searches. Forces and drag are per
unit vortex length, so their ratio gives the in-plane vortex velocity.
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import k1

from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM, VACUUM_PERMEABILITY


@dataclass(frozen=True)
class ReducedVortexMaterial:
    critical_temperature_K: float
    coherence_length_m: float
    penetration_depth_m: float
    normal_resistivity_ohm_m: float

    def validate(self):
        values = (
            self.critical_temperature_K, self.coherence_length_m,
            self.penetration_depth_m, self.normal_resistivity_ohm_m,
        )
        if not all(np.isfinite(value) and value > 0 for value in values):
            raise ValueError("Reduced-vortex material values must be positive and finite.")
        if self.penetration_depth_m <= self.coherence_length_m:
            raise ValueError("Reduced model requires a type-II material with lambda > xi.")

    @property
    def viscosity_N_s_per_m2(self):
        """Bardeen-Stephen viscous drag coefficient per vortex length."""
        return SUPERCONDUCTING_FLUX_QUANTUM**2 / (
            2 * np.pi * self.coherence_length_m**2 * self.normal_resistivity_ohm_m
        )

    @property
    def thermal_force_coefficient_N_per_K(self):
        """Magnitude of -d(vortex energy/length)/dT in the GL estimate."""
        return SUPERCONDUCTING_FLUX_QUANTUM**2 * np.log(
            self.penetration_depth_m / self.coherence_length_m
        ) / (
            4 * np.pi * VACUUM_PERMEABILITY
            * self.penetration_depth_m**2 * self.critical_temperature_K
        )


@dataclass(frozen=True)
class ReducedPinningSite:
    x_m: float
    y_m: float
    energy_per_length_J_per_m: float
    sigma_m: float


@dataclass
class ReducedVortexState:
    positions_m: np.ndarray
    charges: np.ndarray
    thermal_center_m: np.ndarray
    thermal_delta_K: float = 0.0

    def __post_init__(self):
        self.positions_m = np.asarray(self.positions_m, dtype=float)
        self.charges = np.asarray(self.charges, dtype=int)
        self.thermal_center_m = np.asarray(self.thermal_center_m, dtype=float)
        if self.positions_m.ndim != 2 or self.positions_m.shape[1] != 2:
            raise ValueError("Vortex positions must have shape (count, 2).")
        if self.charges.shape != (len(self.positions_m),):
            raise ValueError("One charge is required for every vortex.")
        if np.any(np.abs(self.charges) != 1):
            raise ValueError("Reduced vortices must have charge +1 or -1.")


def gaussian_temperature(position_m, center_m, bath_K, delta_K, sigma_m):
    displacement = np.asarray(position_m) - np.asarray(center_m)
    exponent = -np.sum(displacement**2, axis=-1) / (2 * sigma_m**2)
    return bath_K + delta_K * np.exp(exponent)


def reduced_vortex_forces(
    state, material, *, thermal_sigma_m, pinning_sites=(),
    current_density_A_per_m2=(0.0, 0.0), minimum_separation_m=None,
):
    """Return total and component forces per unit vortex length [N/m]."""
    material.validate()
    positions, charges = state.positions_m, state.charges
    count = len(positions)
    thermal = np.zeros_like(positions)
    pinning = np.zeros_like(positions)
    interaction = np.zeros_like(positions)
    lorentz = np.zeros_like(positions)

    offset = positions - state.thermal_center_m
    profile = state.thermal_delta_K * np.exp(
        -np.sum(offset**2, axis=1) / (2 * thermal_sigma_m**2)
    )
    gradient_temperature = -offset * profile[:, None] / thermal_sigma_m**2
    thermal[:] = material.thermal_force_coefficient_N_per_K * gradient_temperature

    for site in pinning_sites:
        offset = positions - (site.x_m, site.y_m)
        profile = np.exp(-np.sum(offset**2, axis=1) / (2 * site.sigma_m**2))
        pinning -= (
            site.energy_per_length_J_per_m * profile[:, None]
            * offset / site.sigma_m**2
        )

    cutoff = minimum_separation_m or 0.25 * material.coherence_length_m
    coefficient = SUPERCONDUCTING_FLUX_QUANTUM**2 / (
        2 * np.pi * VACUUM_PERMEABILITY * material.penetration_depth_m**3
    )
    for first in range(count):
        for second in range(first + 1, count):
            offset = positions[first] - positions[second]
            distance = max(float(np.linalg.norm(offset)), cutoff)
            direction = offset / distance
            pair_force = (
                charges[first] * charges[second] * coefficient
                * k1(distance / material.penetration_depth_m) * direction
            )
            interaction[first] += pair_force
            interaction[second] -= pair_force

    jx, jy = map(float, current_density_A_per_m2)
    lorentz[:, 0] = charges * SUPERCONDUCTING_FLUX_QUANTUM * jy
    lorentz[:, 1] = -charges * SUPERCONDUCTING_FLUX_QUANTUM * jx
    components = {
        "thermal": thermal, "pinning": pinning,
        "interaction": interaction, "lorentz": lorentz,
    }
    return sum(components.values(), np.zeros_like(positions)), components


def reduced_vortex_step(
    state, material, dt_s, *, commanded_center_m, commanded_delta_K,
    thermal_sigma_m, thermal_response_time_s, pinning_sites=(),
    current_density_A_per_m2=(0.0, 0.0), bounds_m=None,
    maximum_displacement_m=None,
):
    """Advance the thermal lag and overdamped vortex positions by one step."""
    if not np.isfinite(dt_s) or dt_s <= 0:
        raise ValueError("Reduced-vortex timestep must be positive and finite.")
    if thermal_sigma_m <= 0 or thermal_response_time_s <= 0:
        raise ValueError("Thermal width and response time must be positive.")
    relaxation = 1.0 - np.exp(-dt_s / thermal_response_time_s)
    state.thermal_center_m += relaxation * (
        np.asarray(commanded_center_m, dtype=float) - state.thermal_center_m
    )
    state.thermal_delta_K += relaxation * (
        float(commanded_delta_K) - state.thermal_delta_K
    )
    force, components = reduced_vortex_forces(
        state, material, thermal_sigma_m=thermal_sigma_m,
        pinning_sites=pinning_sites,
        current_density_A_per_m2=current_density_A_per_m2,
    )
    displacement = dt_s * force / material.viscosity_N_s_per_m2
    if maximum_displacement_m is not None:
        lengths = np.linalg.norm(displacement, axis=1)
        factors = np.minimum(1.0, maximum_displacement_m / np.maximum(lengths, 1e-300))
        displacement *= factors[:, None]
    state.positions_m += displacement
    if bounds_m is not None:
        xmin, xmax, ymin, ymax = map(float, bounds_m)
        state.positions_m[:, 0] = np.clip(state.positions_m[:, 0], xmin, xmax)
        state.positions_m[:, 1] = np.clip(state.positions_m[:, 1], ymin, ymax)
    return force, components
