"""Lumped SIS Josephson junction dynamics using the RCSJ model."""
from dataclasses import dataclass
import numpy as np
from shs.utils.constants import (SUPERCONDUCTING_FLUX_QUANTUM,
                                 BOLTZMANN_CONSTANT)

ELEMENTARY_CHARGE = 1.602176634e-19


def bcs_gap(T, Tc):
    """Approximate weak-coupling BCS energy gap in joules."""
    T = np.asarray(T, dtype=float)
    if Tc <= 0 or np.any(T < 0): raise ValueError("Temperatures must be valid.")
    reduced = np.clip(T / Tc, 0.0, 1.0)
    argument = np.zeros_like(reduced)
    positive = reduced > 0
    argument[positive] = 1.74 * np.sqrt(np.maximum(1.0 / reduced[positive] - 1.0, 0.0))
    argument[~positive] = np.inf
    return 1.764 * BOLTZMANN_CONSTANT * Tc * np.tanh(argument)


def ambegaokar_baratoff_critical_current(T, Tc, normal_resistance):
    if normal_resistance <= 0: raise ValueError("Normal resistance must be positive.")
    gap = bcs_gap(T, Tc)
    thermal = 1.0 if T == 0 else np.tanh(gap / (2 * BOLTZMANN_CONSTANT * T))
    return float(np.pi * gap * thermal / (2 * ELEMENTARY_CHARGE * normal_resistance))


@dataclass(frozen=True)
class RCSJJunction:
    critical_current_A: float
    normal_resistance_ohm: float
    capacitance_F: float
    magnetic_flux_Wb: float = 0.0
    interference_amplitude: float = 1.0
    phase_offset_rad: float = 0.0

    @property
    def effective_critical_current_A(self):
        x = np.pi * self.magnetic_flux_Wb / SUPERCONDUCTING_FLUX_QUANTUM
        return self.critical_current_A * self.interference_amplitude * (1.0 if x == 0 else np.sin(x) / x)

    @property
    def josephson_energy_J(self):
        return SUPERCONDUCTING_FLUX_QUANTUM * abs(self.effective_critical_current_A) / (2*np.pi)

    def validate(self):
        if self.critical_current_A < 0 or self.normal_resistance_ohm <= 0 or self.capacitance_F < 0 or not 0 <= self.interference_amplitude <= 1:
            raise ValueError("RCSJ parameters require Ic>=0, R>0, and C>=0.")


def extended_junction_interference(width_m, junction_x_m, magnetic_field_T,
                                   effective_magnetic_length_m, vortices=(), samples=401):
    """Return complex aperture average for field and electrode-vortex phase.

    Vortex entries contain x_m, y_m, charge, and electrode ('left'/'right').
    This is a London phase-texture approximation; vortex cores are not resolved.
    """
    if width_m <= 0 or effective_magnetic_length_m < 0 or samples < 3:
        raise ValueError("Extended-junction dimensions and sampling are invalid.")
    y=np.linspace(-width_m/2,width_m/2,int(samples))
    theta=2*np.pi*magnetic_field_T*effective_magnetic_length_m*y/SUPERCONDUCTING_FLUX_QUANTUM
    for vortex in vortices:
        sign=1.0 if vortex.get("electrode","right")=="right" else -1.0
        theta += sign*int(vortex["charge"])*np.arctan2(y-float(vortex["y_m"]),
                                                        junction_x_m-float(vortex["x_m"]))
    factor=np.trapezoid(np.exp(1j*theta),y)/width_m
    return complex(factor), y, theta


def current_biased_rhs(junction, phase, phase_rate, bias_current_A):
    """Return dphase/dt and d2phase/dt2 for a current-biased RCSJ."""
    junction.validate(); scale = SUPERCONDUCTING_FLUX_QUANTUM/(2*np.pi)
    if junction.capacitance_F == 0:
        rate = (bias_current_A-junction.effective_critical_current_A*np.sin(phase+junction.phase_offset_rad))*junction.normal_resistance_ohm/scale
        return rate, 0.0
    acceleration = (bias_current_A-junction.effective_critical_current_A*np.sin(phase+junction.phase_offset_rad)
                    -scale*phase_rate/junction.normal_resistance_ohm)/(junction.capacitance_F*scale)
    return phase_rate, acceleration


def junction_components(junction, phase, phase_rate, phase_acceleration=0.0):
    scale=SUPERCONDUCTING_FLUX_QUANTUM/(2*np.pi); voltage=scale*phase_rate
    return {"voltage_V": voltage,
            "supercurrent_A": junction.effective_critical_current_A*np.sin(phase+junction.phase_offset_rad),
            "resistive_current_A": voltage/junction.normal_resistance_ohm,
            "capacitive_current_A": junction.capacitance_F*scale*phase_acceleration,
            "phase_rad": phase, "phase_rate_rad_per_s": phase_rate,
            "joule_power_W": voltage**2/junction.normal_resistance_ohm}
