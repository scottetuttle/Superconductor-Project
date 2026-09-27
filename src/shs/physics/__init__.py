"""
SHS Physics Package

Contains the fundamental physical models:

- Thermal transport
- Superconductivity
- Electromagnetics
- Vortex dynamics

These modules define equations and physical relationships.
Numerical solution is handled separately in shs.solvers.
"""

from .fields import Fields
from .thermal import ThermalModel
from .electromagnetics import ElectromagneticModel
from .reduced_vortex import (
    ReducedPinningSite, ReducedVortexMaterial, ReducedVortexState,
    reduced_vortex_forces, reduced_vortex_step,
)

__all__ = [
    "Fields",
    "ThermalModel",
    "ReducedPinningSite",
    "ReducedVortexMaterial",
    "ReducedVortexState",
    "reduced_vortex_forces",
    "reduced_vortex_step",
]
