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

__all__ = [
    "Fields",
    "ThermalModel",
]