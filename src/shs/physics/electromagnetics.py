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


@dataclass
class ElectromagneticModel:
    """
    Global electromagnetic solver parameters.
    """

    reference_voltage: float = 1.0

    ground_voltage: float = 0.0


    magnetic_permeability: float = (
        4.0 *
        np.pi *
        1e-7
    )


    include_self_field: bool = False

    include_displacement_current: bool = False