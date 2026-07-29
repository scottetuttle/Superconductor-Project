"""
Boundary condition definitions.

These classes describe the physical environment at the edge of the device.

They do NOT implement numerical algorithms.
Those belong inside the individual solvers.
"""

from dataclasses import dataclass, field
from enum import Enum


class BoundarySide(Enum):
    LEFT = "left"
    RIGHT = "right"
    TOP = "top"
    BOTTOM = "bottom"


class BoundaryType(Enum):
    """
    Physical boundary models.
    """

    FIXED_TEMPERATURE = "fixed_temperature"

    INSULATING = "insulating"

    CONVECTION = "convection"

    FIXED_HEAT_FLUX = "fixed_heat_flux"


@dataclass
class BoundaryCondition:
    """
    Physical description of one boundary.
    """

    side: BoundarySide

    model: BoundaryType

    # Temperature of surrounding bath (K)
    temperature: float | None = None

    # Heat transfer coefficient (W/m²K)
    heat_transfer: float | None = None

    # Applied heat flux (W/m²)
    heat_flux: float | None = None


@dataclass
class BoundarySet:

    boundaries: dict[BoundarySide, BoundaryCondition] = field(
        default_factory=dict
    )

    def add(
        self,
        boundary: BoundaryCondition,
    ):
        """
        Add or replace a boundary.
        """

        self.boundaries[boundary.side] = boundary

    def get(
        self,
        side: BoundarySide,
    ) -> BoundaryCondition:
        """
        Retrieve a boundary.
        """

        return self.boundaries[side]
    def __iter__(self):
        return iter(self.boundaries.values())