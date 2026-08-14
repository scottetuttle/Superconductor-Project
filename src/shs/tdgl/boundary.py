"""
TDGL boundary condition definitions.

These classes describe physical boundary conditions
for the superconducting order parameter.
"""

from dataclasses import dataclass, field
from enum import Enum


class TDGLBoundarySide(Enum):
    LEFT = "left"
    RIGHT = "right"
    TOP = "top"
    BOTTOM = "bottom"


class TDGLBoundaryType(Enum):
    """
    Physical TDGL boundary models.
    """

    INSULATING = "insulating"


@dataclass
class TDGLBoundaryCondition:
    """
    Physical boundary condition for the TDGL order parameter.
    """

    side: TDGLBoundarySide

    type: TDGLBoundaryType


@dataclass
class TDGLBoundarySet:
    """
    Collection of TDGL boundary conditions.
    """

    boundaries: dict[
        TDGLBoundarySide,
        TDGLBoundaryCondition
    ] = field(default_factory=dict)

    def add(
        self,
        boundary: TDGLBoundaryCondition,
    ) -> None:

        self.boundaries[boundary.side] = boundary

    def get(
        self,
        side: TDGLBoundarySide,
    ) -> TDGLBoundaryCondition:

        return self.boundaries[side]

    def __contains__(
        self,
        side: TDGLBoundarySide,
    ) -> bool:

        return side in self.boundaries

    def __iter__(self):

        return iter(
            self.boundaries.values()
        )
def apply_insulating_boundary(
    psi,
    vector_potential_x,
    vector_potential_y,
    dx,
    dy,
):
    """
    Enforce the insulating TDGL boundary condition:

        n · D psi = 0

    with

        D = nabla - i A.
    """

    psi = psi.copy()

    # Left boundary
    #
    # (psi[:, 1] - psi[:, 0]) / dx
    #     - i Ax psi[:, 0] = 0

    Ax = vector_potential_x[:, 0]

    psi[:, 0] = (
        psi[:, 1]
        /
        (1.0 + 1j * Ax * dx)
    )

    # Right boundary
    #
    # (psi[:, -1] - psi[:, -2]) / dx
    #     - i Ax psi[:, -1] = 0

    Ax = vector_potential_x[:, -1]

    psi[:, -1] = (
        psi[:, -2]
        /
        (1.0 - 1j * Ax * dx)
    )

    # Bottom boundary
    #
    # (psi[1, :] - psi[0, :]) / dy
    #     - i Ay psi[0, :] = 0

    Ay = vector_potential_y[0, :]

    psi[0, :] = (
        psi[1, :]
        /
        (1.0 + 1j * Ay * dy)
    )

    # Top boundary
    #
    # (psi[-1, :] - psi[-2, :]) / dy
    #     - i Ay psi[-1, :] = 0

    Ay = vector_potential_y[-1, :]

    psi[-1, :] = (
        psi[-2, :]
        /
        (1.0 - 1j * Ay * dy)
    )

    return psi