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
    NORMAL_CONTACT = "normal_contact"
    NORMAL_CONTACT_MASK = "normal_contact_mask"


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
def apply_insulating_boundary(psi, vector_potential_x, vector_potential_y,
                              dx, dy, sides=None):
    """Zero normal link derivative on selected insulating edges.

    Link phases preserve amplitude under a pure gauge. At corners the y
    boundary is applied last; for nonzero plaquette flux, both inward link
    constraints cannot in general be imposed independently at one corner.
    """
    import numpy as np
    sides = set(TDGLBoundarySide if sides is None else sides)
    result = psi.copy()
    ux = np.exp(-1j*vector_potential_x*dx)
    uy = np.exp(-1j*vector_potential_y*dy)
    if TDGLBoundarySide.LEFT in sides:
        result[:, 0] = ux[:, 0]*result[:, 1]
    if TDGLBoundarySide.RIGHT in sides:
        result[:, -1] = np.conjugate(ux[:, -2])*result[:, -2]
    if TDGLBoundarySide.BOTTOM in sides:
        result[0] = uy[0]*result[1]
    if TDGLBoundarySide.TOP in sides:
        result[-1] = np.conjugate(uy[-2])*result[-2]
    return result


def apply_normal_contact_boundary(psi, sides):
    """Apply psi=0 on superconducting-to-normal terminal interfaces."""
    result = psi.copy()
    sides = set(sides)
    if TDGLBoundarySide.LEFT in sides:
        result[:, 0] = 0.0
    if TDGLBoundarySide.RIGHT in sides:
        result[:, -1] = 0.0
    if TDGLBoundarySide.BOTTOM in sides:
        result[0, :] = 0.0
    if TDGLBoundarySide.TOP in sides:
        result[-1, :] = 0.0
    return result


def apply_normal_contact_mask(psi, contact_mask):
    """Impose ``psi = 0`` only on nodes belonging to a normal terminal.

    The mask must describe nodes on the outer mesh boundary. Keeping this
    operation separate from the insulating edge update lets a single device
    edge contain both a normal contact segment and an insulating segment.
    """
    import numpy as np

    result = psi.copy()
    mask = np.asarray(contact_mask, dtype=bool)
    if mask.shape != result.shape:
        raise ValueError("Normal-contact mask must match the order-parameter grid.")
    boundary = np.zeros_like(mask)
    boundary[[0, -1], :] = True
    boundary[:, [0, -1]] = True
    if np.any(mask & ~boundary):
        raise ValueError("TDGL normal contacts must lie on the outer mesh boundary.")
    result[mask] = 0.0
    return result
