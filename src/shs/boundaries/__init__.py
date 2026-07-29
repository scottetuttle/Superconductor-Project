"""
Boundary condition package.

Defines physical boundary conditions used by all SHS solvers.

Thermal
Electrical
TDGL
Electromagnetic

Each solver interprets the boundary according to its governing equations.
"""

from .boundary import (
    BoundaryCondition,
    BoundarySide,
    BoundaryType,
    BoundarySet,
)