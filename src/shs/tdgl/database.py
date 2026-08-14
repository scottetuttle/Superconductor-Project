"""
TDGL boundary configuration loader.
"""

import json
from pathlib import Path

from .boundary import (
    TDGLBoundaryCondition,
    TDGLBoundarySet,
    TDGLBoundarySide,
    TDGLBoundaryType,
)


def load_tdgl_boundary_set(
    filepath: str | Path,
) -> TDGLBoundarySet:

    filepath = Path(filepath)

    with open(filepath, "r") as file:
        data = json.load(file)

    boundary_data = data["tdgl_boundaries"]

    boundary_set = TDGLBoundarySet()

    for side_name, values in boundary_data.items():

        side = TDGLBoundarySide(
            side_name
        )

        boundary_type = TDGLBoundaryType(
            values["type"]
        )

        boundary = TDGLBoundaryCondition(
            side=side,
            type=boundary_type,
        )

        boundary_set.add(boundary)

    return boundary_set