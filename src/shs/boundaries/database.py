"""
Boundary configuration loader.

Converts JSON boundary definitions into BoundarySet objects.
"""

import json
from pathlib import Path

from .boundary import (
    BoundaryCondition,
    BoundarySet,
    BoundarySide,
    BoundaryType,
)


def load_boundary_set(filepath: str | Path) -> BoundarySet:
    """
    Load boundary conditions from a simulation JSON file.
    """

    filepath = Path(filepath)

    with open(filepath, "r") as file:
        data = json.load(file)

    boundary_data = data["boundaries"]

    boundary_set = BoundarySet()

    for side_name, values in boundary_data.items():

        side = BoundarySide(side_name)

        boundary_type = BoundaryType(
            values["type"]
        )

        boundary = BoundaryCondition(
            side=side,
            type=boundary_type,
            temperature=values.get(
                "temperature"
            ),
            heat_transfer=values.get(
                "heat_transfer"
            ),
            heat_flux=values.get(
                "heat_flux"
            ),
        )

        boundary_set.add(boundary)

    return boundary_set