"""
Geometry configuration loader.

Loads geometry definitions from JSON files.
"""

import json
from pathlib import Path

from .geometry import RectangularFilm


def load_geometry(filepath: str | Path):
    """
    Load geometry from JSON configuration.

    Parameters
    ----------
    filepath:
        Path to geometry JSON file.

    Returns
    -------
    Geometry object
    """

    filepath = Path(filepath)

    with open(filepath, "r") as file:
        data = json.load(file)

    geometry_type = data.pop("type")


    if geometry_type == "rectangle":

        return RectangularFilm(**data)


    raise ValueError(
        f"Unknown geometry type: {geometry_type}"
    )