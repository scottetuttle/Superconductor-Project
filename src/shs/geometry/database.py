"""
Geometry configuration loader.

Loads geometry definitions from JSON files.
"""

import json
from pathlib import Path
from .regions import Region
from .geometry import (
    RectangularFilm,
    Geometry,
    VoltageProbe
)
from .contacts import Contact


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
    if "regions" in data:

        data["regions"] = [
            Region(**region)
            for region in data["regions"]
    ]
    if "contacts" in data:

        data["contacts"] = [
            Contact(**contact)
            for contact in data["contacts"]
    ]
    if "voltage_probes" in data:

        data["voltage_probes"] = [
            VoltageProbe(**probe)
            for probe in data["voltage_probes"]
    ]
        
    geometry_type = data.pop("type")

    if geometry_type == "rectangle":

        film = RectangularFilm(
            name=data["name"],
            width=data["width"],
            height=data["height"],
            thickness=data["thickness"],
            nx=data["nx"],
            ny=data["ny"],
            regions=data.get("regions"),
            contacts=data.get("contacts")
    )

        return Geometry(
            film=film,
            contacts=data.get("contacts", []),
            voltage_probes=data.get("voltage_probes", []),
            regions=data.get("regions", [])
    )


    raise ValueError(
        f"Unknown geometry type: {geometry_type}"
    )