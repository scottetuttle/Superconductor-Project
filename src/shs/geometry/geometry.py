"""
SHS Geometry Module

Defines physical device geometries.

Currently supported:
- Rectangular superconducting films

Future:
- Arbitrary polygons
- Multiple layers
- Contacts
- Defects
- Patterned devices
"""

from dataclasses import dataclass
from typing import List
from .regions import Region


@dataclass
class Contact:
    """
    Represents an electrical contact attached to the device.
    """

    name: str

    x: float
    y: float

    width: float
    height: float


@dataclass
class VoltageProbe:
    """
    Represents a voltage measurement location.
    """

    name: str

    x: float
    y: float


@dataclass(frozen=True)
class Hole:
    """An insulating cutout removed from the superconducting film.

    Coordinates and dimensions use SI metres. Circular holes use ``radius``;
    rectangular holes use ``width`` and ``height``.
    """

    name: str
    shape: str
    x: float
    y: float
    radius: float | None = None
    width: float | None = None
    height: float | None = None


@dataclass
class RectangularFilm:

    name: str

    width: float
    height: float
    thickness: float

    nx: int
    ny: int

    regions: list[Region] | None = None
    contacts: list[Contact] | None = None
    holes: list[Hole] | None = None


@dataclass
class Geometry:
    """
    Complete description of a simulated device.
    """

    film: RectangularFilm

    contacts: List[Contact]

    voltage_probes: List[VoltageProbe]

    regions: List[Region]

    holes: List[Hole] | None = None
