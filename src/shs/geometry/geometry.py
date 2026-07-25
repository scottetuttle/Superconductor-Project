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
class RectangularFilm:

    name: str

    width: float
    height: float
    thickness: float

    nx: int
    ny: int

    regions: list[Region] = None



film = RectangularFilm(
    name = 'NbN_film',
    width=5e-6,
    height=5e-6,
    thickness=100e-9,
    nx=100,
    ny=100
)



device = RectangularFilm(
    name="NbN device",
    width=5e-6,
    height=5e-6,
    thickness=100e-9,
    nx=100,
    ny=100,
    regions=[
        Region(
            name="film",
            region_type="superconductor",
            material="NbN"
        )
    ]
)