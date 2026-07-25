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


@dataclass
class RectangularFilm:
    """
    Represents a rectangular superconducting thin film.

    Parameters
    ----------
    width:
        Physical width in meters.

    height:
        Physical height in meters.

    thickness:
        Film thickness in meters.

    nx:
        Number of grid cells in x direction.

    ny:
        Number of grid cells in y direction.
    """

    name: str

    width: float
    height: float
    thickness: float

    nx: int
    ny: int


    @property
    def dx(self):
        """
        Grid spacing in x direction.
        """
        return self.width / self.nx


    @property
    def dy(self):
        """
        Grid spacing in y direction.
        """
        return self.height / self.ny

film = RectangularFilm(
    name = 'NbN_film',
    width=5e-6,
    height=5e-6,
    thickness=100e-9,
    nx=100,
    ny=100
)