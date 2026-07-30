from dataclasses import dataclass

from shs.config.simulation import SimulationConfig

from shs.geometry.geometry import Geometry
from shs.geometry.mesh import Mesh

from shs.mapping.region_map import RegionMap
from shs.mapping.material_map import MaterialMap

from shs.boundaries.boundary import BoundarySet

from .validator import validate_simulation

from shs.physics.fields import Fields

from shs.mapping.contact_map import ContactMap


@dataclass
class Simulation:
    """
    Fully initialized SHS simulation.
    """

    config: SimulationConfig

    geometry: Geometry

    mesh: Mesh

    region_map: RegionMap

    material_map: MaterialMap

    boundaries: BoundarySet

    fields: Fields

    contact_map: ContactMap


    def validate(self):
        """
        Validate the simulation before running any solver.
        """

        return validate_simulation(self)