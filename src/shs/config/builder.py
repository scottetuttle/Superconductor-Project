"""
Simulation construction pipeline.

Converts a simulation configuration file into
a complete SHS simulation object.
"""

from shs.config.simulation import load_simulation
from shs.config.simulation_state import Simulation

from shs.geometry.database import load_geometry
from shs.geometry.mesh import create_mesh

from shs.materials.database import get_material

from shs.mapping import (
    build_region_map,
    build_material_map,
)

from shs.boundaries import BoundarySet
from shs.boundaries.boundary import (
    BoundaryCondition,
    BoundarySide,
    BoundaryType,
)

from pathlib import Path

from shs.physics.fields import Fields


def build_simulation(filepath):
    CONFIG_ROOT = Path("configs")
    config = load_simulation(filepath)


    # Geometry

    geometry_path = (
        CONFIG_ROOT /
        "geometry" /
        config.geometry
    )

    geometry = load_geometry(
        geometry_path
    )


    mesh = create_mesh(
        geometry
    )


    # Region mapping

    region_map = build_region_map(
        geometry,
        mesh
    )


    # Material mapping

    material_path = (
        CONFIG_ROOT /
        "materials" /
        config.material
    )

    material_name = (
        config.material.replace(
            ".json",
            ""
        )
    )

    material = get_material(
        material_name
    )


    material_map = build_material_map(
        region_map,
        material
    )


    # Boundaries

    boundaries = BoundarySet()

    for side_name, values in config.boundaries.items():

        # Runtime fields

        fields = Fields.create(
        mesh=mesh,
        initial_temperature=config.temperature,
)

        boundary = BoundaryCondition(

            side=BoundarySide(side_name),

            type=BoundaryType(
                values["type"]
            ),

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

        boundaries.add(
            boundary
        )



    return Simulation(

        config=config,

        geometry=geometry,

        mesh=mesh,

        region_map=region_map,

        material_map=material_map,

        boundaries=boundaries,

        fields=fields,
    )