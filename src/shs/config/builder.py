"""
Simulation construction pipeline.

Converts a simulation configuration file into a complete
SHS simulation object.

The builder is responsible for constructing the static
simulation infrastructure and initial physical state.

Physics models and numerical solvers remain responsible
for evolving that state.
"""

from pathlib import Path

from shs.config.simulation import load_simulation
from shs.config.simulation_state import Simulation

from shs.geometry.database import load_geometry
from shs.geometry.mesh import create_mesh

from shs.materials.database import get_material

from shs.mapping import (
    build_region_map,
    build_material_map,
    build_contact_map,
)

from shs.boundaries import (
    BoundarySet,
    BoundaryCondition,
    BoundarySide,
    BoundaryType,
)

from shs.physics.fields import Fields

from shs.tdgl.boundary import (
    TDGLBoundarySet,
    TDGLBoundaryCondition,
    TDGLBoundarySide,
    TDGLBoundaryType,
)

from shs.tdgl.initialization import equilibrium_superconducting_state


def build_simulation(filepath):
    """
    Build a complete SHS simulation from a JSON configuration.

    Parameters
    ----------
    filepath:
        Path to the simulation configuration.

    Returns
    -------
    Simulation
        Fully initialized simulation object.
    """

    config_root = Path("configs")

    config = load_simulation(filepath)

    #
    # Geometry.
    #

    geometry_path = (
        config_root /
        "geometry" /
        config.geometry
    )

    geometry = load_geometry(
        geometry_path
    )

    #
    # Mesh.
    #

    mesh = create_mesh(
        geometry
    )

    #
    # Region mapping.
    #

    region_map = build_region_map(
        geometry,
        mesh
    )

    #
    # Material mapping.
    #

    material_name = config.material.replace(
        ".json",
        ""
    )

    material = get_material(
        material_name
    )

    material_map = build_material_map(
        region_map,
        material
    )

    #
    # Contact mapping.
    #

    contact_map = build_contact_map(
        geometry,
        mesh
    )

    #
    # Initial superconducting order parameter.
    #

    reduced_temperature = (
        config.temperature /
        material_map.Tc
    )

    initial_psi = (
        equilibrium_superconducting_state(
            shape=(mesh.ny, mesh.nx),
            reduced_temperature=float(
                reduced_temperature[0, 0]
            ),
        )
    )

    #
    # Simulation fields.
    #

    fields = Fields.create(
        mesh=mesh,
        initial_temperature=config.temperature,
        initial_psi=initial_psi,
    )

    #
    # Thermal boundary conditions.
    #

    boundaries = BoundarySet()

    for side_name, values in (
        config.boundaries.items()
    ):
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

    #
    # TDGL boundary conditions.
    #
    # Unlike the previous implementation, use the
    # actual TDGL boundary configuration supplied by
    # the simulation JSON.
    #

    tdgl_boundaries = TDGLBoundarySet()

    for side_name, values in (
        config.tdgl_boundaries.items()
    ):
        tdgl_boundary = TDGLBoundaryCondition(
            side=TDGLBoundarySide(
                side_name
            ),

            type=TDGLBoundaryType(
                values["type"]
            ),
        )

        tdgl_boundaries.add(
            tdgl_boundary
        )

    #
    # Complete simulation object.
    #

    simulation = Simulation(
        config=config,

        geometry=geometry,

        mesh=mesh,

        region_map=region_map,

        material_map=material_map,

        boundaries=boundaries,

        fields=fields,

        contact_map=contact_map,

        tdgl_boundaries=tdgl_boundaries,
    )

    #
    # Validate the complete constructed simulation.
    #

    simulation.validate()

    return simulation