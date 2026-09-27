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

import numpy as np

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

from shs.tdgl.model import TDGLModel
from shs.tdgl.parameters import TDGLParameters


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
        mesh,
        roundoff_ulps=config.numerics.contact_roundoff_ulps,
    )

    #
    # Initial superconducting order parameter.
    #

    reduced_temperature = (
        config.temperature /
        material_map.Tc
    )

    initial_model = TDGLModel(TDGLParameters(
        u=config.tdgl.u,
        time_integrator=config.tdgl.time_integrator,
        gamma=config.tdgl.gamma,
        kappa=config.tdgl.kappa,
        normalization=config.tdgl.normalization,
        temperature_model=config.tdgl.temperature_model,
        include_scalar_potential=config.tdgl.include_scalar_potential,
        max_normalized_timestep=config.tdgl.max_normalized_timestep,
        stability_safety_factor=config.tdgl.stability_safety_factor,
    ))
    initial_amplitude = initial_model.equilibrium_amplitude(
        float(reduced_temperature[0, 0])
    )
    initial_psi = np.full((mesh.ny, mesh.nx), initial_amplitude, dtype=complex)
    initial_psi[~region_map.active_mask] = 0.0

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
