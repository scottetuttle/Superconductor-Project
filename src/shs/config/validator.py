"""
Simulation validation.

Checks that all runtime objects are internally consistent
before physics solvers begin.
"""

import numpy as np


def validate_simulation(simulation):
    """
    Validate a Simulation object.

    Raises
    ------
    ValueError
        If any inconsistency is detected.
    """

    mesh_shape = (
        simulation.mesh.ny,
        simulation.mesh.nx,
    )

    # ---------- Region Map ----------

    if simulation.region_map.region_ids.shape != mesh_shape:
        raise ValueError(
            "RegionMap shape does not match mesh."
        )

    # ---------- Material Map ----------

    material_map = simulation.material_map

    arrays = [
        material_map.material_ids,
        material_map.thermal_conductivity,
        material_map.heat_capacity,
        material_map.normal_resistivity,
        material_map.thickness,
        material_map.Tc,
        material_map.coherence_length,
        material_map.penetration_depth,
    ]

    for array in arrays:

        if array.shape != mesh_shape:
            raise ValueError(
                "MaterialMap array shape does not match mesh."
            )

    # ---------- Material IDs ----------

    unique_ids = np.unique(
        material_map.material_ids
    )

    for material_id in unique_ids:

        if material_id not in material_map.materials:
            raise ValueError(
                f"Unknown material id: {material_id}"
            )
    # ---------- Fields ----------

    field_arrays = [

        simulation.fields.temperature,

        simulation.fields.voltage,

        simulation.fields.current_density_x,

        simulation.fields.current_density_y,

        simulation.fields.heat_source,
]

    for array in field_arrays:

        if array.shape != mesh_shape:

            raise ValueError(
                "Field shape does not match mesh."
        )

    # ---------- Boundaries ----------

    if len(simulation.boundaries.boundaries) == 0:
        raise ValueError(
            "Simulation has no boundary conditions."
        )

    return True