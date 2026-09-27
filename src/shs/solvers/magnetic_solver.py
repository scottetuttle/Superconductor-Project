"""Thin-film self-field update for the coupled TDGL solve."""

import numpy as np

from shs.physics.electromagnetics import (
    induced_vector_potential_from_sheet_current,
    perpendicular_magnetic_field,
)


def magnetic_screening_step(simulation):
    """Apply one damped fixed-point update of the induced vector potential."""
    fields = simulation.fields
    simulation_config = getattr(simulation, "config", None)
    config = getattr(simulation_config, "electromagnetic", None)
    if config is None:
        fields.magnetic_screening_iterations = 0
        fields.magnetic_screening_residual = 0.0
        return fields
    if not config.include_self_field:
        fields.magnetic_screening_iterations = 0
        fields.magnetic_screening_residual = 0.0
        return fields

    target_x, target_y = induced_vector_potential_from_sheet_current(
        simulation.mesh,
        fields.current_density_x * simulation.region_map.active_mask,
        fields.current_density_y * simulation.region_map.active_mask,
        simulation.material_map.thickness,
        config.screening_source_stride,
    )
    old_x = fields.induced_vector_potential_x
    old_y = fields.induced_vector_potential_y
    difference_norm = np.sqrt(np.mean((target_x - old_x) ** 2 + (target_y - old_y) ** 2))
    target_norm = np.sqrt(np.mean(target_x**2 + target_y**2))
    scale = max(target_norm, np.finfo(float).tiny)
    residual = float(difference_norm / scale)
    alpha = config.screening_step_size
    fields.induced_vector_potential_x = old_x + alpha * (target_x - old_x)
    fields.induced_vector_potential_y = old_y + alpha * (target_y - old_y)
    fields.vector_potential_x = (
        fields.applied_vector_potential_x + fields.induced_vector_potential_x
    )
    fields.vector_potential_y = (
        fields.applied_vector_potential_y + fields.induced_vector_potential_y
    )
    fields.magnetic_field_z = perpendicular_magnetic_field(
        fields.vector_potential_x,
        fields.vector_potential_y,
        simulation.mesh.dx,
        simulation.mesh.dy,
    )
    fields.magnetic_screening_iterations = 1
    fields.magnetic_screening_residual = residual
    return fields
