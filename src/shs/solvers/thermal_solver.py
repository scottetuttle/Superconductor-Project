"""Conservative explicit heat diffusion on the rectangular node grid.

Boundary nodes own half-width control volumes. Heat flux is positive inward;
convection uses the boundary's ambient temperature. Missing boundaries default
to insulation for programmatically constructed simulations.
"""
import numpy as np
from shs.boundaries import BoundarySide, BoundaryType
from shs.numerics.iterative import face_coefficients
from shs.utils.defaults import default_section


def thermal_step(simulation, dt, thermal_model, *, initial_temperature=None,
                 initial_phonon_temperature=None):
    if thermal_model.model == "two_temperature":
        from shs.solvers.two_temperature_solver import two_temperature_step
        return two_temperature_step(simulation, dt, thermal_model,
                                    initial_temperature=initial_temperature,
                                    initial_phonon_temperature=initial_phonon_temperature)
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('Thermal timestep must be positive and finite.')
    thermal_model.validate()
    fields, mesh, materials = simulation.fields, simulation.mesh, simulation.material_map
    k, capacity = materials.thermal_conductivity, materials.heat_capacity
    active = np.asarray(simulation.region_map.active_mask, dtype=bool)
    if active.shape != fields.temperature.shape or not np.any(active):
        raise ValueError("Thermal active mask must match the grid and contain active nodes.")
    if np.any(capacity <= 0) or np.any(k < 0) or not np.all(np.isfinite(k+capacity)):
        raise ValueError('Heat capacity must be positive and conductivity nonnegative and finite.')
    temperature = np.array(fields.temperature if initial_temperature is None else initial_temperature, copy=True)
    cx, cy = face_coefficients(k, mesh.dx, mesh.dy)
    cx *= active[:, :-1] & active[:, 1:]
    cy *= active[:-1, :] & active[1:, :]
    wx, wy = np.ones(mesh.nx), np.ones(mesh.ny)
    wx[[0, -1]], wy[[0, -1]] = .5, .5
    diagonal = np.zeros_like(temperature)
    diagonal[:, :-1] += cx/wx[:-1]
    diagonal[:, 1:] += cx/wx[1:]
    diagonal[:-1] += cy/wy[:-1, None]
    diagonal[1:] += cy/wy[1:, None]
    sides = {
        BoundarySide.LEFT: (np.s_[:, 0], mesh.dx/2),
        BoundarySide.RIGHT: (np.s_[:, -1], mesh.dx/2),
        BoundarySide.BOTTOM: (np.s_[0, :], mesh.dy/2),
        BoundarySide.TOP: (np.s_[-1, :], mesh.dy/2),
    }
    fixed_mask = np.zeros_like(temperature, dtype=bool)
    fixed_values = np.zeros_like(temperature)
    boundaries = list(simulation.boundaries) if simulation.boundaries is not None else []
    for boundary in boundaries:
        index, width = sides[boundary.side]
        if boundary.type == BoundaryType.FIXED_TEMPERATURE:
            value = boundary.temperature
            if value is None or not np.isfinite(value) or value < 0:
                raise ValueError('Fixed boundary temperature must be nonnegative and finite.')
            if np.any(fixed_mask[index] & (fixed_values[index] != value)):
                raise ValueError('Conflicting fixed temperatures at a corner.')
            fixed_mask[index], fixed_values[index] = True, value
        elif boundary.type == BoundaryType.CONVECTION:
            if (boundary.heat_transfer is None or boundary.heat_transfer < 0
                    or boundary.temperature is None or boundary.temperature < 0
                    or not np.isfinite(boundary.heat_transfer + boundary.temperature)):
                raise ValueError('Convection requires finite nonnegative h and ambient temperature.')
            diagonal[index] += boundary.heat_transfer/width
        elif boundary.type == BoundaryType.FIXED_HEAT_FLUX:
            if boundary.heat_flux is None or not np.isfinite(boundary.heat_flux):
                raise ValueError('Heat flux must be finite (positive inward).')
        elif boundary.type != BoundaryType.INSULATING:
            raise ValueError(f'Unsupported thermal boundary: {boundary.type}')
    fixed_mask &= active
    rate = diagonal/capacity + thermal_model.thermal_relaxation_rate
    max_rate = np.max(rate[active & ~fixed_mask], initial=0)
    stable_dt = thermal_model.stability_safety_factor/max_rate if max_rate else np.inf
    limit = min(thermal_model.max_substep, stable_dt)
    count = max(1, int(np.ceil(dt/limit)))
    numerics = getattr(getattr(simulation, "config", None), "numerics", None)
    max_substeps = (numerics.max_internal_substeps if numerics is not None
                    else default_section("numerics")["max_internal_substeps"])
    if count > max_substeps:
        raise ValueError('Thermal step requires over one million substeps; check dt and max_substep units.')
    sub_dt = dt/count
    heat = fields.heat_source
    if fields.external_heat_source is not None:
        heat = (
            fields.external_heat_source
            + (0 if fields.joule_heat_source is None else fields.joule_heat_source)
            + (0 if fields.laser_heat_source is None else fields.laser_heat_source)
        )
    if not np.all(np.isfinite(temperature)) or not np.all(np.isfinite(heat)):
        raise ValueError('Thermal inputs must be finite.')
    temperature[fixed_mask] = fixed_values[fixed_mask]
    inactive_temperature = temperature[~active].copy()
    for _ in range(count):
        power = np.zeros_like(temperature)
        fx = cx*np.diff(temperature, axis=1)
        fy = cy*np.diff(temperature, axis=0)
        power[:, :-1] += fx/wx[:-1]
        power[:, 1:] -= fx/wx[1:]
        power[:-1] += fy/wy[:-1, None]
        power[1:] -= fy/wy[1:, None]
        for boundary in boundaries:
            index, width = sides[boundary.side]
            if boundary.type == BoundaryType.CONVECTION:
                power[index] += boundary.heat_transfer*(boundary.temperature-temperature[index])/width
            elif boundary.type == BoundaryType.FIXED_HEAT_FLUX:
                power[index] += boundary.heat_flux/width
        derivative = (power+heat)/capacity - thermal_model.thermal_relaxation_rate*(temperature-thermal_model.bath_temperature)
        temperature[active] += sub_dt*derivative[active]
        temperature[fixed_mask] = fixed_values[fixed_mask]
        temperature[~active] = inactive_temperature
        if not np.all(np.isfinite(temperature)):
            raise RuntimeError('Thermal solver produced nonfinite temperature.')
    fields.temperature = temperature
    fields.thermal_solver_substeps = count
    return fields
