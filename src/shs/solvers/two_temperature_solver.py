"""Conservative two-bath film heat evolution with exact local e-ph exchange.

Material heat capacity and conductivity are partitioned by configured fractions.
Those fractions are exploratory until calibrated material-specific values exist.
"""
import numpy as np

from shs.boundaries import BoundarySide, BoundaryType
from shs.numerics.iterative import face_coefficients
from shs.utils.defaults import default_section


def two_temperature_step(simulation, dt, model, *, initial_temperature=None,
                         initial_phonon_temperature=None):
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("Thermal timestep must be positive and finite.")
    model.validate()
    fields, mesh, material = simulation.fields, simulation.mesh, simulation.material_map
    active = np.asarray(simulation.region_map.active_mask, dtype=bool)
    ce = material.heat_capacity * model.electron_heat_capacity_fraction
    cp = material.heat_capacity * (1-model.electron_heat_capacity_fraction)
    ke = material.thermal_conductivity * model.electron_thermal_conductivity_fraction
    kp = material.thermal_conductivity * (1-model.electron_thermal_conductivity_fraction)
    if (np.any(ce <= 0) or np.any(cp <= 0) or np.any(ke < 0) or np.any(kp < 0)
            or not np.all(np.isfinite(ce+cp+ke+kp))):
        raise ValueError("Two-temperature capacities and conductivities must be physical.")
    te = np.array(fields.temperature if initial_temperature is None
                  else initial_temperature, dtype=float, copy=True)
    phonon_source = fields.phonon_temperature
    if phonon_source is None:
        phonon_source = fields.temperature
    tp = np.array(phonon_source if initial_phonon_temperature is None
                  else initial_phonon_temperature, dtype=float, copy=True)
    if te.shape != active.shape or tp.shape != active.shape or not np.any(active):
        raise ValueError("Two-temperature fields must match the active grid.")
    heat = fields.heat_source
    if fields.external_heat_source is not None:
        heat = (fields.external_heat_source
                + (0 if fields.joule_heat_source is None else fields.joule_heat_source)
                + (0 if fields.laser_heat_source is None else fields.laser_heat_source))
    if not all(np.all(np.isfinite(value)) for value in (te, tp, heat)):
        raise ValueError("Two-temperature inputs must be finite.")
    wx, wy = np.ones(mesh.nx), np.ones(mesh.ny)
    wx[[0, -1]] = .5
    wy[[0, -1]] = .5
    faces = []
    max_rates = []
    for conductivity, capacity in ((ke, ce), (kp, cp)):
        cx, cy = face_coefficients(conductivity, mesh.dx, mesh.dy)
        cx *= active[:, :-1] & active[:, 1:]
        cy *= active[:-1, :] & active[1:, :]
        diagonal = np.zeros_like(te)
        diagonal[:, :-1] += cx/wx[:-1]
        diagonal[:, 1:] += cx/wx[1:]
        diagonal[:-1] += cy/wy[:-1, None]
        diagonal[1:] += cy/wy[1:, None]
        max_rates.append(float(np.max((diagonal/capacity)[active])))
        faces.append((cx, cy))
    sides = {
        BoundarySide.LEFT: (np.s_[:, 0], mesh.dx/2),
        BoundarySide.RIGHT: (np.s_[:, -1], mesh.dx/2),
        BoundarySide.BOTTOM: (np.s_[0, :], mesh.dy/2),
        BoundarySide.TOP: (np.s_[-1, :], mesh.dy/2),
    }
    boundaries = list(simulation.boundaries) if simulation.boundaries is not None else []
    fixed = np.zeros_like(active)
    fixed_values = np.zeros_like(te)
    for boundary in boundaries:
        index, width = sides[boundary.side]
        if boundary.type == BoundaryType.FIXED_TEMPERATURE:
            if boundary.temperature is None or not np.isfinite(boundary.temperature) or boundary.temperature < 0:
                raise ValueError("Fixed thermal boundary must have finite nonnegative temperature.")
            if np.any(fixed[index] & (fixed_values[index] != boundary.temperature)):
                raise ValueError("Conflicting fixed temperatures at a corner.")
            fixed[index] = True
            fixed_values[index] = boundary.temperature
        elif boundary.type == BoundaryType.CONVECTION:
            if (boundary.heat_transfer is None or boundary.temperature is None
                    or not np.isfinite(boundary.heat_transfer + boundary.temperature)
                    or boundary.heat_transfer < 0 or boundary.temperature < 0):
                raise ValueError("Convection requires finite nonnegative coefficients.")
            max_rates[1] += boundary.heat_transfer/width/float(np.min(cp[active]))
        elif boundary.type == BoundaryType.FIXED_HEAT_FLUX:
            if boundary.heat_flux is None or not np.isfinite(boundary.heat_flux):
                raise ValueError("Fixed heat flux must be finite.")
        elif boundary.type != BoundaryType.INSULATING:
            raise ValueError("Unsupported thermal boundary.")
    fixed &= active
    max_rate = max(max_rates[0], max_rates[1]+model.phonon_escape_rate_per_s)
    stable_dt = model.stability_safety_factor/max_rate if max_rate else np.inf
    count = max(1, int(np.ceil(dt/min(model.max_substep, stable_dt))))
    numerics = getattr(getattr(simulation, "config", None), "numerics", None)
    maximum = (numerics.max_internal_substeps if numerics is not None
               else default_section("numerics")["max_internal_substeps"])
    if count > maximum:
        raise ValueError("Two-temperature step exceeds internal substep budget.")
    sub_dt = dt/count
    inactive_e, inactive_p = te[~active].copy(), tp[~active].copy()
    te[fixed] = fixed_values[fixed]
    tp[fixed] = fixed_values[fixed]

    def divergence(temperature, faces):
        cx, cy = faces
        power = np.zeros_like(temperature)
        fx, fy = cx*np.diff(temperature, axis=1), cy*np.diff(temperature, axis=0)
        power[:, :-1] += fx/wx[:-1]
        power[:, 1:] -= fx/wx[1:]
        power[:-1] += fy/wy[:-1, None]
        power[1:] -= fy/wy[1:, None]
        return power

    coupling = model.electron_phonon_coupling_W_m3_K
    exchange_rate = coupling*(1/ce+1/cp)
    for _ in range(count):
        electron_power = divergence(te, faces[0]) + heat
        phonon_power = divergence(tp, faces[1])
        for boundary in boundaries:
            index, width = sides[boundary.side]
            if boundary.type == BoundaryType.CONVECTION:
                phonon_power[index] += boundary.heat_transfer*(boundary.temperature-tp[index])/width
            elif boundary.type == BoundaryType.FIXED_HEAT_FLUX:
                phonon_power[index] += boundary.heat_flux/width
        te[active] += sub_dt*(electron_power/ce)[active]
        tp[active] += sub_dt*(phonon_power/cp)[active]
        # Local exchange is exact for constant capacities: no stiff e-ph step limit.
        energy = ce*te + cp*tp
        delta = (te-tp)*np.exp(-exchange_rate*sub_dt)
        te[active] = ((energy+cp*delta)/(ce+cp))[active]
        tp[active] = ((energy-ce*delta)/(ce+cp))[active]
        tp[active] = (model.bath_temperature + (tp[active]-model.bath_temperature)
                      * np.exp(-model.phonon_escape_rate_per_s*sub_dt))
        te[fixed], tp[fixed] = fixed_values[fixed], fixed_values[fixed]
        te[~active], tp[~active] = inactive_e, inactive_p
        if not np.all(np.isfinite(te)) or not np.all(np.isfinite(tp)):
            raise RuntimeError("Two-temperature solver produced nonfinite temperature.")
    fields.temperature = te
    fields.phonon_temperature = tp
    fields.thermal_solver_substeps = count
    return fields
