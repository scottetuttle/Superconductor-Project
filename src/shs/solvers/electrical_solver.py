"""Electrical continuity with SI fields and insulating non-contact faces.

Directional arrays store outgoing x/y link quantities at each node; the
last x column and last y row have no outgoing link. This matches the TDGL
link-current convention. Conductivity follows the existing phenomenological
normal-fraction model, not a microscopic quasiparticle transport model.
"""
import numpy as np
from shs.numerics.iterative import (
    red_black_sor, sparse_direct, sparse_neumann, face_coefficients,
)
from shs.utils.defaults import default_section


_DEFAULTS = default_section("electrical")


def _terminal_side(mask):
    """Return the single outer side containing an edge-terminal mask."""
    mask = np.asarray(mask, dtype=bool)
    candidates = []
    for name, edge in (
        ("left", np.s_[:, 0]), ("right", np.s_[:, -1]),
        ("bottom", np.s_[0, :]), ("top", np.s_[-1, :]),
    ):
        allowed = np.zeros_like(mask)
        allowed[edge] = True
        if np.any(mask) and np.all(~mask | allowed):
            candidates.append(name)
    if len(candidates) != 1:
        raise ValueError("A current terminal must span at least two nodes on one outer edge.")
    side = candidates[0]
    if side == "left":
        edge_values = mask[:, 0]
    elif side == "right":
        edge_values = mask[:, -1]
    elif side == "bottom":
        edge_values = mask[0, :]
    else:
        edge_values = mask[-1, :]
    coordinates = np.flatnonzero(edge_values)
    if coordinates.size < 2 or np.any(np.diff(coordinates) != 1):
        raise ValueError("A current terminal must contain a contiguous edge segment of at least two nodes.")
    return side, coordinates


def current_flux_source(mesh, material_map, contact_map, applied_current,
                        source_contact="left_current", sink_contact="right_current"):
    """Create the Poisson source for balanced terminal current fluxes.

    Positive current enters the device through ``source_contact`` and leaves
    through ``sink_contact``. Contact area uses trapezoidal edge integration;
    the collocated Neumann equation applies the resulting flux density at each
    contact node.
    """
    if not np.isfinite(applied_current):
        raise ValueError("Applied current must be finite.")
    if source_contact == sink_contact:
        raise ValueError("Current source and sink contacts must differ.")
    try:
        source_mask = contact_map.contact_masks[source_contact]
        sink_mask = contact_map.contact_masks[sink_contact]
    except KeyError as error:
        raise ValueError(f"Missing current terminal {error.args[0]}.") from error
    if np.any(source_mask & sink_mask):
        raise ValueError("Current source and sink contacts must not overlap.")
    thickness = np.asarray(material_map.thickness, dtype=float)
    if not np.all(np.isfinite(thickness)) or np.any(thickness <= 0):
        raise ValueError("Current drive requires positive finite film thickness.")
    result = np.zeros_like(thickness)
    terminal_data = []
    for mask, outward_sign in ((source_mask, -1.0), (sink_mask, 1.0)):
        side, coordinates = _terminal_side(mask)
        edge_index = ((coordinates, np.zeros_like(coordinates))
                      if side == "left" else
                      (coordinates, np.full_like(coordinates, mesh.nx - 1))
                      if side == "right" else
                      (np.zeros_like(coordinates), coordinates)
                      if side == "bottom" else
                      (np.full_like(coordinates, mesh.ny - 1), coordinates))
        area_weights = np.ones(coordinates.size)
        area_weights[[0, -1]] = 0.5
        spacing = mesh.dy if side in {"left", "right"} else mesh.dx
        normal_spacing = mesh.dx if side in {"left", "right"} else mesh.dy
        areas = thickness[edge_index] * area_weights * spacing
        area = float(np.sum(areas))
        current_density = applied_current / area
        # b has units A/m^3 in div(sigma grad V)=b.
        domain_weights = np.ones(coordinates.size)
        edge_extent = mesh.ny if side in {"left", "right"} else mesh.nx
        domain_weights[coordinates == 0] = 0.5
        domain_weights[coordinates == edge_extent - 1] = 0.5
        result[edge_index] += (
            outward_sign * 2.0 * current_density / normal_spacing
            * area_weights / domain_weights
        )
        terminal_data.append((side, area, current_density))
    return result, terminal_data


def link_divergence(jx, jy, dx, dy):
    result = np.zeros_like(jx, dtype=float)
    result[:, :-1] += jx[:, :-1]/dx
    result[:, 1:] -= jx[:, :-1]/dx
    result[:-1] += jy[:-1]/dy
    result[1:] -= jy[:-1]/dy
    return result


def electrical_step(fields, mesh, material_map, contact_map,
                    voltage_left=_DEFAULTS["voltage_left"], voltage_right=_DEFAULTS["voltage_right"],
                    superconducting_fraction=None,
                    superconducting_current_x=None, superconducting_current_y=None,
                    solver_tolerance=_DEFAULTS["solver"]["tolerance"],
                    solver_max_iterations=_DEFAULTS["solver"]["max_iterations"],
                    solver_omega=_DEFAULTS["solver"]["omega"],
                    residual_check_interval=_DEFAULTS["solver"]["residual_check_interval"],
                    solver_backend=_DEFAULTS["solver"]["backend"],
                    normal_conductivity_model=_DEFAULTS["normal_conductivity_model"],
                    drive_mode=_DEFAULTS["drive_mode"], applied_current=0.0,
                    source_contact=_DEFAULTS["source_contact"],
                    sink_contact=_DEFAULTS["sink_contact"],
                    reference_voltage=_DEFAULTS["reference_voltage"], active_mask=None):
    shape = fields.voltage.shape
    active = np.ones(shape, dtype=bool) if active_mask is None else np.asarray(active_mask, dtype=bool)
    if active.shape != shape or not np.any(active):
        raise ValueError("Electrical active mask must match the grid and contain active nodes.")
    left = contact_map.contact_masks.get(source_contact)
    right = contact_map.contact_masks.get(sink_contact)
    if left is None or right is None:
        raise ValueError('Current contacts missing.')
    if np.any(left & right):
        raise ValueError('Voltage contacts must not overlap.')
    if np.any((left | right) & ~active):
        raise ValueError("Electrical contacts cannot include inactive geometry nodes.")
    if not np.isfinite(voltage_left + voltage_right):
        raise ValueError('Contact voltages must be finite.')
    jsx = np.zeros(shape) if superconducting_current_x is None else np.asarray(superconducting_current_x)
    jsy = np.zeros(shape) if superconducting_current_y is None else np.asarray(superconducting_current_y)
    # Direct electrical calls retain their normal-state default. Coupled callers
    # explicitly supply the condensate fraction through electrical_simulation_step.
    fraction = np.zeros(shape) if superconducting_fraction is None else np.asarray(superconducting_fraction)
    if any(a.shape != shape or not np.all(np.isfinite(a)) for a in (jsx, jsy, fraction)):
        raise ValueError('Electrical transport fields must be finite and match the mesh.')
    if normal_conductivity_model == "constant":
        sigma = material_map.electrical_conductivity
    elif normal_conductivity_model == "condensate_depletion":
        sigma = material_map.electrical_conductivity * np.clip(1-fraction, 0, 1)
    else:
        raise ValueError("Unknown normal conductivity model.")
    sigma = np.where(active, sigma, 0.0)
    jsx = np.where(active, jsx, 0.0)
    jsy = np.where(active, jsy, 0.0)
    jsx[:, :-1] *= active[:, :-1] & active[:, 1:]
    jsy[:-1, :] *= active[:-1, :] & active[1:, :]
    values = np.zeros(shape)
    # div(-sigma grad V + Js) = 0 => div(sigma grad V) = div Js.
    source = link_divergence(jsx, jsy, mesh.dx, mesh.dy)
    if drive_mode == "voltage":
        mask = left | right
        values[left], values[right] = voltage_left, voltage_right
    elif drive_mode == "current":
        if solver_backend != "sparse_direct":
            raise ValueError("Current drive currently requires the sparse_direct backend.")
        flux_source, _ = current_flux_source(
            mesh, material_map, contact_map, applied_current,
            source_contact, sink_contact,
        )
        source = source + flux_source
        mask = np.zeros(shape, dtype=bool)
    else:
        raise ValueError("Electrical drive_mode must be voltage or current.")
    if drive_mode == "current":
        factorization_cache = getattr(material_map, "_neumann_factorization_cache", None)
        if factorization_cache is None:
            factorization_cache = {}
            material_map._neumann_factorization_cache = factorization_cache
        result = sparse_neumann(
            fields.voltage, sigma, source, mesh.dx, mesh.dy,
            reference_voltage, solver_tolerance, material_map.thickness, active,
            factorization_cache=factorization_cache,
        )
    elif solver_backend == "red_black_sor":
        result = red_black_sor(
            fields.voltage, sigma, source, mask, values,
            mesh.dx, mesh.dy, solver_tolerance,
            solver_max_iterations, solver_omega,
            residual_check_interval,
        )
    elif solver_backend == "sparse_direct":
        result = sparse_direct(
            fields.voltage, sigma, source, mask, values,
            mesh.dx, mesh.dy, solver_tolerance,
        )
    else:
        raise ValueError(f"Unknown electrical solver backend: {solver_backend}")
    if not result.converged:
        raise RuntimeError(f'Electrical solve failed after {result.iterations} iterations '
                           f'(equation residual {result.residual:.3g} V).')
    ex, ey = np.zeros(shape), np.zeros(shape)
    sx, sy = face_coefficients(sigma, mesh.dx, mesh.dy)
    ex[:, :-1] = -np.diff(result.field, axis=1)/mesh.dx
    ey[:-1] = -np.diff(result.field, axis=0)/mesh.dy
    ex[:, :-1] *= active[:, :-1] & active[:, 1:]
    ey[:-1, :] *= active[:-1, :] & active[1:, :]
    ex[~active] = 0.0
    ey[~active] = 0.0
    if not np.any(sigma):
        # Existing ideal nondissipative limit: potential is underdetermined.
        ex.fill(0)
        ey.fill(0)
    jnx, jny = np.zeros(shape), np.zeros(shape)
    jnx[:, :-1] = sx*mesh.dx**2*ex[:, :-1]
    jny[:-1] = sy*mesh.dy**2*ey[:-1]
    joule = jnx*ex + jny*ey
    external = fields.external_heat_source
    if external is None:
        # Preserve legacy heat_source assignments on the first electrical solve.
        external = fields.heat_source.copy() if fields.joule_heat_source is None else np.zeros(shape)
    fields.external_heat_source = np.asarray(external).copy()
    fields.joule_heat_source = joule
    laser = (0.0 if fields.laser_heat_source is None else fields.laser_heat_source)
    fields.heat_source = fields.external_heat_source + joule + laser
    fields.voltage = result.field
    fields.electric_field_x, fields.electric_field_y = ex, ey
    fields.normal_current_density_x, fields.normal_current_density_y = jnx, jny
    fields.current_density_x, fields.current_density_y = jnx+jsx, jny+jsy
    fields.electrical_solver_iterations = result.iterations
    fields.electrical_solver_residual = result.residual
    fields.electrical_solver_converged = result.converged
    return fields


def electrical_simulation_step(simulation, voltage_left=None, voltage_right=None,
                               include_superconductivity=True, solver_tolerance=None):
    """Single configuration-aware transport interface for all coupled solvers."""
    config = getattr(simulation.config, 'electrical', None)
    solver = getattr(config, 'solver', None)
    fields = simulation.fields
    return electrical_step(
        fields, simulation.mesh, simulation.material_map, simulation.contact_map,
        voltage_left=getattr(config, 'voltage_left', _DEFAULTS['voltage_left']) if voltage_left is None else voltage_left,
        voltage_right=getattr(config, 'voltage_right', _DEFAULTS['voltage_right']) if voltage_right is None else voltage_right,
        superconducting_fraction=np.abs(fields.psi)**2 if include_superconductivity else None,
        superconducting_current_x=fields.supercurrent_density_x if include_superconductivity else None,
        superconducting_current_y=fields.supercurrent_density_y if include_superconductivity else None,
        solver_tolerance=getattr(solver, 'tolerance', _DEFAULTS['solver']['tolerance']) if solver_tolerance is None else solver_tolerance,
        solver_max_iterations=getattr(solver, 'max_iterations', _DEFAULTS['solver']['max_iterations']),
        solver_omega=getattr(solver, 'omega', _DEFAULTS['solver']['omega']),
        residual_check_interval=getattr(
            solver,
            'residual_check_interval',
            _DEFAULTS['solver']['residual_check_interval'],
        ),
        solver_backend=getattr(solver, 'backend', _DEFAULTS['solver']['backend']),
        normal_conductivity_model=getattr(
            config,
            'normal_conductivity_model',
            _DEFAULTS['normal_conductivity_model'],
        ),
        drive_mode=getattr(config, 'drive_mode', _DEFAULTS['drive_mode']),
        applied_current=getattr(simulation.config, 'current', 0.0),
        source_contact=getattr(config, 'source_contact', _DEFAULTS['source_contact']),
        sink_contact=getattr(config, 'sink_contact', _DEFAULTS['sink_contact']),
        reference_voltage=getattr(config, 'reference_voltage', _DEFAULTS['reference_voltage']),
        active_mask=simulation.region_map.active_mask,
    )
