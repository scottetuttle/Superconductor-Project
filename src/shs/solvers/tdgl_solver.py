"""Explicit TDGL evolution. Public dt is seconds; shared A and J use SI.

The existing SHS time convention tau = pi*hbar/(8*k_B*Tc) is retained.
The optional electrostatic formulation includes scalar potential through a
temporal link. Nonzero gamma enables the Kramer-Watts-Tobin amplitude update.
Self-consistent vector-potential evolution remains outside this solver.
"""
import numpy as np
from shs.tdgl.scaling import TDGLScales
from shs.tdgl.operators import covariant_laplacian
from shs.tdgl.boundary import (
    apply_insulating_boundary,
    apply_normal_contact_boundary,
    apply_normal_contact_mask,
    TDGLBoundarySide,
    TDGLBoundaryType,
)
from shs.physics.pinning import pinning_suppression
from shs.physics.josephson import weak_link_suppression
from shs.utils.defaults import default_section


def tdgl_scales(simulation, model):
    material = simulation.material_map.materials[0]
    conductivity = 1.0 / material.normal_resistivity
    return TDGLScales(material.coherence_length, material.penetration_depth,
                      material.Tc, model.characteristic_time(
                          material.Tc, conductivity, material.penetration_depth
                      ), normalization=model.parameters.normalization)


def generalized_tdgl_update(psi, rhs, scalar_potential, dt, u, gamma):
    """Advance one local generalized-TDGL step using a temporal link.

    This is the stable quadratic update for the Kramer-Watts-Tobin amplitude
    correction. With gamma=0 it reduces to phase-rotated explicit Euler.
    """
    phase = np.exp(-1j * scalar_potential * dt)
    amplitude_squared = np.abs(psi) ** 2
    gamma_squared = gamma**2
    z = 0.5 * gamma_squared * phase * psi
    w = z * amplitude_squared + phase * (
        psi
        + (dt / u)
        * np.sqrt(1.0 + gamma_squared * amplitude_squared)
        * rhs
    )
    coupling = np.real(np.conjugate(z) * w)
    leading = 2.0 * coupling + 1.0
    discriminant = (
        leading**2
        - 4.0 * np.abs(z) ** 2 * np.abs(w) ** 2
    )
    roundoff = 64.0 * np.finfo(float).eps * np.maximum(leading**2, 1.0)
    if np.any(discriminant < -roundoff):
        raise RuntimeError(
            "Generalized TDGL update has a negative amplitude discriminant; "
            "reduce the normalized timestep."
        )
    discriminant = np.maximum(discriminant, 0.0)
    next_amplitude_squared = (
        2.0 * np.abs(w) ** 2
        / (leading + np.sqrt(discriminant))
    )
    return w - z * next_amplitude_squared


def tdgl_step(simulation, dt, tdgl_model, *, initial_psi=None):
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('TDGL timestep must be positive finite seconds.')
    tdgl_model.parameters.validate()
    fields, mesh = simulation.fields, simulation.mesh
    scales = tdgl_scales(simulation, tdgl_model)
    dx, dy = mesh.dx/scales.xi, mesh.dy/scales.xi
    ax = scales.vector_potential_to_dimensionless(fields.vector_potential_x)
    ay = scales.vector_potential_to_dimensionless(fields.vector_potential_y)
    scalar_potential = (
        scales.scalar_potential_to_dimensionless(fields.voltage)
        if tdgl_model.parameters.include_scalar_potential
        else np.zeros_like(fields.voltage)
    )
    reduced = fields.temperature/simulation.material_map.Tc
    epsilon = (
        tdgl_model.temperature_coefficient(reduced)
        - pinning_suppression(simulation)
        - weak_link_suppression(simulation)
    )
    psi = np.array(fields.psi if initial_psi is None else initial_psi, dtype=complex, copy=True)
    active_mask = np.asarray(simulation.region_map.active_mask, dtype=bool)
    psi[~active_mask] = 0.0
    checked = (psi, reduced, ax, ay, scalar_potential)
    if not all(np.all(np.isfinite(a)) for a in checked):
        raise ValueError('TDGL fields must be finite.')
    boundaries = list(simulation.tdgl_boundaries)
    supported = {
        TDGLBoundaryType.INSULATING,
        TDGLBoundaryType.NORMAL_CONTACT,
        TDGLBoundaryType.NORMAL_CONTACT_MASK,
    }
    if any(b.type not in supported for b in boundaries):
        raise ValueError('Unsupported TDGL boundary type.')
    full_contact_sides = [
        b.side for b in boundaries if b.type == TDGLBoundaryType.NORMAL_CONTACT
    ]
    masked_contact_sides = [
        b.side for b in boundaries if b.type == TDGLBoundaryType.NORMAL_CONTACT_MASK
    ]
    insulating_sides = [side for side in TDGLBoundarySide if side not in full_contact_sides]
    contact_mask = np.zeros_like(psi, dtype=bool)
    if masked_contact_sides:
        side_indices = {
            TDGLBoundarySide.LEFT: (slice(None), 0),
            TDGLBoundarySide.RIGHT: (slice(None), -1),
            TDGLBoundarySide.BOTTOM: (0, slice(None)),
            TDGLBoundarySide.TOP: (-1, slice(None)),
        }
        for name, mask in simulation.contact_map.contact_masks.items():
            if simulation.contact_map.contact_types.get(name) != "current":
                continue
            for side in masked_contact_sides:
                index = side_indices[side]
                contact_mask[index] |= np.asarray(mask, dtype=bool)[index]
        if not np.any(contact_mask):
            raise ValueError(
                "A normal-contact TDGL boundary requires a current-contact mask "
                "on the configured mesh side."
            )
    dt_normalized = float(scales.time_to_dimensionless(dt))
    elapsed = 0.0
    count = 0
    while elapsed < dt_normalized:
        psi = apply_insulating_boundary(psi, ax, ay, dx, dy, insulating_sides)
        if full_contact_sides:
            psi = apply_normal_contact_boundary(psi, full_contact_sides)
        if masked_contact_sides:
            psi = apply_normal_contact_mask(psi, contact_mask)
        # Conservative explicit bound including diffusion and local cubic slope.
        gamma_factor = np.sqrt(
            1.0 + tdgl_model.parameters.gamma**2 * np.max(np.abs(psi) ** 2)
        )
        rate = gamma_factor * (
            2/dx**2 + 2/dy**2 + np.max(np.abs(epsilon))
            + 3*np.max(np.abs(psi)**2)
        ) / tdgl_model.parameters.u
        limit = min(tdgl_model.parameters.max_normalized_timestep,
                    tdgl_model.parameters.stability_safety_factor/rate)
        numerics = getattr(getattr(simulation, "config", None), "numerics", None)
        max_substeps = (numerics.max_internal_substeps if numerics is not None
                        else default_section("numerics")["max_internal_substeps"])
        if (dt_normalized-elapsed)/limit > max_substeps-count:
            raise ValueError('TDGL step requires over one million substeps; dt must be in seconds.')
        sub_dt = min(limit, dt_normalized-elapsed)
        rhs = (
            covariant_laplacian(psi, ax, ay, dx, dy, active_mask)
            + epsilon*psi
            - np.abs(psi)**2*psi
        )
        predictor = generalized_tdgl_update(
            psi,
            rhs,
            scalar_potential,
            sub_dt,
            tdgl_model.parameters.u,
            tdgl_model.parameters.gamma,
        )
        if tdgl_model.parameters.time_integrator == "heun":
            predictor[~active_mask] = 0.0
            predictor = apply_insulating_boundary(
                predictor, ax, ay, dx, dy, insulating_sides)
            if full_contact_sides:
                predictor = apply_normal_contact_boundary(predictor, full_contact_sides)
            if masked_contact_sides:
                predictor = apply_normal_contact_mask(predictor, contact_mask)
            predicted_rhs = (
                covariant_laplacian(predictor, ax, ay, dx, dy, active_mask)
                + epsilon*predictor - np.abs(predictor)**2*predictor
            )
            predicted_rhs *= np.exp(1j*scalar_potential*sub_dt)
            psi = generalized_tdgl_update(
                psi, 0.5*(rhs+predicted_rhs), scalar_potential, sub_dt,
                tdgl_model.parameters.u, tdgl_model.parameters.gamma,
            )
        else:
            psi = predictor
        psi[~active_mask] = 0.0
        if not np.all(np.isfinite(psi)):
            raise RuntimeError('TDGL produced a nonfinite order parameter.')
        elapsed += sub_dt
        count += 1
    psi = apply_insulating_boundary(psi, ax, ay, dx, dy, insulating_sides)
    if full_contact_sides:
        psi = apply_normal_contact_boundary(psi, full_contact_sides)
    if masked_contact_sides:
        psi = apply_normal_contact_mask(psi, contact_mask)
    jx, jy = tdgl_model.supercurrent_density(
        psi, ax, ay, dx, dy, active_mask
    )
    fields.psi = psi
    fields.supercurrent_density_x = scales.supercurrent_density_to_physical(jx)
    fields.supercurrent_density_y = scales.supercurrent_density_to_physical(jy)
    fields.tdgl_solver_substeps = count
    return fields
