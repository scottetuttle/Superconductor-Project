"""Configuration-sweep diagnostics for the coupled SHS/TDGL simulator.

This tool runs independent production-solver simulations, records physics and
solver metrics, detects signed phase vortices, and creates CSV/JSON data,
kymographs, comparison plots, snapshots, and optional GIF animations.

Run from the repository root:
    python tools/tdgl_diagnostics.py
    python tools/tdgl_diagnostics.py --config tools/config/tdgl_diagnostics.json
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from copy import deepcopy
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from shs.config.builder import build_simulation
from shs.physics.diagnostics import evaluate_physics_diagnostics
from shs.physics.fluxoid import rectangular_fluxoid
from shs.physics.electromagnetics import (
    magnetic_field_from_sheet_current,
    perpendicular_magnetic_field,
    uniform_perpendicular_vector_potential,
)
from shs.solvers.coupled_solver import coupled_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.solvers.coupled_solver import build_tdgl_model
from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM, VACUUM_PERMEABILITY
from shs.optics.moving_laser import apply_laser_heat_source
from shs.physics.pinning import pinning_suppression
from shs.physics.josephson import junction_observables
from shs.utils.output import reserve_output_directory
from shs.utils.progress import PercentProgress


DEFAULT_CONFIG = Path(__file__).resolve().parent / "config" / "tdgl_diagnostics.json"


@dataclass(frozen=True)
class VortexResult:
    winding: np.ndarray
    positive_count: int
    negative_count: int

    @property
    def total_count(self):
        return self.positive_count + self.negative_count


def _wrapped(angle):
    return np.angle(np.exp(1j * angle))


def relative_phase_change(frame, reference, amplitude_floor=1e-8):
    """Return spatial phase change with unobservable global rotation removed."""
    delta = _wrapped(frame["phase"] - reference["phase"])
    weights = np.minimum(frame["amplitude"], reference["amplitude"])
    valid = weights > amplitude_floor
    if np.any(valid):
        global_rotation = np.angle(np.sum(weights[valid] * np.exp(1j * delta[valid])))
        delta = _wrapped(delta - global_rotation)
    return np.where(valid, delta, np.nan)


def detect_vortices(psi, ax=None, ay=None, dx=1.0, dy=1.0, amplitude_floor=1e-8):
    """Return signed integer winding on each grid plaquette.

    ``ax`` and ``ay`` are dimensionless vector potentials. Restoring their
    plaquette circulation to the gauge-covariant link phases recovers the
    topological phase winding while retaining gauge invariance.
    """
    psi = np.asarray(psi, dtype=complex)
    if psi.ndim != 2 or min(psi.shape) < 2:
        raise ValueError("Vortex detection requires a two-dimensional 2x2 grid.")
    ax = np.zeros_like(psi.real) if ax is None else np.asarray(ax, dtype=float)
    ay = np.zeros_like(psi.real) if ay is None else np.asarray(ay, dtype=float)
    if ax.shape != psi.shape or ay.shape != psi.shape:
        raise ValueError("psi, ax, and ay must have identical shapes.")

    phase = np.angle(psi)
    qx = _wrapped(phase[:, 1:] - phase[:, :-1] - ax[:, :-1] * dx)
    qy = _wrapped(phase[1:, :] - phase[:-1, :] - ay[:-1, :] * dy)
    gauge_circulation = qx[:-1, :] + qy[:, 1:] - qx[1:, :] - qy[:, :-1]
    flux = (
        ax[:-1, :-1] * dx
        + ay[:-1, 1:] * dy
        - ax[1:, :-1] * dx
        - ay[:-1, :-1] * dy
    )
    winding = np.rint((gauge_circulation + flux) / (2.0 * np.pi)).astype(int)
    valid = np.minimum.reduce(
        [
            np.abs(psi[:-1, :-1]), np.abs(psi[:-1, 1:]),
            np.abs(psi[1:, 1:]), np.abs(psi[1:, :-1]),
        ]
    ) > amplitude_floor
    winding = np.where(valid, winding, 0)
    return VortexResult(
        winding=winding,
        positive_count=int(np.sum(winding > 0)),
        negative_count=int(np.sum(winding < 0)),
    )


def track_vortices(frames, rows, max_displacement_cells=12.0):
    """Assign persistent IDs using sign-preserving nearest-neighbor matching."""
    active = {}
    next_id = 1
    records = []
    for frame_index, (frame, row) in enumerate(zip(frames, rows)):
        winding = frame["winding"]
        detections = [
            (int(winding[y, x]), float(x) + 0.5, float(y) + 0.5)
            for y, x in np.argwhere(winding != 0)
        ]
        assigned_ids, used = set(), set()
        candidates = []
        for detection_index, (charge, x, y) in enumerate(detections):
            for track_id, previous in active.items():
                if previous[0] == charge:
                    distance = np.hypot(x - previous[1], y - previous[2])
                    if distance <= max_displacement_cells:
                        candidates.append((distance, detection_index, track_id))
        matches = {}
        for _, detection_index, track_id in sorted(candidates):
            if detection_index not in matches and track_id not in assigned_ids:
                matches[detection_index] = track_id
                assigned_ids.add(track_id)
        updated = {}
        for detection_index, (charge, x, y) in enumerate(detections):
            track_id = matches.get(detection_index)
            if track_id is None:
                track_id, next_id = next_id, next_id + 1
            updated[track_id] = (charge, x, y)
            used.add(track_id)
            records.append({
                "track_id": track_id, "charge": charge,
                "frame": frame_index, "step": row["step"], "time_s": row["time_s"],
                "x_cell": x, "y_cell": y,
                "x_m": x * frame["mesh_dx_m"], "y_m": y * frame["mesh_dy_m"],
            })
        active = updated
    return records


def _set_path(root, dotted_path, value):
    target = root
    parts = dotted_path.split(".")
    for part in parts[:-1]:
        target = getattr(target, part)
    if not hasattr(target, parts[-1]):
        raise ValueError(f"Unknown sweep setting: {dotted_path}")
    setattr(target, parts[-1], value)


def build_cases(sweeps):
    """Expand a mapping of dotted configuration paths into a Cartesian sweep."""
    if not sweeps:
        return [("baseline", {})]
    paths = list(sweeps)
    values = []
    for path in paths:
        choices = sweeps[path]
        if not isinstance(choices, list) or not choices:
            raise ValueError(f"Sweep {path} must contain a nonempty list.")
        values.append(choices)
    cases = []
    for index, combination in enumerate(itertools.product(*values), 1):
        settings = dict(zip(paths, combination))
        label = f"case_{index:03d}_" + "_".join(
            f"{path.split('.')[-1]}={value:g}" if isinstance(value, (int, float))
            else f"{path.split('.')[-1]}={value}"
            for path, value in settings.items()
        )
        cases.append((label.replace("/", "_"), settings))
    return cases


def assess_model_regime(simulation, validity):
    """Report whether mesh, temperature, and thickness fit the 2D TDGL use case."""
    material = simulation.material_map.materials[0]
    reduced_field = simulation.fields.temperature / material.Tc
    reduced_temperature = float(np.mean(reduced_field))
    metrics = {
        "reduced_temperature": reduced_temperature,
        "minimum_reduced_temperature": float(np.min(reduced_field)),
        "maximum_reduced_temperature": float(np.max(reduced_field)),
        "distance_from_Tc": float(np.max(np.abs(1.0 - reduced_field))),
        "thickness_over_xi": simulation.geometry.film.thickness / material.coherence_length,
        "thickness_over_lambda": simulation.geometry.film.thickness / material.penetration_depth,
        "dx_over_xi": simulation.mesh.dx / material.coherence_length,
        "dy_over_xi": simulation.mesh.dy / material.coherence_length,
    }
    warnings = []
    if metrics["distance_from_Tc"] > validity["maximum_distance_from_Tc"]:
        warnings.append("Temperature is outside the configured near-Tc TDGL range.")
    if metrics["thickness_over_xi"] >= validity["maximum_thickness_over_xi"]:
        warnings.append("Film is not thin compared with the coherence length.")
    if metrics["thickness_over_lambda"] >= validity["maximum_thickness_over_lambda"]:
        warnings.append("Film is not thin compared with the penetration depth.")
    if max(metrics["dx_over_xi"], metrics["dy_over_xi"]) > validity["maximum_spacing_over_xi"]:
        warnings.append("Mesh is too coarse relative to the coherence length.")
    return {"metrics": metrics, "warnings": warnings, "passes": not warnings}


def _seed_initial_state(simulation, config):
    fields, mesh = simulation.fields, simulation.mesh
    yy, xx = np.indices(fields.psi.shape)
    hotspot = config.get("temperature_hotspot", {})
    if hotspot.get("enabled", False):
        cx = hotspot.get("x_fraction", 0.5) * (mesh.nx - 1)
        cy = hotspot.get("y_fraction", 0.5) * (mesh.ny - 1)
        radius = float(hotspot.get("radius_cells", 3.0))
        fields.temperature += float(hotspot.get("peak_delta_K", 1.0)) * np.exp(
            -((xx - cx) ** 2 + (yy - cy) ** 2) / (2.0 * radius**2)
        )

    for vortex in config.get("vortices", []):
        cx = float(vortex.get("x_fraction", 0.5)) * (mesh.nx - 1)
        cy = float(vortex.get("y_fraction", 0.5)) * (mesh.ny - 1)
        charge = int(vortex.get("charge", 1))
        radius = max(float(vortex.get("core_radius_cells", 2.0)), 1e-12)
        distance = np.hypot(xx - cx, yy - cy)
        phase = charge * np.arctan2(yy - cy, xx - cx)
        fields.psi *= np.tanh(distance / radius) * np.exp(1j * phase)


def apply_uniform_magnetic_field(simulation, magnetic_field_z, gauge="symmetric"):
    """Install applied Bz and transform psi from the symmetric reference gauge."""
    ax, ay = uniform_perpendicular_vector_potential(
        simulation.mesh, magnetic_field_z, gauge
    )
    simulation.fields.vector_potential_x = ax
    simulation.fields.vector_potential_y = ay
    simulation.fields.applied_vector_potential_x = ax.copy()
    simulation.fields.applied_vector_potential_y = ay.copy()
    simulation.fields.induced_vector_potential_x.fill(0.0)
    simulation.fields.induced_vector_potential_y.fill(0.0)
    simulation.fields.magnetic_field_z = perpendicular_magnetic_field(
        ax, ay, simulation.mesh.dx, simulation.mesh.dy
    )
    if gauge != "symmetric":
        xx, yy = np.meshgrid(
            simulation.mesh.x - np.mean(simulation.mesh.x),
            simulation.mesh.y - np.mean(simulation.mesh.y),
        )
        sign = -1.0 if gauge == "landau_x" else 1.0
        gauge_function = sign * 0.5 * magnetic_field_z * xx * yy
        simulation.fields.psi *= np.exp(
            1j * 2.0 * np.pi * gauge_function / SUPERCONDUCTING_FLUX_QUANTUM
        )


def cross_section_current(simulation, column=None):
    """Integrate total x-current density through a vertical film section."""
    column = simulation.mesh.nx // 2 if column is None else int(column)
    column = min(max(column, 0), simulation.mesh.nx - 2)
    profile = simulation.fields.current_density_x[:, column]
    return float(
        np.trapezoid(profile, dx=simulation.mesh.dy)
        * float(np.mean(simulation.material_map.thickness))
    )


def current_uniformity(simulation, column=None):
    column = simulation.mesh.nx // 2 if column is None else int(column)
    column = min(max(column, 0), simulation.mesh.nx - 2)
    profile = simulation.fields.current_density_x[:, column]
    return float(np.std(profile) / max(np.sqrt(np.mean(profile**2)), np.finfo(float).tiny))


def cross_section_current_profile(simulation):
    """Integrated x-directed current through every inter-column section."""
    return (
        np.trapezoid(simulation.fields.current_density_x[:, :-1],
                     dx=simulation.mesh.dy, axis=0)
        * float(np.mean(simulation.material_map.thickness))
    )


def current_controlled_step(simulation, dt, tdgl_model, target_current,
                            relative_tolerance, max_iterations, voltage_limit,
                            time_s=0.0):
    """Tune terminal voltage until an accepted step carries target current.

    Each voltage trial restarts from the same accepted physical state. This is
    a global current controller around the Dirichlet transport solve; it does
    not silently advance physical time during control iterations.
    """
    if not np.isfinite(target_current):
        raise ValueError("Target current must be finite.")
    accepted = deepcopy(simulation.fields)
    right = simulation.config.electrical.voltage_right
    voltage = simulation.config.electrical.voltage_left - right
    if voltage == 0:
        voltage = np.copysign(min(1e-9, voltage_limit), target_current or 1.0)
    last = None
    for control_iteration in range(1, max_iterations + 1):
        simulation.fields = deepcopy(accepted)
        simulation.config.electrical.voltage_left = right + voltage
        _, result = coupled_step(
            simulation, dt=dt, tdgl_model=tdgl_model, time_s=time_s
        )
        if not result.converged:
            continue
        measured = cross_section_current(simulation)
        error = measured - target_current
        if abs(error) <= relative_tolerance * max(abs(target_current), 1e-30):
            return result, control_iteration, measured, voltage
        if last is not None and measured != last[1]:
            next_voltage = voltage - error * (voltage - last[0]) / (measured - last[1])
        elif measured != 0:
            next_voltage = voltage * target_current / measured
        else:
            next_voltage = 2.0 * voltage
        last = (voltage, measured)
        voltage = float(np.clip(next_voltage, -voltage_limit, voltage_limit))
    simulation.fields = accepted
    raise RuntimeError(
        f"Current controller did not reach {target_current:.6g} A after "
        f"{max_iterations} voltage trials."
    )


def _new_output_directory(base):
    return reserve_output_directory(base)


def _field_snapshot(simulation, time_value, step, scales, vortex_floor,
                    magnetic_config=None, compute_self_field=False,
                    fluxoid_config=None):
    fields, mesh = simulation.fields, simulation.mesh
    vortices = detect_vortices(
        fields.psi,
        scales.vector_potential_to_dimensionless(fields.vector_potential_x),
        scales.vector_potential_to_dimensionless(fields.vector_potential_y),
        simulation.mesh.dx / scales.xi,
        simulation.mesh.dy / scales.xi,
        vortex_floor,
    )
    physics = evaluate_physics_diagnostics(simulation).to_dict()
    calculated_bz = perpendicular_magnetic_field(
        fields.vector_potential_x,
        fields.vector_potential_y,
        mesh.dx,
        mesh.dy,
    )
    applied_bz = perpendicular_magnetic_field(
        fields.applied_vector_potential_x,
        fields.applied_vector_potential_y,
        mesh.dx, mesh.dy,
    )
    induced_bz = perpendicular_magnetic_field(
        fields.induced_vector_potential_x,
        fields.induced_vector_potential_y,
        mesh.dx, mesh.dy,
    )
    total_bz = applied_bz + induced_bz
    field_error = calculated_bz - applied_bz
    measured_current = cross_section_current(simulation)
    center_column = min(simulation.mesh.nx // 2, simulation.mesh.nx - 2)
    thickness = float(np.mean(simulation.material_map.thickness))
    supercurrent_A = float(np.trapezoid(
        fields.supercurrent_density_x[:, center_column], dx=mesh.dy) * thickness)
    normal_current_A = float(np.trapezoid(
        fields.normal_current_density_x[:, center_column], dx=mesh.dy) * thickness)
    current_sections = cross_section_current_profile(simulation)
    longitudinal_imbalance = float(
        np.std(current_sections)
        / max(np.sqrt(np.mean(current_sections**2)), np.finfo(float).tiny)
    )
    cell_area = mesh.dx * mesh.dy
    applied_flux = float(np.sum(applied_bz) * cell_area)
    laser_x = getattr(fields, "laser_position_x_m", None)
    laser_y = getattr(fields, "laser_position_y_m", None)
    laser_x = np.nan if laser_x is None else float(laser_x)
    laser_y = np.nan if laser_y is None else float(laser_y)
    vortex_y, vortex_x = np.nonzero(vortices.winding)
    if vortex_x.size and np.isfinite(laser_x) and np.isfinite(laser_y):
        vortex_x_m = (vortex_x + 0.5) * mesh.dx
        vortex_y_m = (vortex_y + 0.5) * mesh.dy
        nearest_laser_distance = float(np.min(np.hypot(
            vortex_x_m - laser_x, vortex_y_m - laser_y
        )))
    else:
        nearest_laser_distance = float("nan")
    row = {
        "step": step,
        "time_s": time_value,
        "mean_temperature_K": float(np.mean(fields.temperature)),
        "maximum_temperature_K": float(np.max(fields.temperature)),
        "mean_phonon_temperature_K": float(np.mean(fields.phonon_temperature)) if fields.phonon_temperature is not None else float(np.mean(fields.temperature)),
        "maximum_phonon_temperature_K": float(np.max(fields.phonon_temperature)) if fields.phonon_temperature is not None else float(np.max(fields.temperature)),
        "maximum_electron_phonon_difference_K": float(np.max(np.abs(fields.temperature-fields.phonon_temperature))) if fields.phonon_temperature is not None else 0.0,
        "mean_order_parameter": float(np.mean(np.abs(fields.psi))),
        "minimum_order_parameter": float(np.min(np.abs(fields.psi))),
        "mean_voltage_V": float(np.mean(fields.voltage)),
        "positive_vortices": vortices.positive_count,
        "negative_vortices": vortices.negative_count,
        "total_vortices": vortices.total_count,
        "coupling_iterations": 0,
        "coupling_residual": 0.0,
        "tdgl_internal_substeps": int(getattr(fields, "tdgl_solver_substeps", 0)),
        "thermal_internal_substeps": int(getattr(fields, "thermal_solver_substeps", 0)),
        "electrical_iterations": int(getattr(fields, "electrical_solver_iterations", 0)),
        "electrical_residual_V": float(getattr(fields, "electrical_solver_residual", 0.0)),
        "transport_current_A": measured_current,
        "cross_section_supercurrent_A": supercurrent_A,
        "cross_section_normal_current_A": normal_current_A,
        "current_profile_nonuniformity": current_uniformity(simulation),
        "longitudinal_current_imbalance": longitudinal_imbalance,
        "maximum_current_density_A_per_m2": float(np.max(np.hypot(
            fields.current_density_x, fields.current_density_y
        ))),
        "laser_x_m": laser_x,
        "laser_y_m": laser_y,
        "deposited_laser_power_W": float(np.sum(
            (np.zeros_like(fields.temperature) if fields.laser_heat_source is None
             else fields.laser_heat_source)
            * simulation.material_map.thickness * mesh.dx * mesh.dy
            * simulation.region_map.active_mask
        )),
        "nearest_vortex_laser_distance_m": nearest_laser_distance,
        "mean_applied_magnetic_field_T": float(np.mean(applied_bz)),
        "mean_total_magnetic_field_T": float(np.mean(total_bz)),
        "center_total_magnetic_field_T": float(
            total_bz[mesh.ny // 2, mesh.nx // 2]
        ),
        "maximum_induced_magnetic_field_T": float(np.max(np.abs(induced_bz))),
        "maximum_magnetic_curl_error_T": float(np.max(np.abs(field_error))),
        "applied_flux_Wb": applied_flux,
        "applied_flux_quanta": applied_flux / SUPERCONDUCTING_FLUX_QUANTUM,
        "applied_magnetic_energy_J": float(
            np.sum(applied_bz**2) * cell_area
            * float(np.mean(simulation.material_map.thickness))
            / (2.0 * VACUUM_PERMEABILITY)
        ),
        **physics,
    }
    if getattr(simulation.config.josephson, "enabled", False):
        row.update(junction_observables(simulation))
    fluxoid_config = {} if fluxoid_config is None else fluxoid_config
    inset = float(fluxoid_config.get("contour_inset_fraction", 0.25))
    y0 = max(0, min(mesh.ny - 2, int(round(inset * (mesh.ny - 1)))))
    y1 = min(mesh.ny - 1, max(y0 + 1, int(round((1.0 - inset) * (mesh.ny - 1)))))
    x0 = max(0, min(mesh.nx - 2, int(round(inset * (mesh.nx - 1)))))
    x1 = min(mesh.nx - 1, max(x0 + 1, int(round((1.0 - inset) * (mesh.nx - 1)))))
    try:
        fluxoid = rectangular_fluxoid(
            simulation, (y0, y1, x0, x1), amplitude_floor=vortex_floor
        )
        row.update({
            "fluxoid_flux_part_Wb": fluxoid.flux_part_Wb,
            "fluxoid_supercurrent_part_Wb": fluxoid.supercurrent_part_Wb,
            "london_fluxoid_quanta": fluxoid.london_fluxoid_quanta,
            "topological_fluxoid_quanta": fluxoid.topological_fluxoid_quanta,
            "fluxoid_quantization_error": (
                fluxoid.london_fluxoid_quanta - fluxoid.topological_fluxoid_quanta
            ),
        })
    except ValueError:
        row.update({
            "fluxoid_flux_part_Wb": float("nan"),
            "fluxoid_supercurrent_part_Wb": float("nan"),
            "london_fluxoid_quanta": float("nan"),
            "topological_fluxoid_quanta": float("nan"),
            "fluxoid_quantization_error": float("nan"),
        })
    frame = {
        "temperature": fields.temperature.copy(),
        "phonon_temperature": (fields.phonon_temperature.copy() if fields.phonon_temperature is not None else fields.temperature.copy()),
        "amplitude": np.abs(fields.psi).copy(),
        "phase": np.angle(fields.psi).copy(),
        "voltage": fields.voltage.copy(),
        "winding": vortices.winding.copy(),
        "current_x": fields.current_density_x.copy(),
        "current_y": fields.current_density_y.copy(),
        "supercurrent_x": fields.supercurrent_density_x.copy(),
        "supercurrent_y": fields.supercurrent_density_y.copy(),
        "normal_current_x": fields.normal_current_density_x.copy(),
        "normal_current_y": fields.normal_current_density_y.copy(),
        "applied_magnetic_field_z": applied_bz.copy(),
        "induced_magnetic_field_z": induced_bz.copy(),
        "magnetic_field_z": total_bz.copy(),
        "vector_potential_x": fields.vector_potential_x.copy(),
        "vector_potential_y": fields.vector_potential_y.copy(),
        "laser_heat_source": (
            np.zeros_like(fields.temperature) if fields.laser_heat_source is None
            else fields.laser_heat_source.copy()
        ),
        "laser_x_m": laser_x,
        "laser_y_m": laser_y,
        "pinning_suppression": pinning_suppression(simulation).copy(),
        "mesh_dx_m": mesh.dx,
        "mesh_dy_m": mesh.dy,
    }
    junction = getattr(simulation.config, "josephson", None)
    if junction is not None and junction.enabled:
        frame["junction_x_bounds_cells"] = (
            (junction.center_x_m - 0.5 * junction.width_m) / mesh.dx,
            (junction.center_x_m + 0.5 * junction.width_m) / mesh.dx,
        )
    if (compute_self_field and magnetic_config
            and magnetic_config.get("calculate_self_field", False)):
        bx, by, bz = magnetic_field_from_sheet_current(
            mesh,
            fields.current_density_x,
            fields.current_density_y,
            float(np.mean(simulation.material_map.thickness)),
            float(magnetic_config["observation_height_m"]),
            int(magnetic_config["source_stride"]),
        )
        frame["self_field_x"] = bx
        frame["self_field_y"] = by
        frame["self_field_z"] = bz
        row["maximum_self_field_T"] = float(np.max(np.sqrt(bx**2 + by**2 + bz**2)))
    else:
        frame["self_field_x"] = np.zeros_like(applied_bz)
        frame["self_field_y"] = np.zeros_like(applied_bz)
        frame["self_field_z"] = np.zeros_like(applied_bz)
        row["maximum_self_field_T"] = 0.0
    return row, frame


def _write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def _plot_timeseries(rows, path, dpi):
    time_fs = np.array([row["time_s"] for row in rows]) * 1e15
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axes[0, 0].plot(time_fs, [r["mean_order_parameter"] for r in rows], label="mean")
    axes[0, 0].plot(time_fs, [r["minimum_order_parameter"] for r in rows], label="minimum")
    axes[0, 0].set_ylabel("|psi|"); axes[0, 0].legend()
    axes[0, 1].plot(time_fs, [r["maximum_temperature_K"] for r in rows])
    axes[0, 1].set_ylabel("Maximum temperature (K)")
    axes[1, 0].step(time_fs, [r["positive_vortices"] for r in rows], where="post", label="+1")
    axes[1, 0].step(time_fs, [r["negative_vortices"] for r in rows], where="post", label="-1")
    axes[1, 0].set_ylabel("Detected vortices"); axes[1, 0].legend()
    axes[1, 1].semilogy(time_fs, np.maximum([r["current_continuity_relative"] for r in rows], 1e-30))
    axes[1, 1].set_ylabel("Relative current-continuity defect")
    for axis in axes[-1]: axis.set_xlabel("Time (fs)")
    fig.savefig(path, dpi=dpi); plt.close(fig)


def _plot_kymographs(frames, rows, centerline, path, dpi):
    if centerline == "horizontal":
        extract = lambda array: array[array.shape[0] // 2, :]
        spatial_label = "x cell"
    else:
        extract = lambda array: array[:, array.shape[1] // 2]
        spatial_label = "y cell"
    names = [("temperature", "Temperature (K)"), ("amplitude", "|psi|"), ("phase", "Phase (rad)")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    extent = [0, len(extract(frames[0]["amplitude"])) - 1,
              rows[0]["time_s"] * 1e15, rows[-1]["time_s"] * 1e15]
    for axis, (name, title) in zip(axes, names):
        data = np.stack([extract(frame[name]) for frame in frames])
        image = axis.imshow(data, origin="lower", aspect="auto", extent=extent,
                            cmap="twilight" if name == "phase" else "viridis")
        axis.set_title(title); axis.set_xlabel(spatial_label)
        fig.colorbar(image, ax=axis)
    axes[0].set_ylabel("Time (fs)")
    fig.savefig(path, dpi=dpi); plt.close(fig)


def _plot_final_state(frame, path, dpi):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    specs = [("amplitude", "|psi|", "viridis"), ("phase", "Phase", "twilight"),
             ("temperature", "Temperature (K)", "inferno"), ("voltage", "Voltage (V)", "coolwarm")]
    for axis, (name, title, cmap) in zip(axes.flat, specs):
        image = axis.imshow(frame[name], origin="lower", cmap=cmap)
        axis.set_title(title); fig.colorbar(image, ax=axis)
        if "junction_x_bounds_cells" in frame:
            axis.axvspan(*frame["junction_x_bounds_cells"], color="cyan", alpha=0.16)
            for boundary in frame["junction_x_bounds_cells"]:
                axis.axvline(boundary, color="cyan", linewidth=0.8, linestyle="--")
    vortex_y, vortex_x = np.nonzero(frame["winding"])
    if vortex_x.size:
        axes[0, 1].scatter(vortex_x + 0.5, vortex_y + 0.5,
                           c=np.where(frame["winding"][vortex_y, vortex_x] > 0, "white", "black"),
                           s=25, marker="x")
    fig.savefig(path, dpi=dpi); plt.close(fig)


def _plot_current_and_field(frame, path, dpi, vector_stride, self_consistent=False):
    total_magnitude = np.hypot(frame["current_x"], frame["current_y"])
    super_magnitude = np.hypot(frame["supercurrent_x"], frame["supercurrent_y"])
    normal_magnitude = np.hypot(frame["normal_current_x"], frame["normal_current_y"])
    total_bz = frame["magnetic_field_z"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    for axis, magnitude, title in (
        (axes[0, 0], total_magnitude, "Total current density (A/m²)"),
        (axes[0, 1], super_magnitude, "Supercurrent density (A/m²)"),
        (axes[1, 0], normal_magnitude, "Normal current density (A/m²)"),
    ):
        image = axis.imshow(magnitude, origin="lower", cmap="magma")
        fig.colorbar(image, ax=axis)
        axis.set_title(title)
        if "junction_x_bounds_cells" in frame:
            axis.axvspan(*frame["junction_x_bounds_cells"], color="cyan", alpha=0.16)
    stride = max(1, int(vector_stride))
    yy, xx = np.indices(total_magnitude.shape)
    axes[0, 0].streamplot(
        xx[0], yy[:, 0], frame["current_x"], frame["current_y"],
        color="white", density=0.8, linewidth=0.5,
    )
    field_image = axes[1, 1].imshow(total_bz, origin="lower", cmap="coolwarm")
    fig.colorbar(field_image, ax=axes[1, 1], label="Bz (T)")
    axes[1, 1].set_title(
        "Self-consistent total Bz" if self_consistent
        else "Applied + diagnostic self-field"
    )
    axes[1, 1].streamplot(
        xx[0], yy[:, 0], frame["vector_potential_x"],
        frame["vector_potential_y"], color="white", density=0.7,
        linewidth=0.5,
    )
    sampled_bx = frame["self_field_x"][::stride, ::stride]
    sampled_by = frame["self_field_y"][::stride, ::stride]
    if np.any(np.hypot(sampled_bx, sampled_by) > 0):
        axes[1, 1].quiver(
            xx[::stride, ::stride], yy[::stride, ::stride],
            sampled_bx, sampled_by, color="black", pivot="mid",
        )
    fig.savefig(path, dpi=dpi); plt.close(fig)


def _plot_junction_profile(frame, path, dpi):
    if "junction_x_bounds_cells" not in frame:
        return
    x = np.arange(frame["amplitude"].shape[1])
    amplitude = np.mean(frame["amplitude"], axis=0)
    supercurrent = np.mean(frame["supercurrent_x"], axis=0)
    normal = np.mean(frame["normal_current_x"], axis=0)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    axes[0].plot(x, amplitude, color="navy"); axes[0].set_ylabel("Mean |psi|")
    axes[1].plot(x, supercurrent, label="supercurrent")
    axes[1].plot(x, normal, label="normal current"); axes[1].legend()
    axes[1].set(xlabel="x cell", ylabel="Mean Jx (A/m²)")
    for axis in axes:
        axis.axvspan(*frame["junction_x_bounds_cells"], color="cyan", alpha=0.2,
                    label="configured weak link")
        axis.grid(alpha=.25)
    axes[0].set_title("Cross-device weak-link profile")
    fig.savefig(path, dpi=dpi); plt.close(fig)


def _plot_vortex_trajectories(frames, rows, records, path, dpi):
    """Plot signed vortex tracks together with the commanded laser path."""
    fig, axis = plt.subplots(figsize=(8, 7), constrained_layout=True)
    axis.imshow(frames[-1]["amplitude"], origin="lower", cmap="gray", alpha=0.75)
    track_ids = sorted({record["track_id"] for record in records})
    for track_id in track_ids:
        track = [record for record in records if record["track_id"] == track_id]
        charge = track[0]["charge"]
        color = "#ff4d4d" if charge > 0 else "#4da6ff"
        axis.plot([point["x_cell"] for point in track],
                  [point["y_cell"] for point in track], "-o", color=color,
                  linewidth=1.8, markersize=3,
                  label=f"{'+' if charge > 0 else ''}{charge}, track {track_id}")
        axis.scatter(track[-1]["x_cell"], track[-1]["y_cell"], color=color,
                     marker="x", s=55)
    laser_x = [frame["laser_x_m"] / frame["mesh_dx_m"] for frame in frames
               if np.isfinite(frame["laser_x_m"])]
    laser_y = [frame["laser_y_m"] / frame["mesh_dy_m"] for frame in frames
               if np.isfinite(frame["laser_y_m"])]
    if laser_x:
        axis.plot(laser_x, laser_y, "--", color="#ffd43b", linewidth=2.2,
                  label="laser path")
        axis.scatter(laser_x[-1], laser_y[-1], s=100, facecolors="none",
                     edgecolors="#ffd43b", linewidths=2, label="laser final")
    axis.set(title="Vortex trajectories and laser path", xlabel="x cell", ylabel="y cell")
    handles, labels = axis.get_legend_handles_labels()
    if handles:
        axis.legend(loc="upper left", fontsize=8, framealpha=0.85)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def _plot_laser_vortex_distance(rows, path, dpi):
    """Plot the observable used to distinguish capture from incidental motion."""
    time_fs = np.asarray([row["time_s"] for row in rows]) * 1e15
    distance_nm = np.asarray([
        row["nearest_vortex_laser_distance_m"] for row in rows
    ]) * 1e9
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                             constrained_layout=True)
    axes[0].plot(time_fs, distance_nm, "o-", color="#6f42c1")
    axes[0].set_ylabel("Nearest vortex-laser distance (nm)")
    axes[0].grid(alpha=0.25)
    axes[1].plot(time_fs, [row["maximum_temperature_K"] for row in rows],
                 color="#d94801", label="maximum temperature")
    axes[1].set(xlabel="Time (fs)", ylabel="Temperature (K)")
    axes[1].grid(alpha=0.25)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def _save_gif(frames, rows, path, fps, dpi):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    amp_image = axes[0].imshow(frames[0]["amplitude"], origin="lower", cmap="viridis", vmin=0)
    phase_image = axes[1].imshow(frames[0]["phase"], origin="lower", cmap="twilight", vmin=-np.pi, vmax=np.pi)
    initial_delta = relative_phase_change(frames[0], frames[0])
    delta_image = axes[2].imshow(
        initial_delta, origin="lower", cmap="coolwarm", vmin=-1e-6, vmax=1e-6
    )
    axes[0].set_title("|psi|")
    axes[1].set_title("Absolute phase")
    axes[2].set_title("Spatial phase change")
    title = fig.suptitle("")

    def update(index):
        amp_image.set_data(frames[index]["amplitude"])
        phase_image.set_data(frames[index]["phase"])
        phase_delta = relative_phase_change(frames[index], frames[0])
        delta_image.set_data(phase_delta)
        finite_delta = np.abs(phase_delta[np.isfinite(phase_delta)])
        delta_limit = max(float(np.max(finite_delta, initial=0.0)), 1e-6)
        delta_image.set_clim(-delta_limit, delta_limit)
        amp_image.set_clim(0, max(float(np.max(frames[index]["amplitude"])), 1e-12))
        title.set_text(
            f"t = {rows[index]['time_s'] * 1e15:.3f} fs; "
            f"vortices = {rows[index]['total_vortices']}; "
            f"max spatial dphase = {delta_limit:.3g} rad"
        )
        return amp_image, phase_image, delta_image, title

    movie = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
    movie.save(path, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)


def run_case(label, settings, config, output_directory):
    simulation = build_simulation(config["simulation_config"])
    diagnostic_settings = {
        path.removeprefix("diagnostics."): value
        for path, value in settings.items() if path.startswith("diagnostics.")
    }
    for path, value in settings.items():
        if not path.startswith("diagnostics."):
            _set_path(simulation.config, path, value)
    simulation.config.validate()
    magnetic_config = deepcopy(config["magnetic_field"])
    magnetic_config.update({
        key: value for key, value in diagnostic_settings.items()
        if key in magnetic_config
    })
    applied_field = float(magnetic_config["applied_Bz_T"])
    apply_uniform_magnetic_field(simulation, applied_field, magnetic_config["gauge"])
    _seed_initial_state(simulation, config.get("initial_conditions", {}))
    regime = assess_model_regime(simulation, config["validity"])
    tdgl_model = build_tdgl_model(simulation)
    scales = tdgl_scales(simulation, tdgl_model)
    run = config["run"]
    dt = float(run.get("dt", simulation.config.dt))
    duration_multiplier = float(run.get("duration_multiplier", 1.0))
    steps = max(1, int(round(
        simulation.config.duration * duration_multiplier / dt
    )))
    sample_every = int(run["sample_every"])
    resolved_duration = steps * dt
    if not np.isclose(
        resolved_duration, simulation.config.duration,
        rtol=64 * np.finfo(float).eps, atol=0.0,
    ):
        print(
            "  configuration override: diagnostics run duration "
            f"{resolved_duration:.6g} s replaces simulation duration "
            f"{simulation.config.duration:.6g} s",
            flush=True,
        )
    case_directory = output_directory / label
    case_directory.mkdir(parents=True)
    rows, frames = [], []
    apply_laser_heat_source(simulation, 0.0)
    row, frame = _field_snapshot(
        simulation, 0.0, 0, scales,
        config["vortex_detection"]["amplitude_floor"], magnetic_config, False,
        config.get("fluxoid", {}),
    )
    rows.append(row); frames.append(frame)
    start = perf_counter()
    converged = True
    coupling_iterations = 0
    voltage_control_iterations = 0
    transport = config.get("transport", {})
    drive_mode = transport.get(
        "drive_mode", simulation.config.electrical.drive_mode
    )
    if "target_current_A" in diagnostic_settings:
        drive_mode = "current"
    target_current = float(
        diagnostic_settings.get("target_current_A", simulation.config.current)
    )
    if drive_mode == "current":
        simulation.config.electrical.drive_mode = "current"
        simulation.config.current = target_current
    else:
        simulation.config.electrical.drive_mode = "voltage"
    progress_reporter = PercentProgress(steps, label=label)
    for step in range(1, steps + 1):
        if drive_mode == "voltage_search":
            result, control_work, _, _ = current_controlled_step(
                simulation, dt, tdgl_model, target_current,
                float(transport["relative_tolerance"]),
                int(transport["max_voltage_iterations"]),
                float(transport["voltage_limit_V"]),
                time_s=(step - 1) * dt,
            )
            voltage_control_iterations += control_work
        else:
            _, result = coupled_step(
                simulation, dt=dt, tdgl_model=tdgl_model,
                time_s=(step - 1) * dt,
            )
        coupling_iterations += result.iterations
        progress_reporter.update(step)
        if not result.converged:
            converged = False
            break
        if step % sample_every == 0 or step == steps:
            row, frame = _field_snapshot(simulation, step * dt, step, scales,
                                         config["vortex_detection"]["amplitude_floor"],
                                         magnetic_config, step == steps,
                                         config.get("fluxoid", {}))
            row["coupling_iterations"] = result.iterations
            row["coupling_residual"] = result.final_residual
            rows.append(row); frames.append(frame)
    runtime = perf_counter() - start
    _write_csv(case_directory / "diagnostics.csv", rows)
    trajectory_records = track_vortices(
        frames, rows,
        float(config.get("vortex_tracking", {}).get("max_displacement_cells", 12.0)),
    )
    _write_csv(case_directory / "vortex_trajectories.csv", trajectory_records)
    output = config["output"]
    dpi = int(output["dpi"])
    if output["plots"]:
        _plot_timeseries(rows, case_directory / "timeseries.png", dpi)
        _plot_kymographs(frames, rows, output["kymograph_centerline"],
                         case_directory / "kymographs.png", dpi)
        _plot_final_state(frames[-1], case_directory / "final_state.png", dpi)
        _plot_current_and_field(
            frames[-1], case_directory / "current_and_field.png", dpi,
            output["vector_stride"],
            simulation.config.electromagnetic.include_self_field,
        )
        _plot_junction_profile(frames[-1], case_directory / "junction_profile.png", dpi)
        _plot_vortex_trajectories(
            frames, rows, trajectory_records,
            case_directory / "vortex_trajectories.png", dpi,
        )
        _plot_laser_vortex_distance(
            rows, case_directory / "laser_vortex_distance.png", dpi,
        )
    if output["gif"]["enabled"]:
        _save_gif(frames, rows, case_directory / "order_parameter.gif",
                  int(output["gif"]["fps"]), dpi)
    finite_distances = np.asarray([
        row["nearest_vortex_laser_distance_m"] for row in rows
    ], dtype=float)
    finite_distances = finite_distances[np.isfinite(finite_distances)]
    summary = {
        "label": label, "settings": settings, "converged": converged,
        "completed_steps": rows[-1]["step"], "runtime_seconds": runtime,
        "coupling_iterations": coupling_iterations, "final": rows[-1],
        "voltage_control_iterations": voltage_control_iterations,
        "target_current_A": target_current if drive_mode != "voltage" else None,
        "current_control_method": drive_mode,
        "model_regime": regime,
        "configuration_provenance": {
            "simulation_config": config["simulation_config"],
            "sweep_overrides": settings,
            "resolved_dt_s": dt,
            "diagnostics_run_steps": steps,
            "resolved_run_duration_s": resolved_duration,
            "simulation_duration_ignored_by_diagnostics": not np.isclose(
                resolved_duration, simulation.config.duration,
                rtol=64 * np.finfo(float).eps, atol=0.0,
            ),
            "transport_drive_mode": drive_mode,
            "resolved_target_current_A": (
                target_current if drive_mode != "voltage" else None
            ),
            "resolved_applied_Bz_T": applied_field,
        },
        "effective_simulation_config": asdict(simulation.config),
        "laser_vortex_metrics": {
            "minimum_separation_m": (
                float(np.min(finite_distances)) if finite_distances.size else None
            ),
            "final_separation_m": (
                float(finite_distances[-1]) if finite_distances.size else None
            ),
            "capture_claimed": False,
            "capture_note": (
                "Separation is diagnostic only. A capture claim requires a configured "
                "threshold, sustained following, and successful release controls."
            ),
        },
        "physics_models": {
            "tdgl_normalization": simulation.config.tdgl.normalization,
            "temperature_model": simulation.config.tdgl.temperature_model,
            "normal_conductivity_model": simulation.config.electrical.normal_conductivity_model,
            "magnetic_screening": (
                "self_consistent"
                if simulation.config.electromagnetic.include_self_field
                else "diagnostic_only"
            ),
        },
    }
    with (case_directory / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def _plot_comparison(summaries, path, dpi):
    labels = [entry["label"] for entry in summaries]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 2, figsize=(max(10, len(labels) * 1.5), 8), constrained_layout=True)
    metrics = [("maximum_temperature_K", "Final maximum T (K)"),
               ("mean_order_parameter", "Final mean |psi|"),
               ("total_vortices", "Final vortex count")]
    for axis, (key, title) in zip(axes.flat, metrics):
        axis.bar(x, [entry["final"][key] for entry in summaries]); axis.set_ylabel(title)
    axes[1, 1].bar(x, [entry["runtime_seconds"] for entry in summaries])
    axes[1, 1].set_ylabel("Runtime (s)")
    for axis in axes.flat:
        axis.set_xticks(x, labels, rotation=35, ha="right")
    fig.savefig(path, dpi=dpi); plt.close(fig)


def load_config(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    if (float(config["run"].get("duration_multiplier", 1.0)) <= 0
            or int(config["run"]["sample_every"]) < 1):
        raise ValueError(
            "run.duration_multiplier and run.sample_every must be positive."
        )
    if config["output"]["kymograph_centerline"] not in {"horizontal", "vertical"}:
        raise ValueError("kymograph_centerline must be horizontal or vertical.")
    transport = config.get("transport", {})
    if transport.get("drive_mode", "current") not in {"current", "voltage", "voltage_search"}:
        raise ValueError(
            "transport.drive_mode must be current, voltage, or voltage_search."
        )
    if float(transport.get("relative_tolerance", 1.0e-3)) <= 0:
        raise ValueError("Current-control relative_tolerance must be positive.")
    if int(transport.get("max_voltage_iterations", 8)) < 1:
        raise ValueError("max_voltage_iterations must be positive.")
    if float(transport.get("voltage_limit_V", 0.01)) <= 0:
        raise ValueError("voltage_limit_V must be positive.")
    magnetic = config["magnetic_field"]
    if magnetic["gauge"] not in {"symmetric", "landau_x", "landau_y"}:
        raise ValueError("Unknown magnetic vector-potential gauge.")
    if float(magnetic["observation_height_m"]) <= 0:
        raise ValueError("Magnetic observation height must be positive.")
    if int(magnetic["source_stride"]) < 1:
        raise ValueError("Magnetic source_stride must be positive.")
    inset = float(config.get("fluxoid", {}).get("contour_inset_fraction", 0.25))
    if not 0.0 <= inset < 0.5:
        raise ValueError("fluxoid.contour_inset_fraction must lie in [0, 0.5).")
    return config


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args(argv)
    with args.config.open("r", encoding="utf-8") as handle:
        backend_probe = json.load(handle)
    if backend_probe.get("physics_backend", "full_tdgl") == "reduced_vortex":
        from reduced_vortex_diagnostics import load_config as load_reduced_config
        from reduced_vortex_diagnostics import run as run_reduced
        reduced_config = load_reduced_config(args.config)
        summary = run_reduced(reduced_config)
        output_directory = Path(summary["output_directory"])
        print(f"Reduced-vortex diagnostics written to {output_directory}")
        return output_directory
    config = load_config(args.config)
    output_directory = _new_output_directory(config["output"]["directory"])
    cases = build_cases(config.get("sweeps", {}))
    summaries = []
    for index, (label, settings) in enumerate(cases, 1):
        print(f"[{index}/{len(cases)}] {label}", flush=True)
        summaries.append(run_case(label, settings, config, output_directory))
    with (output_directory / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump({"config": config, "cases": summaries}, handle, indent=2)
    _plot_comparison(summaries, output_directory / "comparison.png", int(config["output"]["dpi"]))
    print(f"Diagnostics written to {output_directory}")
    return output_directory


if __name__ == "__main__":
    main()
