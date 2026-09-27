"""Short-time optical vortex response with optional Gaussian pinning.

This is a diagnostic sweep, not a calibrated vortex-force calculation.
"""

import argparse
import json
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shs.config.builder import build_simulation
from shs.config.simulation import PinningSiteConfig
from shs.physics.pinning import pinning_suppression
from shs.solvers.coupled_solver import build_tdgl_model, coupled_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.utils.output import reserve_output_directory
from tdgl_diagnostics import apply_uniform_magnetic_field, _seed_initial_state
from optical_gradient_probe_sweep import _run_branch, _sample, _write_csv


def _position(row):
    return np.array([row["vortex_x_zero_m"], row["vortex_y_zero_m"]], dtype=float)


def _trajectory(laser_rows, dark_rows, direction, sim):
    """Attach local drive proxies and dark-subtracted finite-window velocity."""
    dark_by_step = {row["step"]: row for row in dark_rows}
    pin = pinning_suppression(sim)
    pin_y, pin_x = np.gradient(pin, sim.mesh.dy, sim.mesh.dx)
    tc = float(np.min(sim.material_map.Tc))
    annotated = []
    previous = None
    for row in laser_rows:
        dark = dark_by_step.get(row["step"])
        if dark is None:
            continue
        row = dict(row)
        point = _position(row)
        dark_point = _position(dark)
        valid = np.all(np.isfinite(point)) and np.all(np.isfinite(dark_point))
        delta = float(np.dot(point - dark_point, direction) * 1e9) if valid else float("nan")
        x_cell, y_cell = point / [sim.mesh.dx, sim.mesh.dy] if valid else (float("nan"),) * 2
        pin_gradient = (np.array([_sample(pin_x, x_cell, y_cell),
                                  _sample(pin_y, x_cell, y_cell)])
                        if valid else np.full(2, np.nan))
        thermal_gradient = np.array([
            row["electron_gradient_x_at_core_K_per_m"],
            row["electron_gradient_y_at_core_K_per_m"],
        ])
        # Pinning enters epsilon as -suppression. Its gradient is a TDGL
        # coefficient diagnostic, not a vortex force in N/m.
        row["pinning_suppression_at_core"] = _sample(pin, x_cell, y_cell) if valid else float("nan")
        row["thermal_gradient_toward_beam_K_per_nm"] = float(np.dot(thermal_gradient, direction) * 1e-9)
        row["pinning_coefficient_gradient_toward_beam_per_nm"] = float(np.dot(-pin_gradient, direction) * 1e-9)
        row["dark_subtracted_displacement_toward_beam_nm"] = delta
        row["dark_subtracted_velocity_toward_beam_m_per_s"] = (
            (delta - previous["dark_subtracted_displacement_toward_beam_nm"]) * 1e-9
            / (row["time_s"] - previous["time_s"])
            if previous is not None and np.isfinite(delta)
            and np.isfinite(previous["dark_subtracted_displacement_toward_beam_nm"])
            else float("nan")
        )
        row["subcritical"] = bool(row["maximum_electron_temperature_K"] < tc)
        row["single_vortex"] = bool(row["positive_vortices"] == 1 and row["negative_vortices"] == 0)
        annotated.append(row)
        previous = row
    return annotated


def _plot(rows, cases, output):
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True, constrained_layout=True)
    for case in cases:
        subset = [r for r in rows if r["case_id"] == case["case_id"]]
        if not subset:
            continue
        label = f'{case["pinning_strength"]:g} pin, {case["offset_nm"]:g} nm, {case["power_nW"]:g} nW'
        x = np.array([r["time_s"] for r in subset]) * 1e12
        axes[0].plot(x, [r["dark_subtracted_displacement_toward_beam_nm"] for r in subset], label=label)
        axes[1].plot(x, [r["thermal_gradient_toward_beam_K_per_nm"] for r in subset])
        axes[2].plot(x, [r["dark_subtracted_velocity_toward_beam_m_per_s"] for r in subset])
    axes[0].set_ylabel("Laser-on minus dark motion (nm)")
    axes[1].set_ylabel("Electron gradient toward beam (K/nm)")
    axes[2].set_ylabel("Sampled velocity toward beam (m/s)")
    axes[2].set_xlabel("Time after laser turn-on (ps)")
    for axis in axes:
        axis.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.savefig(output / "response.png", dpi=140)
    plt.close(fig)


def _prepare(config, strength, dark_steps):
    sim = build_simulation(config["simulation_config"])
    sim.config.current = 0.0
    sim.config.electrical.drive_mode = "current"
    sim.config.thermal.model = "two_temperature"
    sim.config.laser.enabled = False
    sim.config.pinning.enabled = strength > 0
    sim.config.pinning.maximum_suppression = float(config["pinning"]["maximum_suppression"])
    seed = np.asarray(config["seed_position_m"], dtype=float)
    site_offset = np.asarray(config["pinning"]["site_offset_m"], dtype=float)
    sim.config.pinning.sites = [PinningSiteConfig(
        *(seed + site_offset), float(config["pinning"]["sigma_m"]), strength)]
    sim.config.pinning.validate()
    apply_uniform_magnetic_field(sim, float(config["applied_field_T"]))
    _seed_initial_state(sim, {"vortices": [config["initial_vortex"]]})
    model = build_tdgl_model(sim)
    scales = tdgl_scales(sim, model)
    dt = sim.config.dt
    for step in range(dark_steps):
        _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=step * dt)
        if not result.converged:
            raise RuntimeError(f"Dark preparation failed at step {step + 1}.")
    return sim, model, scales, deepcopy(sim.fields)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path("tools/config/optical_pinning_response_sweep.json"))
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--dark-steps", type=int, default=None)
    parser.add_argument("--exposure-steps", type=int, default=None)
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    output = reserve_output_directory(config["output_directory"])
    strengths = [float(value) for value in config["pinning"]["strengths"]]
    if not config["pinning"]["enabled"]:
        strengths = [0.0]
    if not strengths or any(not np.isfinite(value) or value < 0 for value in strengths):
        raise ValueError("Pinning strengths must be nonempty, finite and nonnegative.")
    offsets = [float(value) for value in config["sweep"]["offsets_m"]]
    powers = [float(value) for value in config["sweep"]["powers_W"]]
    if not offsets or not powers or any(not np.isfinite(x) or x <= 0 for x in offsets + powers):
        raise ValueError("Offsets and powers must be nonempty, finite and positive.")
    dark_steps = config["timing"]["dark_steps"] if args.dark_steps is None else args.dark_steps
    steps = config["timing"]["exposure_steps"] if args.exposure_steps is None else args.exposure_steps
    sample_every = int(config["timing"]["sample_every"])
    if dark_steps < 0 or steps < 1 or sample_every < 1:
        raise ValueError("Invalid step counts or sampling interval.")
    # Keep unpinned/pinned controls adjacent so --max-cases can make a paired run.
    planned = [(strength, offset, power) for power in powers
               for offset in offsets for strength in strengths]
    if args.max_cases is not None:
        if args.max_cases < 1:
            raise ValueError("--max-cases must be positive.")
        planned = planned[:args.max_cases]
    print(f"Planned {len(planned)} cases; each has laser-on and dark branches.", flush=True)
    trajectories, summaries = [], []
    prepared_by_strength = {}
    dark_by_strength = {}
    seed = np.asarray(config["seed_position_m"], dtype=float)
    direction = np.asarray(config["sweep"]["direction"], dtype=float)
    if direction.shape != (2,) or not np.all(np.isfinite(direction)) or np.linalg.norm(direction) == 0:
        raise ValueError("sweep.direction must be a nonzero [dx,dy] vector.")
    direction /= np.linalg.norm(direction)
    for number, (strength, offset, power) in enumerate(planned, 1):
        if strength not in prepared_by_strength:
            prepared_by_strength[strength] = _prepare(config, strength, dark_steps)
        sim, model, scales, prepared = prepared_by_strength[strength]
        dt = sim.config.dt
        displacement = direction * offset
        if strength not in dark_by_strength:
            dark_by_strength[strength], _ = _run_branch(
                sim, prepared, model, scales, seed, dt, steps, sample_every,
                "dark", None, displacement, power,
                float(config["sweep"]["spot_sigma_m"]))
        dark_rows = dark_by_strength[strength]
        laser_rows, _ = _run_branch(sim, prepared, model, scales, seed, dt, steps,
                                    sample_every, "laser", 1, displacement, power,
                                    float(config["sweep"]["spot_sigma_m"]))
        measured = _trajectory(laser_rows, dark_rows, direction, sim)
        for row in measured:
            row["case_id"] = number
            row["pinning_strength"] = strength
            row["offset_m"] = offset
            row["power_W"] = power
        trajectories.extend(measured)
        valid = bool(measured and all(r["subcritical"] and r["single_vortex"] for r in measured)
                     and measured[-1]["step"] == steps)
        velocities = [r["dark_subtracted_velocity_toward_beam_m_per_s"] for r in measured
                      if np.isfinite(r["dark_subtracted_velocity_toward_beam_m_per_s"])]
        summary = {"case_id": number, "pinning_enabled": strength > 0,
                   "pinning_strength": strength, "offset_nm": offset * 1e9,
                   "power_nW": power * 1e9, "completed_steps": measured[-1]["step"],
                   "valid_subcritical_single_vortex": valid,
                   "final_dark_subtracted_motion_nm": measured[-1]["dark_subtracted_displacement_toward_beam_nm"],
                   "mean_sampled_velocity_m_per_s": float(np.mean(velocities)) if velocities else None,
                   "final_gradient_K_per_nm": measured[-1]["thermal_gradient_toward_beam_K_per_nm"],
                   "peak_electron_temperature_K": max(r["maximum_electron_temperature_K"] for r in measured)}
        summaries.append(summary)
        _write_csv(output / "trajectories.csv", trajectories)
        _write_csv(output / "sweep_summary.csv", summaries)
        (output / "summary.json").write_text(json.dumps({
            "config": config, "effective_dark_steps": dark_steps,
            "effective_exposure_steps": steps, "completed_cases": len(summaries),
            "planned_cases": len(planned), "cases": summaries,
            "interpretation": "Gradients are TDGL drive proxies, not calibrated forces in N/m. "
            "Velocity is a finite difference of dark-subtracted, interpolated vortex positions; "
            "short-window values can be noisy. Pinning strength is dimensionless GL suppression."
        }, indent=2), encoding="utf-8")
        print(f"case {number}: {summary['final_dark_subtracted_motion_nm']:.4g} nm, "
              f"gradient {summary['final_gradient_K_per_nm']:.4g} K/nm, "
              f"valid={valid}", flush=True)
    _plot(trajectories, summaries, output)
    print(output)
    return output


if __name__ == "__main__":
    main()
