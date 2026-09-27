"""Restartable long-duration optical vortex transport with thermal hysteresis."""

import argparse
import csv
import json
import os
import time
from dataclasses import fields as dataclass_fields
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shs.config.builder import build_simulation
from shs.config.simulation import LaserWaypointConfig, PinningSiteConfig
from shs.solvers.coupled_solver import build_tdgl_model, coupled_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.utils.output import reserve_output_directory
from shs.utils.progress import PercentProgress
from tdgl_diagnostics import (apply_uniform_magnetic_field, detect_vortices,
                              _seed_initial_state)
from optical_feedback_tweezer import bounded_command, vortex_center
from optical_gradient_probe import _sample, _snapshot
from shs.utils.constants import (SUPERCONDUCTING_FLUX_QUANTUM,
                                 VACUUM_PERMEABILITY)
from shs.physics.josephson import junction_observables


def thermal_switch(on, peak_K, *, off_K, on_K, dwell, last_change, step):
    """Apply hysteresis only after the configured minimum dwell time."""
    if step - last_change < dwell:
        return on, last_change
    if on and peak_K >= off_K:
        return False, step
    if not on and peak_K <= on_K:
        return True, step
    return on, last_change


def drive_current(time_s, *, enabled, dc_A, ac_A, frequency_Hz, phase_rad=0.0):
    """Return the imposed contact current at the requested physical time."""
    if not enabled:
        return 0.0
    if ac_A != 0 and frequency_Hz <= 0:
        raise ValueError("A nonzero AC current requires a positive frequency_Hz.")
    return float(dc_A + ac_A*np.sin(2*np.pi*frequency_Hz*time_s + phase_rad))


def stalled(progress_history, step, window, minimum_nm, consecutive):
    """Return true only after consecutive completed windows lack progress."""
    if step < window * consecutive or step % window:
        return False
    recent = progress_history[-(consecutive + 1):]
    return len(recent) == consecutive + 1 and all(
        recent[i + 1] - recent[i] < minimum_nm for i in range(consecutive))


def _resume_compatible_config(config):
    """Remove documentation and observational controls from resume checks."""
    return {key: value for key, value in config.items()
            if key not in {"_documentation", "live_preview", "profiling"}}


def _checkpoint(output, sim, state):
    arrays = {}
    scalar_fields = {}
    absent = []
    for item in dataclass_fields(sim.fields):
        value = getattr(sim.fields, item.name)
        if value is None:
            absent.append(item.name)
        elif isinstance(value, np.ndarray):
            arrays[item.name] = value
        else:
            scalar_fields[item.name] = value
    previous_state = (json.loads((output / "state.json").read_text(encoding="utf-8"))
                      if (output / "state.json").exists() else None)
    filename = f'checkpoint_{state["step"]:08d}.npz'
    path = output / filename
    temporary = output / (filename + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, checkpoint_step=np.asarray(state["step"]), **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)
    state = dict(state, absent_fields=absent, scalar_fields=scalar_fields,
                 checkpoint_file=filename)
    temporary = output / "state.json.tmp"
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(state, stream, indent=2)
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(output / "state.json")
    if previous_state is not None:
        old_name = previous_state.get("checkpoint_file")
        if old_name and old_name != filename:
            (output / Path(old_name).name).unlink(missing_ok=True)
    return state


def _restore(output, sim):
    state = json.loads((output / "state.json").read_text(encoding="utf-8"))
    with np.load(output / state.get("checkpoint_file", "checkpoint.npz")) as archive:
        if int(archive["checkpoint_step"]) != int(state["step"]):
            raise RuntimeError("Checkpoint and state metadata have different steps; keep the prior files for recovery.")
        for item in dataclass_fields(sim.fields):
            if item.name in archive:
                value = archive[item.name]
                old = getattr(sim.fields, item.name)
                if isinstance(old, np.ndarray) and value.shape != old.shape:
                    raise ValueError(f"Checkpoint shape differs for {item.name}.")
                setattr(sim.fields, item.name, value.copy())
            elif item.name in state["absent_fields"]:
                setattr(sim.fields, item.name, None)
            else:
                setattr(sim.fields, item.name, state["scalar_fields"][item.name])
    return state


def _truncate_diagnostics(path, step):
    """Discard rows written after the committed checkpoint before resuming."""
    if not path.exists():
        return
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames
        rows = [row for row in reader if int(row["step"]) <= step]
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _append_row(path, row):
    exists = path.exists() and path.stat().st_size > 0
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        if not exists:
            writer.writeheader()
        writer.writerow(row)
        stream.flush()
        os.fsync(stream.fileno())


def _atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2)
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def _timed(timings, name, started):
    entry = timings.setdefault(name, {"seconds": 0.0, "calls": 0})
    entry["seconds"] += time.perf_counter() - started
    entry["calls"] += 1


def _performance_report(timings, started, step, physical_time_s):
    elapsed = time.monotonic() - started
    sections = {}
    accounted = 0.0
    for name, entry in timings.items():
        seconds = float(entry["seconds"])
        accounted += seconds
        calls = int(entry["calls"])
        sections[name] = {
            "seconds": seconds, "calls": calls,
            "mean_seconds_per_call": seconds/calls if calls else 0.0,
            "fraction_of_wall_time": seconds/elapsed if elapsed else 0.0,
        }
    return {"step": int(step), "physical_time_s": float(physical_time_s),
            "wall_seconds": elapsed, "seconds_per_accepted_step":
            elapsed/max(int(step), 1), "accounted_seconds": accounted,
            "unaccounted_seconds": max(elapsed-accounted, 0.0),
            "sections": sections}


def _write_preview_page(directory, refresh_seconds):
    directory.mkdir(exist_ok=True)
    html = f"""<!doctype html><html><head><meta charset=\"utf-8\">
<meta http-equiv=\"refresh\" content=\"{float(refresh_seconds):g}\">
<title>SHS live simulation preview</title>
<style>body{{background:#111;color:#eee;font:16px sans-serif;margin:1rem}}
img{{max-width:100%;height:auto}} code{{color:#9ee}}</style></head><body>
<h1>SHS live simulation preview</h1><p>Automatically refreshes every
{float(refresh_seconds):g} seconds. Create <code>STOP_REQUESTED</code> in the
run directory to stop cleanly.</p><img src=\"latest_state.png\"></body></html>"""
    temporary = directory / "progress.html.tmp"
    temporary.write_text(html, encoding="utf-8")
    temporary.replace(directory / "progress.html")


def _write_live_preview(output, sim, scales, *, step, time_s, center, beam,
                        goal, laser_on, history_frames, status):
    """Atomically publish a compact spatial view for an active run."""
    directory = output / "live"
    directory.mkdir(exist_ok=True)
    f, mesh = sim.fields, sim.mesh
    extent = [mesh.x[0]*1e9, mesh.x[-1]*1e9,
              mesh.y[0]*1e9, mesh.y[-1]*1e9]
    vortices = detect_vortices(
        f.psi, scales.vector_potential_to_dimensionless(f.vector_potential_x),
        scales.vector_potential_to_dimensionless(f.vector_potential_y),
        mesh.dx/scales.xi, mesh.dy/scales.xi)
    positive = np.argwhere(vortices.winding == 1)
    negative = np.argwhere(vortices.winding == -1)
    current_magnitude = np.hypot(f.current_density_x, f.current_density_y)
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    entries = ((f.temperature, "inferno", "Electron temperature (K)", None),
               (np.abs(f.psi), "viridis", "|psi|", None),
               (np.angle(f.psi), "twilight", "Phase", (-np.pi, np.pi)),
               (current_magnitude, "magma", "Current density magnitude", None))
    for axis, (field, cmap, title, limits) in zip(axes.flat, entries):
        image = axis.imshow(field, origin="lower", extent=extent, cmap=cmap,
                            vmin=None if limits is None else limits[0],
                            vmax=None if limits is None else limits[1])
        if len(positive):
            axis.scatter((positive[:, 1]+.5)*mesh.dx*1e9,
                         (positive[:, 0]+.5)*mesh.dy*1e9,
                         marker="o", facecolors="none", edgecolors="cyan",
                         linewidths=1.2, label="+ vortex")
        if len(negative):
            axis.scatter((negative[:, 1]+.5)*mesh.dx*1e9,
                         (negative[:, 0]+.5)*mesh.dy*1e9,
                         marker="x", color="lime", label="- vortex")
        axis.plot(center[0]*1e9, center[1]*1e9, "wx", markersize=8,
                  markeredgewidth=2, label="tracked")
        axis.plot(beam[0]*1e9, beam[1]*1e9, "c+", markersize=9,
                  markeredgewidth=2, label="laser")
        axis.plot(goal[0]*1e9, goal[1]*1e9, "w*", markersize=8, label="target")
        junction = getattr(sim.config, "josephson", None)
        if junction is not None and junction.enabled:
            axis.axvspan((junction.center_x_m-.5*junction.width_m)*1e9,
                         (junction.center_x_m+.5*junction.width_m)*1e9,
                         color="cyan", alpha=.14)
        axis.set(xlabel="x (nm)", ylabel="y (nm)", title=title)
        fig.colorbar(image, ax=axis, shrink=.83)
    axes[0, 0].legend(fontsize=7, loc="best")
    fig.suptitle(f"step={step}  t={time_s*1e12:.3f} ps  "
                 f"vortices=+{vortices.positive_count}/-{vortices.negative_count}  "
                 f"laser={'on' if laser_on else 'off'}")
    frame = directory / f"preview_{step:08d}.png"
    temporary = frame.with_name(frame.name + ".tmp")
    fig.savefig(temporary, format="png", dpi=110)
    plt.close(fig)
    temporary.replace(frame)
    latest = directory / "latest_state.png"
    latest_tmp = directory / "latest_state.png.tmp"
    latest_tmp.write_bytes(frame.read_bytes())
    latest_tmp.replace(latest)
    retained = sorted(directory.glob("preview_*.png"))
    for old in retained[:-history_frames]:
        old.unlink(missing_ok=True)
    _atomic_json(directory / "latest_status.json", status)


def _plot(output):
    path = output / "diagnostics.csv"
    if not path.exists():
        return
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < 2:
        return
    t = np.array([float(r["time_s"]) * 1e12 for r in rows])
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    axes[0].plot(t, [float(r["progress_nm"]) for r in rows])
    axes[0].set_ylabel("Targetward motion (nm)")
    axes[1].plot(t, [float(r["maximum_electron_temperature_K"]) for r in rows])
    axes[1].plot(t, [float(r["electron_temperature_at_core_K"]) for r in rows])
    axes[1].set_ylabel("Electron T (K)")
    axes[2].plot(t, [float(r["gradient_toward_target_K_per_nm"]) for r in rows])
    axes[2].set_ylabel("Gradient at vortex (K/nm)")
    axes[2].set_xlabel("Physical time (ps)")
    for axis in axes:
        axis.grid(alpha=0.3)
    fig.savefig(output / "transport_diagnostics.png", dpi=140)
    plt.close(fig)
    if "thermal_force_london_toward_target_N_per_m" in rows[0]:
        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True,
                                 constrained_layout=True)
        axes[0].plot(t, [float(r["sampled_velocity_m_per_s"]) for r in rows])
        axes[0].set_ylabel("Sampled targetward velocity (m/s)")
        if "sampled_acceleration_m_per_s2" in rows[0]:
            acceleration_axis = axes[0].twinx()
            acceleration_axis.plot(
                t, [float(r["sampled_acceleration_m_per_s2"]) for r in rows],
                color="tab:purple", alpha=.55)
            acceleration_axis.set_ylabel("Acceleration (m/s²)", color="tab:purple")
        axes[1].plot(t, [float(r["thermal_force_london_toward_target_N_per_m"])
                         for r in rows], label="London thermal estimate")
        axes[1].plot(t, [float(r["lorentz_force_toward_target_N_per_m"])
                         for r in rows], label="local J x Phi0")
        axes[1].set_ylabel("Force/length estimate (N/m)"); axes[1].legend(fontsize=8)
        axes[2].plot(t, [float(r["target_distance_nm"]) for r in rows], label="target distance")
        axes[2].plot(t, [float(r["nearest_boundary_distance_nm"]) for r in rows],
                     label="nearest boundary")
        axes[2].set_ylabel("Distance (nm)"); axes[2].set_xlabel("Physical time (ps)")
        axes[2].legend(fontsize=8)
        for axis in axes:
            axis.grid(alpha=.3)
        fig.savefig(output / "motion_and_force_proxies.png", dpi=140)
        plt.close(fig)
    if "junction_voltage_V" in rows[0]:
        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True,
                                 constrained_layout=True)
        axes[0].plot(t, [float(r["applied_current_A"])*1e6 for r in rows])
        axes[0].set_ylabel("Applied current (µA)")
        axes[1].plot(t, [float(r["junction_voltage_V"])*1e6 for r in rows],
                     label="junction voltage")
        axes[1].set_ylabel("Voltage (µV)")
        axes[2].plot(t, [float(r["junction_phase_difference_rad"]) for r in rows],
                     label="phase difference")
        axes[2].plot(t, [float(r["junction_mean_amplitude"]) for r in rows],
                     label="mean weak-link |psi|")
        axes[2].set_ylabel("Phase (rad) / amplitude")
        axes[2].set_xlabel("Physical time (ps)")
        axes[2].legend(fontsize=8)
        for axis in axes:
            axis.grid(alpha=.3)
        fig.savefig(output / "junction_dynamics.png", dpi=140)
        plt.close(fig)


def _temperature_at(sim, position_m):
    cell = np.asarray(position_m) / [sim.mesh.dx, sim.mesh.dy]
    return _sample(sim.fields.temperature, cell[0], cell[1])


def _save_field_frame(output, sim, step, time_s, center, beam, laser_on):
    directory = output / "field_frames"
    directory.mkdir(exist_ok=True)
    path = directory / f"frame_{step:08d}.npz"
    temporary = directory / (path.name + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(
            stream, step=np.asarray(step), time_s=np.asarray(time_s),
            temperature=sim.fields.temperature,
            phonon_temperature=sim.fields.phonon_temperature,
            psi=sim.fields.psi, vortex_center_m=np.asarray(center),
            vector_potential_x=sim.fields.vector_potential_x,
            vector_potential_y=sim.fields.vector_potential_y,
            beam_m=np.asarray(beam), laser_on=np.asarray(laser_on))
        stream.flush(); os.fsync(stream.fileno())
    temporary.replace(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path("tools/config/optical_transport_overnight.json"))
    parser.add_argument("--resume", type=Path, help="Resume an existing checkpoint directory.")
    parser.add_argument("--max-steps", type=int, help="Smoke-test maximum physical step.")
    parser.add_argument("--dark-steps", type=int, help="Smoke-test dark preparation override.")
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    output = args.resume if args.resume else reserve_output_directory(config["output_directory"])
    sim = build_simulation(config["simulation_config"])
    drive = config.get("electrical_drive", {})
    drive_enabled = bool(drive.get("enabled", False))
    dc_current = float(drive.get("dc_current_A", 0.0))
    ac_current = float(drive.get("ac_amplitude_A", 0.0))
    drive_frequency = float(drive.get("frequency_Hz", 0.0))
    drive_phase = float(drive.get("phase_rad", 0.0))
    if ac_current != 0 and drive_frequency <= 0:
        raise ValueError("A nonzero AC current requires a positive frequency_Hz.")
    sim.config.current = dc_current if drive_enabled else 0.0
    sim.config.electrical.drive_mode = "current"
    sim.config.thermal.model = "two_temperature"
    sim.config.laser.enabled = False
    sim.config.laser.absorbed_power_W = float(config["laser"]["absorbed_power_W"])
    sim.config.laser.sigma_m = float(config["laser"]["spot_sigma_m"])
    pin = config["pinning"]
    sim.config.pinning.enabled = bool(pin["enabled"])
    sim.config.pinning.sites = [PinningSiteConfig(
        *site["position_m"], site["sigma_m"], site["strength"])
        for site in pin["sites"]] if pin["enabled"] else []
    sim.config.pinning.validate()
    apply_uniform_magnetic_field(sim, float(config["applied_field_T"]))
    model = build_tdgl_model(sim)
    scales = tdgl_scales(sim, model)
    dt = sim.config.dt
    target_displacement = np.asarray(config["target"]["displacement_m"], dtype=float)
    offset = np.asarray(config["laser"]["offset_m"], dtype=float)
    if target_displacement.shape != (2,) or np.linalg.norm(target_displacement) <= 0:
        raise ValueError("target.displacement_m must be a nonzero 2-vector.")
    if offset.shape != (2,) or np.linalg.norm(offset) <= 0:
        raise ValueError("laser.offset_m must be a nonzero 2-vector.")
    direction = target_displacement / np.linalg.norm(target_displacement)
    thermal = config["thermal_control"]
    tc = float(np.min(sim.material_map.Tc))
    hard_stop = thermal.get("hard_stop_K")
    if not (sim.config.temperature < thermal["restart_below_K"]
            < thermal["switch_off_K"] < tc):
        raise ValueError("Require base T < restart < switch-off < Tc.")
    if hard_stop is not None and hard_stop <= thermal["switch_off_K"]:
        raise ValueError("An optional hard stop must exceed the switch-off temperature.")
    sensor = thermal.get("sensor", "peak")
    if sensor not in {"peak", "vortex"}:
        raise ValueError("thermal_control.sensor must be peak or vortex.")
    timing = config["timing"]
    tracking = config.get("tracking", {})
    require_unique = bool(tracking.get("require_unique_topology", True))
    maximum_jump_m = tracking.get("maximum_jump_m")
    initial_search_radius_m = tracking.get("initial_search_radius_m")
    max_steps = int(timing["max_steps"] if args.max_steps is None else args.max_steps)
    if max_steps < 1 or timing["sample_every"] < 1 or timing["checkpoint_every"] < 1:
        raise ValueError("Step counts and output intervals must be positive.")
    preview = config.get("live_preview", {})
    preview_enabled = bool(preview.get("enabled", False))
    preview_every = int(preview.get("every_steps", 200))
    preview_history = int(preview.get("history_frames", 20))
    refresh_seconds = float(preview.get("page_refresh_seconds", 10))
    if preview_every < 1 or preview_history < 1 or refresh_seconds <= 0:
        raise ValueError("Live-preview cadence, history and refresh must be positive.")
    profiling_enabled = bool(config.get("profiling", {}).get("enabled", True))
    timings = {}
    start_wall = time.monotonic()
    if args.resume:
        state = _restore(output, sim)
        saved_config = json.loads((output / "config.json").read_text(encoding="utf-8"))
        if (_resume_compatible_config(saved_config) != _resume_compatible_config(config)
                or state["dt_s"] != dt):
            raise ValueError("Resume config or timestep differs from checkpoint.")
        step0 = int(state["step"])
        _truncate_diagnostics(output / "diagnostics.csv", step0)
        initial = np.asarray(state["initial_center_m"])
        center = np.asarray(state["last_center_m"])
        beam = np.asarray(state["beam_m"])
        target_beam = np.asarray(state.get("target_beam_m", center + offset))
        on = bool(state["laser_on"])
        last_change = int(state["last_switch_step"])
        progress_history = list(state["progress_window_nm"])
        last_sample_step = int(state.get("last_sample_step", step0))
        last_sample_progress = float(state.get("last_sample_progress_nm", state["progress_nm"]))
        last_sample_velocity = float(state.get("last_sample_velocity_m_per_s", float("nan")))
    else:
        (output / "config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
        initial_vortices = config.get("initial_vortices")
        if initial_vortices is None:
            initial_vortices = [config["initial_vortex"]]
        if not initial_vortices:
            raise ValueError("At least one initial vortex is required.")
        _seed_initial_state(sim, {"vortices": initial_vortices})
        dark_steps = int(timing["dark_steps"] if args.dark_steps is None else args.dark_steps)
        if dark_steps < 0:
            raise ValueError("Dark preparation steps cannot be negative.")
        for step in range(dark_steps):
            operation_started = time.perf_counter()
            _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=step * dt)
            _timed(timings, "dark_solver", operation_started)
            if not result.converged:
                raise RuntimeError(f"Dark preparation failed at step {step+1}.")
        initial = vortex_center(
            sim, scales, np.asarray(config["seed_position_m"]),
            require_unique=require_unique,
            maximum_jump_m=initial_search_radius_m)
        if initial is None:
            seed = np.asarray(config["seed_position_m"], dtype=float)
            vortices = detect_vortices(
                sim.fields.psi,
                scales.vector_potential_to_dimensionless(sim.fields.vector_potential_x),
                scales.vector_potential_to_dimensionless(sim.fields.vector_potential_y),
                sim.mesh.dx/scales.xi, sim.mesh.dy/scales.xi)
            _save_field_frame(output, sim, dark_steps, dark_steps*dt,
                              seed, seed+offset, False)
            failure = {
                "reason": "prepared_vortex_not_found_within_initial_search_radius",
                "dark_steps": dark_steps,
                "positive_vortices": vortices.positive_count,
                "negative_vortices": vortices.negative_count,
                "seed_position_m": seed.tolist(),
                "initial_search_radius_m": initial_search_radius_m,
                "applied_field_T": float(config["applied_field_T"]),
            }
            _atomic_json(output / "preparation_failure.json", failure)
            raise RuntimeError(
                "Prepared tracked vortex was not found within the initial "
                f"search radius; detected +{vortices.positive_count}/"
                f"-{vortices.negative_count}. See preparation_failure.json.")
        center = initial.copy()
        beam = center + offset
        target_beam = beam.copy()
        on, last_change, step0 = True, -int(thermal["minimum_dwell_steps"]), 0
        progress_history = [0.0]
        last_sample_step, last_sample_progress = 0, 0.0
        last_sample_velocity = float("nan")
    goal = initial + target_displacement
    snapshot_every = int(config.get("output", {}).get("field_snapshot_every", 0))
    if snapshot_every < 0:
        raise ValueError("output.field_snapshot_every cannot be negative.")
    if snapshot_every and not args.resume:
        _save_field_frame(output, sim, 0, 0.0, center, beam, on)
    if max_steps <= step0:
        raise ValueError("Maximum step must exceed the checkpoint step.")
    progress_reporter = PercentProgress(
        max_steps, label="optical transport", completed_steps=step0)
    state = {}
    reason = "maximum_steps"
    run_start_step = step0
    if preview_enabled:
        _write_preview_page(output / "live", refresh_seconds)
    for step in range(step0 + 1, max_steps + 1):
        peak_before = float(np.max(sim.fields.temperature))
        control_temperature = (_temperature_at(sim, center)
                               if sensor == "vortex" else peak_before)
        on, last_change = thermal_switch(
            on, control_temperature, off_K=float(thermal["switch_off_K"]),
            on_K=float(thermal["restart_below_K"]),
            dwell=int(thermal["minimum_dwell_steps"]), last_change=last_change, step=step)
        sim.config.laser.enabled = on
        start = (step - 1) * dt
        midpoint = start + .5*dt
        applied_current = drive_current(
            midpoint, enabled=drive_enabled, dc_A=dc_current, ac_A=ac_current,
            frequency_Hz=drive_frequency, phase_rad=drive_phase)
        sim.config.current = float(applied_current)
        if (step - 1) % int(config["laser"]["update_every_steps"]) == 0:
            target_beam = center + offset
        next_beam = bounded_command(beam, target_beam,
                                    float(config["laser"]["max_speed_m_s"]) * dt)
        sim.config.laser.waypoints = [LaserWaypointConfig(start, *beam),
                                      LaserWaypointConfig(start + dt, *next_beam)]
        if on:
            sim.config.laser.validate()
        operation_started = time.perf_counter()
        _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=start)
        _timed(timings, "coupled_solver", operation_started)
        progress_reporter.update(step)
        if not result.converged:
            reason = "solver_nonconvergence"
            break
        beam = next_beam
        operation_started = time.perf_counter()
        measured = vortex_center(
            sim, scales, center, require_unique=require_unique,
            maximum_jump_m=maximum_jump_m)
        _timed(timings, "vortex_tracking", operation_started)
        if measured is None:
            reason = "tracked_vortex_unresolved"
            if snapshot_every:
                _save_field_frame(output, sim, step, step*dt, center, beam, on)
            break
        center = measured
        peak = float(np.max(sim.fields.temperature))
        if hard_stop is not None and peak >= float(hard_stop):
            reason = "hard_temperature_limit"
        progress = float(np.dot(center - initial, direction) * 1e9)
        remaining = float(np.linalg.norm(goal - center) * 1e9)
        if remaining <= float(config["target"]["tolerance_nm"]):
            reason = "target_reached"
        write_row = (step % int(timing["sample_every"]) == 0 or
                     step % int(timing["checkpoint_every"]) == 0 or
                     reason != "maximum_steps" or step == max_steps)
        if write_row:
            operation_started = time.perf_counter()
            sample_dt = (step-last_sample_step)*dt
            sample_velocity = ((progress-last_sample_progress)*1e-9/sample_dt
                               if step % int(timing["sample_every"]) == 0
                               and sample_dt > 0 else float("nan"))
            sample_acceleration = ((sample_velocity-last_sample_velocity)/sample_dt
                                   if np.isfinite(sample_velocity)
                                   and np.isfinite(last_sample_velocity) else float("nan"))
            row = _snapshot(sim, model, scales, center, step * dt, step,
                            "transport", "feedback")
            # The continuation tracker has robust fallbacks for an
            # ill-conditioned complex-zero fit.  Use that accepted position
            # as the trajectory source of truth while retaining the raw fit's
            # validity flag from _snapshot.
            row["vortex_x_zero_m"] = float(center[0])
            row["vortex_y_zero_m"] = float(center[1])
            row["tracking_core_valid"] = True
            row.update(laser_on=on, commanded_beam_x_m=float(beam[0]),
                       commanded_beam_y_m=float(beam[1]),
                       applied_current_A=float(applied_current),
                       microwave_phase_rad=float((2*np.pi*drive_frequency*step*dt
                                                  + drive_phase) % (2*np.pi)),
                       progress_nm=progress, target_distance_nm=remaining,
                       gradient_toward_target_K_per_nm=float(np.dot(direction, [
                           row["electron_gradient_x_at_core_K_per_m"],
                           row["electron_gradient_y_at_core_K_per_m"]]) * 1e-9),
                       sampled_velocity_m_per_s=sample_velocity,
                       sampled_acceleration_m_per_s2=sample_acceleration)
            cell = center / [sim.mesh.dx, sim.mesh.dy]
            current = np.array([
                _sample(sim.fields.current_density_x, cell[0], cell[1]),
                _sample(sim.fields.current_density_y, cell[0], cell[1])])
            lorentz = SUPERCONDUCTING_FLUX_QUANTUM * np.array([current[1], -current[0]])
            penetration = _sample(sim.material_map.penetration_depth, cell[0], cell[1])
            coherence = _sample(sim.material_map.coherence_length, cell[0], cell[1])
            local_tc = _sample(sim.material_map.Tc, cell[0], cell[1])
            london_coefficient = (SUPERCONDUCTING_FLUX_QUANTUM**2
                * np.log(max(penetration/coherence, 1.0))
                / (4*np.pi*VACUUM_PERMEABILITY*penetration**2*local_tc))
            row.update(
                topology_is_unique=bool(row["positive_vortices"] == 1
                                        and row["negative_vortices"] == 0),
                tracking_allows_additional_vortices=not require_unique,
                thermal_control_temperature_K=control_temperature,
                thermal_force_london_toward_target_N_per_m=london_coefficient
                * row["gradient_toward_target_K_per_nm"] * 1e9,
                lorentz_force_toward_target_N_per_m=float(np.dot(lorentz, direction)),
                gl_gradient_toward_target_per_nm=float(np.dot(direction, [
                    row["gl_coefficient_gradient_x_at_core_per_m"],
                    row["gl_coefficient_gradient_y_at_core_per_m"]]) * 1e-9),
                nearest_boundary_distance_nm=float(min(
                    center[0], sim.mesh.x[-1]-center[0], center[1],
                    sim.mesh.y[-1]-center[1]) * 1e9))
            junction = getattr(sim.config, "josephson", None)
            if junction is not None and junction.enabled:
                row.update(junction_observables(sim))
                voltage_scale = SUPERCONDUCTING_FLUX_QUANTUM*drive_frequency
                row["shapiro_voltage_normalized"] = (
                    row["junction_voltage_V"]/voltage_scale
                    if drive_enabled and drive_frequency > 0 else float("nan"))
            _append_row(output / "diagnostics.csv", row)
            _timed(timings, "diagnostics", operation_started)
            if step % int(timing["sample_every"]) == 0:
                last_sample_progress, last_sample_step = progress, step
                last_sample_velocity = sample_velocity
        if snapshot_every and (step % snapshot_every == 0 or reason != "maximum_steps"
                               or step == max_steps):
            operation_started = time.perf_counter()
            _save_field_frame(output, sim, step, step*dt, center, beam, on)
            _timed(timings, "field_snapshot", operation_started)
        if step % int(timing["stall_window_steps"]) == 0:
            progress_history.append(progress)
            if stalled(progress_history, step, int(timing["stall_window_steps"]),
                       float(timing["minimum_progress_per_window_nm"]),
                       int(timing["stall_consecutive_windows"])):
                reason = "stalled"
        if time.monotonic() - start_wall >= float(timing["max_wall_hours"]) * 3600:
            reason = "wall_time_limit"
        stop_request = output / "STOP_REQUESTED"
        if stop_request.exists():
            reason = "user_stop_requested"
            stop_request.unlink(missing_ok=True)
        state = {"step": step, "dt_s": dt, "initial_center_m": initial.tolist(),
                 "last_center_m": center.tolist(), "beam_m": beam.tolist(),
                 "target_beam_m": target_beam.tolist(),
                 "laser_on": on, "last_switch_step": last_change,
                 "progress_window_nm": progress_history, "reason": reason,
                 "last_sample_step": last_sample_step,
                 "last_sample_progress_nm": last_sample_progress,
                 "last_sample_velocity_m_per_s": last_sample_velocity,
                 "peak_electron_temperature_K": peak, "progress_nm": progress,
                 "target_distance_nm": remaining}
        if preview_enabled and (step % preview_every == 0
                                or reason != "maximum_steps" or step == max_steps):
            operation_started = time.perf_counter()
            preview_status = {
                "step": step, "physical_time_s": step*dt,
                "reason": reason, "progress_nm": progress,
                "target_distance_nm": remaining,
                "peak_electron_temperature_K": peak,
                "laser_on": on,
                "stop_file": str(output / "STOP_REQUESTED"),
            }
            if profiling_enabled:
                preview_status["performance"] = _performance_report(
                    timings, start_wall, step-run_start_step, step*dt)
            _write_live_preview(
                output, sim, scales, step=step, time_s=step*dt,
                center=center, beam=beam, goal=goal, laser_on=on,
                history_frames=preview_history, status=preview_status)
            _timed(timings, "live_preview", operation_started)
        if step % int(timing["checkpoint_every"]) == 0 or reason != "maximum_steps" or step == max_steps:
            operation_started = time.perf_counter()
            state = _checkpoint(output, sim, state)
            _timed(timings, "checkpoint", operation_started)
            if profiling_enabled:
                _atomic_json(output / "performance.json", _performance_report(
                    timings, start_wall, step-run_start_step, step*dt))
            print(f"step={step} progress={progress:.4g} nm remaining={remaining:.4g} nm "
                  f"Te_max={peak:.4g} K laser={'on' if on else 'off'}", flush=True)
        if reason != "maximum_steps":
            break
    if not state:
        state = {"step": step0, "reason": reason}
    status = {"reason": reason, "last_committed_step": state["step"],
              "checkpoint_file": state.get("checkpoint_file"),
              "progress_nm": state.get("progress_nm"),
              "target_distance_nm": state.get("target_distance_nm"),
              "elapsed_wall_seconds": time.monotonic() - start_wall}
    if profiling_enabled:
        status["performance_file"] = "performance.json"
        _atomic_json(output / "performance.json", _performance_report(
            timings, start_wall, max(int(state["step"])-run_start_step, 0),
            int(state["step"])*dt))
    _atomic_json(output / "run_status.json", status)
    _plot(output)
    print(f"{reason}: {output}", flush=True)
    return output


if __name__ == "__main__":
    main()
