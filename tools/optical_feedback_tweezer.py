"""Feedback-guided stationary-offset optical vortex transport experiment.

The beam command uses only the previously accepted TDGL state. It is an
external feedback protocol, not an additional vortex force in the TDGL PDE.
"""
import argparse
import csv
import json
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from shs.config.builder import build_simulation
from shs.config.simulation import LaserWaypointConfig
from shs.solvers.coupled_solver import build_tdgl_model, coupled_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.utils.output import reserve_output_directory
from tdgl_diagnostics import apply_uniform_magnetic_field, detect_vortices, _seed_initial_state
from optical_gradient_probe import _snapshot, complex_zero_core, subcell_core


def robust_core_position(psi, plaquette_x, plaquette_y):
    """Locate a vortex core without discarding a valid winding on a bad fit."""
    x_cell, y_cell = plaquette_x + .5, plaquette_y + .5
    cx, cy = complex_zero_core(psi, x_cell, y_cell)
    if np.all(np.isfinite([cx, cy])):
        return np.array([cx, cy]), "complex_zero"
    cx, cy = subcell_core(np.abs(psi), x_cell, y_cell, radius=3)
    if (np.all(np.isfinite([cx, cy]))
            and abs(cx-x_cell) <= 4 and abs(cy-y_cell) <= 4):
        return np.array([cx, cy]), "amplitude_minimum"
    return np.array([x_cell, y_cell]), "plaquette_center"


def vortex_center(sim, scales, reference_m, *, require_unique=True,
                  maximum_jump_m=None):
    """Return the nearest tracked +1 complex zero.

    Strict uniqueness remains the default. Identity-tracking experiments may
    allow other vortices while limiting how far the selected core can jump.
    """
    f, mesh = sim.fields, sim.mesh
    vortices = detect_vortices(
        f.psi, scales.vector_potential_to_dimensionless(f.vector_potential_x),
        scales.vector_potential_to_dimensionless(f.vector_potential_y),
        mesh.dx/scales.xi, mesh.dy/scales.xi)
    if require_unique and (vortices.positive_count != 1 or vortices.negative_count):
        return None
    candidates = []
    for y, x in np.argwhere(vortices.winding == 1):
        core, _ = robust_core_position(f.psi, x, y)
        candidates.append(core*np.array([mesh.dx, mesh.dy]))
    if not candidates:
        return None
    center = min(candidates, key=lambda point: np.linalg.norm(point-reference_m))
    limit = (0.25*min(mesh.x[-1], mesh.y[-1]) if maximum_jump_m is None
             else float(maximum_jump_m))
    if not np.isfinite(limit) or limit <= 0:
        raise ValueError("maximum_jump_m must be positive and finite.")
    if np.linalg.norm(center-reference_m) > limit:
        return None
    return center


def bounded_command(current, target, maximum_step_m):
    delta = np.asarray(target)-np.asarray(current)
    distance = float(np.linalg.norm(delta))
    if distance <= maximum_step_m:
        return np.asarray(target, dtype=float)
    return np.asarray(current, dtype=float) + delta*(maximum_step_m/distance)


def _write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def _frame(sim, row):
    fields = sim.fields
    return {"time_s": row["time_s"], "temperature": fields.temperature.copy(),
            "phonon_temperature": fields.phonon_temperature.copy(),
            "amplitude": np.abs(fields.psi).copy(),
            "phase": np.angle(fields.psi).copy(),
            "beam_x_m": row["laser_x_m"], "beam_y_m": row["laser_y_m"],
            "vortex_x_m": row["vortex_x_zero_m"],
            "vortex_y_m": row["vortex_y_zero_m"]}


def _plot_motion(branches, output, offset_m, dpi):
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True,
                             constrained_layout=True)
    for mode, rows in branches.items():
        time = np.array([row["time_s"] for row in rows])*1e12
        axes[0].plot(time, [row["vortex_y_zero_m"]*1e9 for row in rows], label=mode)
        axes[1].plot(time, [row["electron_temperature_at_core_K"] for row in rows],
                     label=mode)
        axes[2].plot(time, [row["electron_gradient_y_at_core_K_per_m"]*1e-9
                              for row in rows], label=mode)
        if mode != "dark":
            axes[0].plot(time, [row["laser_y_m"]*1e9 for row in rows],
                         linestyle="--", alpha=.6, label=f"{mode} beam")
    axes[0].set_ylabel("y position (nm)")
    axes[1].set_ylabel("Electron T at core (K)")
    axes[2].set_ylabel("dTe/dy at core (K/nm)")
    axes[2].set_xlabel("Time after dark preparation (ps)")
    axes[0].set_title(f"Beam commanded {np.linalg.norm(offset_m)*1e9:g} nm ahead of vortex")
    for axis in axes:
        axis.grid(alpha=.25)
        axis.legend(fontsize=8)
    fig.savefig(output/"motion_and_gradient.png", dpi=dpi)
    plt.close(fig)
    fig, axis = plt.subplots(figsize=(9, 4), constrained_layout=True)
    for mode, rows in branches.items():
        if mode == "dark":
            continue
        axis.plot([row["time_s"]*1e12 for row in rows],
                  [row["beam_vortex_separation_m"]*1e9 for row in rows], label=mode)
    axis.axhline(np.linalg.norm(offset_m)*1e9, color="black", linestyle=":",
                 label="requested offset")
    axis.set(xlabel="Time (ps)", ylabel="Beam-vortex separation (nm)")
    axis.grid(alpha=.25); axis.legend()
    fig.savefig(output/"beam_offset_error.png", dpi=dpi)
    plt.close(fig)


def _plot_relative_motion(branches, output, dpi):
    dark = {row["step"]: row for row in branches["dark"]}
    fig, axis = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for mode in ("dynamic", "stationary"):
        aligned = [(row, dark[row["step"]]) for row in branches[mode]
                   if row["step"] in dark]
        axis.plot([row["time_s"]*1e12 for row, _ in aligned],
                  [(row["vortex_y_zero_m"]-background["vortex_y_zero_m"])*1e9
                   for row, background in aligned], label=f"{mode} minus dark")
    axis.axhline(0, color="black", linewidth=.8)
    axis.set(xlabel="Time after preparation (ps)",
             ylabel="Vortex y displacement relative to dark (nm)",
             title="Subcell vortex response")
    axis.grid(alpha=.25); axis.legend()
    fig.savefig(output/"vortex_motion_zoom.png", dpi=dpi)
    plt.close(fig)


def _plot_film(frames, output, mesh, dpi, zoom=None):
    selected = sorted(set([0, len(frames)//2, len(frames)-1]))
    temperature_limits = (min(float(np.min(frame["temperature"])) for frame in frames),
                          max(float(np.max(frame["temperature"])) for frame in frames))
    amplitude_max = max(float(np.max(frame["amplitude"])) for frame in frames)
    fig, axes = plt.subplots(3, len(selected), figsize=(4*len(selected), 10),
                             constrained_layout=True)
    for column, index in enumerate(selected):
        frame = frames[index]
        for row, (key, cmap, label) in enumerate(
            (("temperature", "inferno", "Electron T (K)"),
             ("amplitude", "viridis", "|psi|"),
             ("phase", "twilight", "Phase"))):
            axis = axes[row, column]
            limits = (temperature_limits if key == "temperature" else
                      (0, amplitude_max) if key == "amplitude" else (-np.pi, np.pi))
            image = axis.imshow(frame[key], origin="lower", cmap=cmap,
                                vmin=limits[0], vmax=limits[1])
            if frame["beam_x_m"] is not None and np.isfinite(frame["beam_x_m"]):
                axis.plot(frame["beam_x_m"]/mesh.dx,
                          frame["beam_y_m"]/mesh.dy, "c+", markersize=8)
            if np.isfinite(frame["vortex_x_m"]):
                axis.plot(frame["vortex_x_m"]/mesh.dx,
                          frame["vortex_y_m"]/mesh.dy, "rx", markersize=6)
            axis.set_title(f"{label}, t={frame['time_s']*1e12:.2f} ps")
            if zoom is not None:
                center, radius = zoom
                axis.set_xlim(center[0]-radius, center[0]+radius)
                axis.set_ylim(center[1]-radius, center[1]+radius)
            fig.colorbar(image, ax=axis)
    fig.savefig(output/("film_snapshots_zoom.png" if zoom else "film_snapshots.png"),
                dpi=dpi)
    plt.close(fig)


def _save_gif(frames, output, mesh, fps, dpi, zoom=None):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    keys = ("temperature", "amplitude", "phase")
    cmaps = ("inferno", "viridis", "twilight")
    limits = [(min(float(np.min(frame[key])) for frame in frames),
               max(float(np.max(frame[key])) for frame in frames)) for key in keys]
    images, beams, cores = [], [], []
    for axis, key, cmap, (minimum, maximum) in zip(axes, keys, cmaps, limits):
        image = axis.imshow(frames[0][key], origin="lower", cmap=cmap,
                            vmin=minimum, vmax=maximum)
        beam, = axis.plot([], [], "c+", markersize=8)
        core, = axis.plot([], [], "rx", markersize=6)
        axis.set_title(key)
        if zoom is not None:
            center, radius = zoom
            axis.set_xlim(center[0]-radius, center[0]+radius)
            axis.set_ylim(center[1]-radius, center[1]+radius)
        fig.colorbar(image, ax=axis)
        images.append(image); beams.append(beam); cores.append(core)
    title = fig.suptitle("")
    def update(index):
        frame = frames[index]
        for image, key, beam, core in zip(images, keys, beams, cores):
            image.set_data(frame[key])
            if frame["beam_x_m"] is not None and np.isfinite(frame["beam_x_m"]):
                beam.set_data([frame["beam_x_m"]/mesh.dx],
                              [frame["beam_y_m"]/mesh.dy])
            else:
                beam.set_data([], [])
            if np.isfinite(frame["vortex_x_m"]):
                core.set_data([frame["vortex_x_m"]/mesh.dx],
                              [frame["vortex_y_m"]/mesh.dy])
            else:
                core.set_data([], [])
        title.set_text(f"t = {frame['time_s']*1e12:.3f} ps; red = vortex, cyan = beam")
        return images+beams+cores+[title]
    movie = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
    movie.save(output/("feedback_motion_zoom.gif" if zoom else "feedback_motion.gif"),
               writer=animation.PillowWriter(fps=fps),
               dpi=dpi)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path("tools/config/optical_feedback_tweezer.json"))
    parser.add_argument("--max-steps", type=int, help="Override exposure steps for a smoke test.")
    parser.add_argument("--dark-steps", type=int, help="Override dark preparation steps for a smoke test.")
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    output = reserve_output_directory(config["output_directory"])
    sim = build_simulation(config["simulation_config"])
    sim.config.current = 0.0
    sim.config.electrical.drive_mode = "current"
    sim.config.thermal.model = "two_temperature"
    sim.config.thermal.validate()
    sim.config.laser.enabled = False
    apply_uniform_magnetic_field(sim, config["applied_field_T"])
    _seed_initial_state(sim, {"vortices": [config["initial_vortex"]]})
    model = build_tdgl_model(sim)
    scales = tdgl_scales(sim, model)
    dt = sim.config.dt
    reference = np.asarray(config["seed_position_m"], dtype=float)
    offset = np.asarray(config["controller"]["offset_m"], dtype=float)
    if offset.shape != (2,) or np.linalg.norm(offset) <= 0:
        raise ValueError("controller.offset_m must be a nonzero [dx,dy] vector.")
    if config["controller"]["max_speed_m_s"] <= 0:
        raise ValueError("controller.max_speed_m_s must be positive.")
    dark_changes = []
    dark_steps = (int(config["timing"]["dark_steps"]) if args.dark_steps is None
                  else args.dark_steps)
    if dark_steps < 0:
        raise ValueError("Dark preparation steps cannot be negative.")
    for step in range(dark_steps):
        prior = sim.fields.psi.copy()
        _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=step*dt)
        if not result.converged:
            raise RuntimeError(f"Dark preparation failed at step {step+1}.")
        dark_changes.append(float(np.max(np.abs(sim.fields.psi-prior))))
    prepared = deepcopy(sim.fields)
    prepared_center = vortex_center(sim, scales, reference)
    if prepared_center is None:
        raise RuntimeError("No unique +1 vortex after dark preparation.")
    initial_beam = prepared_center+offset
    steps = int(config["timing"]["exposure_steps"] if args.max_steps is None
                else args.max_steps)
    if steps < 1:
        raise ValueError("Exposure steps must be positive.")
    sample_every = int(config["timing"]["sample_every"])
    frame_every = int(config["output"]["frame_every"])
    if sample_every < 1 or frame_every < 1:
        raise ValueError("Sample and frame intervals must be positive.")
    branches, frames = {}, []
    completed = {}
    termination = {}
    for mode in ("dynamic", "stationary", "dark"):
        sim.fields = deepcopy(prepared)
        sim.config.laser.enabled = mode != "dark"
        sim.config.laser.absorbed_power_W = float(config["absorbed_power_W"])
        sim.config.laser.sigma_m = float(config["spot_sigma_m"])
        beam = initial_beam.copy()
        target = beam.copy()
        last_center = prepared_center.copy()
        measurements = [prepared_center.copy()]
        rows = []
        reason = "completed"
        for step in range(steps+1):
            if step:
                start = (step-1)*dt
                if mode == "dynamic" and (step-1) % int(config["controller"]["update_every_steps"]) == 0:
                    current_center = vortex_center(sim, scales, last_center)
                    if current_center is None:
                        reason = "unique vortex center lost"
                        break
                    last_center = current_center
                    measurements.append(current_center.copy())
                    delay = int(config["controller"]["measurement_delay_updates"])
                    target = measurements[max(0, len(measurements)-1-delay)]+offset
                next_beam = (bounded_command(beam, target,
                    float(config["controller"]["max_speed_m_s"])*dt)
                    if mode == "dynamic" else beam.copy())
                if mode != "dark":
                    sim.config.laser.waypoints = [
                        LaserWaypointConfig(start, *beam),
                        LaserWaypointConfig(start+dt, *next_beam)]
                    sim.config.laser.validate()
                _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=start)
                if not result.converged:
                    raise RuntimeError(f"{mode} failed at step {step}.")
                beam = next_beam
            if step % sample_every == 0 or step == steps:
                row = _snapshot(sim, model, scales, reference, step*dt, step,
                                "exposure", mode)
                if mode != "dark" and step == 0:
                    row["laser_x_m"], row["laser_y_m"] = beam
                row["commanded_beam_x_m"] = beam[0] if mode != "dark" else None
                row["commanded_beam_y_m"] = beam[1] if mode != "dark" else None
                row["beam_vortex_separation_m"] = (
                    float(np.linalg.norm(beam-np.asarray([
                        row["vortex_x_zero_m"], row["vortex_y_zero_m"]])))
                    if mode != "dark" and row["complex_zero_valid"] else float("nan"))
                row["temperature_peak_beam_lag_m"] = (
                    float(np.linalg.norm(np.asarray([
                        row["electron_temperature_peak_x_m"],
                        row["electron_temperature_peak_y_m"]])-beam))
                    if mode != "dark" else float("nan"))
                rows.append(row)
                if mode == "dynamic" and (step % frame_every == 0 or step == steps):
                    frames.append(_frame(sim, row))
            if step and np.max(sim.fields.temperature) >= float(config["maximum_electron_temperature_K"]):
                if rows[-1]["step"] != step:
                    row = _snapshot(sim, model, scales, reference, step*dt, step,
                                    "thermal_cutoff", mode)
                    row["commanded_beam_x_m"] = beam[0] if mode != "dark" else None
                    row["commanded_beam_y_m"] = beam[1] if mode != "dark" else None
                    row["beam_vortex_separation_m"] = (
                        float(np.linalg.norm(beam-np.asarray([
                            row["vortex_x_zero_m"], row["vortex_y_zero_m"]])))
                        if mode != "dark" and row["complex_zero_valid"] else float("nan"))
                    row["temperature_peak_beam_lag_m"] = float("nan")
                    rows.append(row)
                    if mode == "dynamic":
                        frames.append(_frame(sim, row))
                reason = "electron-temperature cutoff"
                break
            if step and step % 50 == 0:
                print(f"{mode}: {step}/{steps}", flush=True)
        branches[mode] = rows
        completed[mode] = rows[-1]["step"]
        termination[mode] = reason
        _write_csv(output/f"{mode}_diagnostics.csv", rows)
    common_step = min(completed.values())
    comparable = {mode: max((row for row in rows if row["step"] <= common_step),
                            key=lambda row: row["step"])
                  for mode, rows in branches.items()}
    dark_y = comparable["dark"]["vortex_y_zero_m"]
    response = {
        "common_step": common_step,
        "dynamic_y_minus_dark_nm": (comparable["dynamic"]["vortex_y_zero_m"]-dark_y)*1e9,
        "stationary_y_minus_dark_nm": (comparable["stationary"]["vortex_y_zero_m"]-dark_y)*1e9,
        "dynamic_beam_vortex_distance_nm":
            comparable["dynamic"]["beam_vortex_separation_m"]*1e9,
        "dynamic_beam_displacement_nm":
            (comparable["dynamic"]["commanded_beam_y_m"]-initial_beam[1])*1e9,
        "dynamic_peak_electron_temperature_K":
            max(row["maximum_electron_temperature_K"] for row in branches["dynamic"]),
    }
    response["extra_motion_vs_stationary_nm"] = (
        response["dynamic_y_minus_dark_nm"]-response["stationary_y_minus_dark_nm"])
    _plot_motion(branches, output, offset, int(config["output"]["dpi"]))
    _plot_relative_motion(branches, output, int(config["output"]["dpi"]))
    _plot_film(frames, output, sim.mesh, int(config["output"]["dpi"]))
    zoom = (prepared_center/np.asarray([sim.mesh.dx, sim.mesh.dy]),
            int(config["output"]["zoom_radius_cells"]))
    _plot_film(frames, output, sim.mesh, int(config["output"]["dpi"]), zoom)
    if config["output"]["gif_enabled"]:
        _save_gif(frames, output, sim.mesh, int(config["output"]["gif_fps"]),
                  int(config["output"]["dpi"]))
        _save_gif(frames, output, sim.mesh, int(config["output"]["gif_fps"]),
                  int(config["output"]["dpi"]), zoom)
    summary = {"config": config, "effective_exposure_steps": steps,
               "effective_dark_steps": dark_steps,
               "dark_preparation_last_psi_change": dark_changes[-1] if dark_changes else None,
               "prepared_vortex_center_m": prepared_center.tolist(),
               "initial_beam_m": initial_beam.tolist(),
               "completed_steps": completed, "termination": termination,
               "response": response,
               "interpretation": "Feedback holds a commanded offset using the previous accepted vortex state; compare dynamic, stationary, and dark branches at the common step. Subcell motion and temperature-peak lag are exploratory; reaching the thermal cutoff is not sustained capture."}
    (output/"summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(output)
    return output


if __name__ == "__main__":
    main()
