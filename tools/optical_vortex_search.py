"""Search TDGL laser protocols for controlled vortex motion, with paired controls."""
import argparse
import csv
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shs.config.builder import build_simulation
from shs.config.simulation import LaserWaypointConfig
from shs.solvers.coupled_solver import build_tdgl_model, coupled_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.utils.output import reserve_output_directory
from tdgl_diagnostics import (
    _field_snapshot, _seed_initial_state, apply_uniform_magnetic_field,
    track_vortices, _plot_final_state, _plot_current_and_field,
    _plot_timeseries, _plot_kymographs,
    _plot_vortex_trajectories, _plot_laser_vortex_distance, _save_gif,
)


def _write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def score_motion(trajectory, laser, start, direction, thresholds):
    """Require one same-sign continuous track; score following rather than proximity."""
    if len(trajectory) != len(laser) or len(laser) < 3:
        reason = "insufficient transport samples" if len(laser) < 3 else "vortex track lost"
        return {"qualified": False, "reason": reason, "track_fraction":
                len(trajectory)/max(len(laser), 1)}
    vortex = np.asarray(trajectory, dtype=float)
    beam = np.asarray(laser, dtype=float)
    unit = np.asarray(direction, dtype=float)
    unit /= np.linalg.norm(unit)
    displacement = float(np.dot(vortex[-1]-vortex[0], unit))
    beam_displacement = float(np.dot(beam[-1]-beam[0], unit))
    lag = np.linalg.norm(vortex-beam, axis=1)
    lag_rms = float(np.sqrt(np.mean(lag**2)))
    projection = (vortex-vortex[0])@unit
    commanded = (beam-beam[0])@unit
    correlation = (float(np.corrcoef(projection, commanded)[0, 1])
                   if np.std(projection) > 0 and np.std(commanded) > 0 else 0.0)
    qualified = bool(
        beam_displacement >= thresholds["minimum_beam_displacement_m"]
        and displacement >= thresholds["minimum_vortex_displacement_m"]
        and correlation >= thresholds["minimum_position_correlation"]
        and lag_rms <= thresholds["maximum_rms_lag_m"]
        and np.linalg.norm(vortex[0]-start) <= thresholds["maximum_initial_offset_m"]
    )
    return {"qualified": qualified, "reason": "kinematic threshold" if qualified else "following threshold failed",
            "track_fraction": 1.0, "displacement_along_path_m": displacement,
            "beam_displacement_m": beam_displacement,
            "position_correlation": correlation, "rms_lag_m": lag_rms,
            "final_lag_m": float(lag[-1])}


def controls_confirm(forward, controls, thresholds):
    """Separate commanded following from spontaneous or static-beam drift."""
    if not forward["qualified"] or not forward["converged"]:
        return False
    if not all(controls[mode]["converged"] for mode in
               ("no_laser", "stationary", "reverse")):
        return False
    minimum = thresholds["minimum_control_difference_m"]
    displacement = forward["displacement_along_path_m"]
    return bool(
        displacement > controls["no_laser"].get("displacement_along_path_m", 0) + minimum
        and displacement > controls["stationary"].get("displacement_along_path_m", 0) + minimum
        and controls["reverse"]["qualified"]
    )


def _run_protocol(config, field_t, power_w, speed_m_s, mode, output, label,
                  full_diagnostics=False):
    sim = build_simulation(config["simulation_config"])
    dt = sim.config.dt
    sim.config.current = 0.0
    sim.config.electrical.drive_mode = "current"
    sim.config.thermal.model = "two_temperature"
    sim.config.thermal.validate()
    path = config["path"]
    start = np.asarray(path["start_m"], dtype=float)
    axis = np.asarray(path["direction"], dtype=float)
    axis /= np.linalg.norm(axis)
    if mode == "reverse":
        axis = -axis
    distance = float(path["distance_m"])
    travel = distance/speed_m_s
    total = float(path["initial_hold_s"] + travel + path["final_hold_s"])
    steps = int(np.ceil(total/dt))
    if steps > int(config["search"]["maximum_steps_per_trial"]):
        raise ValueError(f"Trial requires {steps} steps; increase speed or step budget.")
    t0 = float(path["initial_hold_s"])
    t1 = t0+travel
    end = start + axis*distance
    sim.config.laser.enabled = mode != "no_laser"
    sim.config.laser.absorbed_power_W = 0.0 if mode == "no_laser" else power_w
    sim.config.laser.sigma_m = float(config["search"]["spot_sigma_m"])
    sim.config.laser.waypoints = [
        LaserWaypointConfig(0.0, *start),
        LaserWaypointConfig(t0, *start),
        LaserWaypointConfig(t1, *(start if mode == "stationary" else end)),
        LaserWaypointConfig(total, *(start if mode == "stationary" else end)),
    ]
    sim.config.laser.validate()
    apply_uniform_magnetic_field(sim, field_t, config["magnetic_field"]["gauge"])
    _seed_initial_state(sim, {"vortices": [config["initial_vortex"]]})
    model = build_tdgl_model(sim)
    scales = tdgl_scales(sim, model)
    rows, frames = [], []
    sample_every = (1 if full_diagnostics or travel/dt < 8
                    else int(config["search"]["sample_every"]))
    critical_temperature = float(np.min(sim.material_map.Tc))
    termination_reason = None
    for step in range(steps+1):
        if step:
            _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=(step-1)*dt)
            if not result.converged:
                termination_reason = "coupled solve did not converge"
                break
            if steps >= 200 and step % 50 == 0:
                print(f"  {label} {mode}: {step}/{steps} steps", flush=True)
        overheated = bool(np.max(sim.fields.temperature) >= critical_temperature)
        if step % sample_every == 0 or step == steps or overheated:
            row, frame = _field_snapshot(sim, step*dt, step, scales,
                                          float(config["search"]["vortex_amplitude_floor"]),
                                          None, False)
            rows.append(row)
            frames.append(frame)
        if overheated:
            termination_reason = "electron temperature reached Tc"
            break
    records = track_vortices(frames, rows, float(config["search"]["max_track_step_cells"]))
    _write_csv(output/f"{label}_diagnostics.csv", rows)
    _write_csv(output/f"{label}_vortex_trajectories.csv", records)
    if full_diagnostics:
        detail = output/"candidate_case"
        detail.mkdir(exist_ok=False)
        _write_csv(detail/"diagnostics.csv", rows)
        _write_csv(detail/"vortex_trajectories.csv", records)
        dpi = int(config["output"]["dpi"])
        _plot_final_state(frames[-1], detail/"final_state.png", dpi)
        _plot_current_and_field(frames[-1], detail/"current_and_field.png",
                                dpi, vector_stride=6)
        _plot_timeseries(rows, detail/"timeseries.png", dpi)
        _plot_kymographs(frames, rows, "vertical", detail/"kymographs.png", dpi)
        _plot_vortex_trajectories(frames, rows, records,
                                  detail/"vortex_trajectories.png", dpi)
        _plot_laser_vortex_distance(rows, detail/"laser_vortex_distance.png", dpi)
        fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
        for axis, key, title in zip(axes, ("temperature", "phonon_temperature", "amplitude"),
                                    ("Electrons (K)", "Phonons (K)", "Order amplitude")):
            image = axis.imshow(frames[-1][key], origin="lower", cmap="inferno")
            axis.set_title(title)
            fig.colorbar(image, ax=axis)
        fig.savefig(detail/"electron_phonon_film.png", dpi=dpi)
        plt.close(fig)
        if config["output"]["gif_enabled"]:
            _save_gif(frames, rows, detail/"order_parameter.gif",
                      int(config["output"]["gif_fps"]), dpi)
    first = [r for r in records if r["time_s"] <= t0
             and r["charge"] == int(config["initial_vortex"]["charge"])]
    if first:
        earliest = min(r["frame"] for r in first)
        seed = min((r for r in first if r["frame"] == earliest),
                   key=lambda r: np.linalg.norm(np.asarray([r["x_m"], r["y_m"]])-start))
        track = [r for r in records if r["track_id"] == seed["track_id"]]
    else:
        track = []
    by_frame = {r["frame"]: np.asarray([r["x_m"], r["y_m"]]) for r in track}
    selected = [i for i, row in enumerate(rows) if t0 <= row["time_s"] <= t1]
    vortex = [by_frame[i] for i in selected if i in by_frame]
    beam = [start + axis*min(max((rows[i]["time_s"]-t0)/travel, 0), 1)*distance
            if mode not in {"stationary", "no_laser"} else start for i in selected]
    result = score_motion(vortex, beam, start, axis, config["thresholds"])
    result.update({"mode": mode, "converged": len(rows) and rows[-1]["step"] == steps,
                   "completed_steps": rows[-1]["step"] if rows else 0,
                   "initial_vortex_count": rows[0]["total_vortices"] if rows else None,
                   "final_vortex_count": rows[-1]["total_vortices"] if rows else None,
                   "peak_electron_temperature_K": max(r["maximum_temperature_K"] for r in rows),
                   "peak_phonon_temperature_K": max(r["maximum_phonon_temperature_K"] for r in rows),
                   "trajectory_file": f"{label}_vortex_trajectories.csv"})
    result["termination_reason"] = termination_reason
    if result["peak_electron_temperature_K"] >= critical_temperature:
        result["qualified"] = False
        result["reason"] = "electron temperature reached Tc; vortex identity is not reliable"
    if len(vortex) != len(selected):
        result["qualified"] = False
        result["reason"] = "vortex track lost or charge changed"
    if any(rows[i]["total_vortices"] != 1 for i in selected):
        result["qualified"] = False
        result["reason"] = "additional or missing vortices during transport"
    return result


def _plot_summary(results, output):
    if not results:
        return
    fig, axis = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for field in sorted({r["field_T"] for r in results}):
        subset = [r for r in results if r["field_T"] == field]
        axis.scatter([r["speed_m_s"] for r in subset],
                     [r["forward"].get("displacement_along_path_m", 0)*1e9 for r in subset],
                     c=[r["power_W"] for r in subset], cmap="viridis", s=60,
                     label=f"B={field:g} T")
    axis.set(xscale="log", xlabel="Laser speed (m/s)",
             ylabel="Vortex displacement along commanded path (nm)",
             title="Optical vortex search; inspect controls before claiming following")
    axis.legend()
    fig.savefig(output/"search_displacement.png", dpi=160)
    plt.close(fig)


def _write_success_report(path, result):
    forward = result["detailed_run"]
    controls = result["controls"]
    lines = [
        "# Optical vortex-following candidate",
        "",
        f"Applied field: {result['field_T']:g} T; absorbed power: {result['power_W']:g} W; beam speed: {result['speed_m_s']:g} m/s.",
        f"Vortex displacement along path: {forward['displacement_along_path_m']*1e9:.3g} nm; position correlation: {forward['position_correlation']:.3g}; RMS beam lag: {forward['rms_lag_m']*1e9:.3g} nm.",
        f"Peak electron temperature: {forward['peak_electron_temperature_K']:.3g} K; peak phonon temperature: {forward['peak_phonon_temperature_K']:.3g} K.",
        "",
        "## Paired controls",
        "",
    ]
    for mode, control in controls.items():
        lines.append(f"- {mode}: {control.get('displacement_along_path_m', float('nan'))*1e9:.3g} nm along its commanded direction; converged={control['converged']}.")
    lines.extend(["", "This is a reproducible response within the configured model and window, not experimental validation of optical capture.", "",
                  "See the CSV files, figures, GIF, and summary.json here for the full trajectory and state diagnostics."])
    path.write_text("\n".join(lines)+"\n", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path("tools/config/optical_vortex_two_temperature_search.json"))
    parser.add_argument("--max-trials", type=int,
                        help="Override the configured trial budget for a short smoke run.")
    parser.add_argument("--resume", type=Path,
                        help="Continue an existing result directory from its saved summary.json.")
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if args.resume is None:
        output = reserve_output_directory(config["output_directory"])
        results = []
    else:
        output = args.resume
        previous = json.loads((output/"summary.json").read_text(encoding="utf-8"))
        if previous["config"] != config:
            raise ValueError("Resume config differs from the saved search configuration.")
        results = previous["results"]
    grid = list(itertools.product(
        config["search"]["magnetic_fields_T"],
        config["search"]["laser_powers_W"],
        config["search"]["laser_speeds_m_s"],
    ))
    # Seeded permutation spreads a bounded search across all three axes and
    # makes the same configuration reproducible.
    candidates = (grid[index] for index in np.random.default_rng(
        int(config["search"]["sampling_seed"])).permutation(len(grid)))
    for index, (field_t, power_w, speed) in enumerate(candidates):
        if results and results[-1]["correlated_motion"] and config["search"]["stop_at_first_correlation"]:
            break
        if index < len(results):
            continue
        if index >= (args.max_trials if args.max_trials is not None
                     else int(config["search"]["maximum_trials"])):
            break
        label = f"trial_{index+1:03d}"
        print(f"{label}: B={field_t:g} T, P={power_w:g} W, speed={speed:g} m/s", flush=True)
        forward = _run_protocol(config, field_t, power_w, speed, "moving", output, label)
        result = {"field_T": field_t, "power_W": power_w, "speed_m_s": speed,
                  "forward": forward, "controls": {}, "correlated_motion": False}
        if forward["qualified"] and forward["converged"]:
            for mode in ("no_laser", "stationary", "reverse"):
                result["controls"][mode] = _run_protocol(
                    config, field_t, power_w, speed, mode, output, f"{label}_{mode}")
            result["correlated_motion"] = controls_confirm(
                forward, result["controls"], config["thresholds"])
            if result["correlated_motion"]:
                detailed = _run_protocol(config, field_t, power_w, speed,
                                         "moving", output, f"{label}_detailed",
                                         full_diagnostics=True)
                result["detailed_run"] = detailed
                result["correlated_motion"] = bool(detailed["qualified"] and detailed["converged"])
                detail_directory = output/("successful_case" if result["correlated_motion"]
                                           else "candidate_case")
                if result["correlated_motion"]:
                    (output/"candidate_case").rename(detail_directory)
                (detail_directory/"summary.json").write_text(
                    json.dumps({"parameters": {"field_T": field_t,
                          "power_W": power_w, "speed_m_s": speed},
                          "search_result": result, "search_config": config},
                          indent=2), encoding="utf-8")
                if result["correlated_motion"]:
                    _write_success_report(detail_directory/"report.md", result)
        results.append(result)
        (output/"summary.json").write_text(json.dumps({"config": config,
            "results": results}, indent=2), encoding="utf-8")
        if result["correlated_motion"] and config["search"]["stop_at_first_correlation"]:
            break
    _plot_summary(results, output)
    print(f"Search written to {output}; correlated trials: "
          f"{sum(r['correlated_motion'] for r in results)}")
    return output


if __name__ == "__main__":
    main()
