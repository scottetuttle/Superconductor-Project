"""Run an inexpensive optical-manipulation force-balance experiment."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from shs.physics.reduced_vortex import (
    ReducedPinningSite, ReducedVortexMaterial, ReducedVortexState,
    gaussian_temperature, reduced_vortex_step,
)
from shs.utils.output import reserve_output_directory


DEFAULT_CONFIG = Path(__file__).parent / "config" / "reduced_optical_tweezers.json"


def protocol_value(protocol, time_s):
    """Linearly interpolate laser position and requested peak temperature rise."""
    times = np.asarray([point["time_s"] for point in protocol], dtype=float)
    index = int(np.searchsorted(times, time_s, side="right"))
    if index == 0:
        point = protocol[0]
        return np.array([point["x_m"], point["y_m"]]), float(point["delta_temperature_K"])
    if index >= len(protocol):
        point = protocol[-1]
        return np.array([point["x_m"], point["y_m"]]), float(point["delta_temperature_K"])
    before, after = protocol[index - 1], protocol[index]
    fraction = (time_s - before["time_s"]) / (after["time_s"] - before["time_s"])
    center = np.array([before["x_m"], before["y_m"]], dtype=float) + fraction * (
        np.array([after["x_m"], after["y_m"]], dtype=float)
        - np.array([before["x_m"], before["y_m"]], dtype=float)
    )
    delta = before["delta_temperature_K"] + fraction * (
        after["delta_temperature_K"] - before["delta_temperature_K"]
    )
    return center, float(delta)


def load_config(path):
    with Path(path).open(encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("physics_backend") != "reduced_vortex":
        raise ValueError("physics_backend must be reduced_vortex for this runner.")
    protocol = config["laser"]["protocol"]
    times = [float(point["time_s"]) for point in protocol]
    if not protocol or any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("Laser protocol times must be nonempty and strictly increasing.")
    if float(config["run"]["dt_s"]) <= 0 or int(config["run"]["sample_every"]) < 1:
        raise ValueError("Run timestep and sample interval must be positive.")
    critical = float(config["material"]["critical_temperature_K"])
    bath = float(config["bath_temperature_K"])
    if bath < 0 or any(float(point["delta_temperature_K"]) < 0 for point in protocol):
        raise ValueError("Bath temperature and laser temperature rises must be nonnegative.")
    if max(bath + float(point["delta_temperature_K"]) for point in protocol) >= critical:
        raise ValueError(
            "Reduced-vortex backend requires the prescribed hotspot to remain below Tc; "
            "use full TDGL for moving normal regions."
        )
    for site in config.get("pinning_sites", []):
        if float(site["energy_per_length_J_per_m"]) < 0 or float(site["sigma_m"]) <= 0:
            raise ValueError("Pinning energy must be nonnegative and sigma_m positive.")
    return config


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)


def run(config):
    output = reserve_output_directory(config["output"]["directory"])
    material = ReducedVortexMaterial(**config["material"])
    material.validate()
    pins = [ReducedPinningSite(**site) for site in config.get("pinning_sites", [])]
    vortices = config["initial_vortices"]
    first_center, _ = protocol_value(config["laser"]["protocol"], 0.0)
    state = ReducedVortexState(
        positions_m=np.asarray([[item["x_m"], item["y_m"]] for item in vortices]),
        charges=np.asarray([item["charge"] for item in vortices]),
        thermal_center_m=first_center.copy(), thermal_delta_K=0.0,
    )
    run_config, laser = config["run"], config["laser"]
    dt = float(run_config["dt_s"])
    duration = float(run_config["duration_s"])
    steps = int(np.ceil(duration / dt))
    sample_every = int(run_config["sample_every"])
    bounds = config["domain"]
    bounds_tuple = (0.0, bounds["width_m"], 0.0, bounds["height_m"])
    rows = []
    limiter_activations = 0
    maximum_unlimited_displacement = 0.0

    def sample(step, force_components):
        time_s = step * dt
        commanded, commanded_delta = protocol_value(laser["protocol"], time_s)
        for index, (position, charge) in enumerate(zip(state.positions_m, state.charges)):
            components = force_components or {
                name: np.zeros_like(state.positions_m)
                for name in ("thermal", "pinning", "interaction", "lorentz")
            }
            rows.append({
                "step": step, "time_s": time_s, "vortex_id": index + 1,
                "charge": int(charge), "x_m": position[0], "y_m": position[1],
                "commanded_laser_x_m": commanded[0], "commanded_laser_y_m": commanded[1],
                "thermal_center_x_m": state.thermal_center_m[0],
                "thermal_center_y_m": state.thermal_center_m[1],
                "commanded_delta_temperature_K": commanded_delta,
                "thermal_delta_temperature_K": state.thermal_delta_K,
                "local_temperature_K": float(gaussian_temperature(
                    position, state.thermal_center_m, config["bath_temperature_K"],
                    state.thermal_delta_K, laser["sigma_m"]
                )),
                "laser_vortex_distance_m": float(np.linalg.norm(position-state.thermal_center_m)),
                **{f"{name}_force_N_per_m": float(np.linalg.norm(values[index]))
                   for name, values in components.items()},
            })

    sample(0, None)
    force_components = None
    for step in range(1, steps + 1):
        time_s = min(step * dt, duration)
        commanded, delta = protocol_value(laser["protocol"], time_s)
        force, force_components = reduced_vortex_step(
            state, material, min(dt, duration-(step-1)*dt),
            commanded_center_m=commanded, commanded_delta_K=delta,
            thermal_sigma_m=float(laser["sigma_m"]),
            thermal_response_time_s=float(laser["thermal_response_time_s"]),
            pinning_sites=pins,
            current_density_A_per_m2=config.get("current_density_A_per_m2", [0.0, 0.0]),
            bounds_m=bounds_tuple,
            maximum_displacement_m=float(run_config["maximum_displacement_per_step_m"]),
        )
        raw_displacement = float(np.max(
            np.linalg.norm(force, axis=1) * min(dt, duration-(step-1)*dt)
            / material.viscosity_N_s_per_m2
        ))
        maximum_unlimited_displacement = max(maximum_unlimited_displacement, raw_displacement)
        limiter_activations += int(
            raw_displacement > float(run_config["maximum_displacement_per_step_m"])
        )
        if step % sample_every == 0 or step == steps:
            sample(step, force_components)

    _write_csv(output / "reduced_vortex_trajectories.csv", rows)
    dpi = int(config["output"].get("dpi", 150))
    fig, axis = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for vortex_id in sorted({row["vortex_id"] for row in rows}):
        track = [row for row in rows if row["vortex_id"] == vortex_id]
        axis.plot(np.asarray([r["x_m"] for r in track])*1e6,
                  np.asarray([r["y_m"] for r in track])*1e6, "-", linewidth=2,
                  label=f"vortex {vortex_id}")
    axis.plot(np.asarray([point["x_m"] for point in laser["protocol"]])*1e6,
              np.asarray([point["y_m"] for point in laser["protocol"]])*1e6,
              "--", color="gold", linewidth=2, label="commanded laser")
    axis.scatter([site.x_m*1e6 for site in pins], [site.y_m*1e6 for site in pins],
                 marker="x", color="black", label="pinning sites")
    axis.set(xlabel="x (micrometres)", ylabel="y (micrometres)",
             xlim=(0, bounds["width_m"]*1e6), ylim=(0, bounds["height_m"]*1e6),
             title="Reduced vortex manipulation trajectory")
    axis.legend(); axis.set_aspect("equal")
    fig.savefig(output / "reduced_vortex_trajectory.png", dpi=dpi); plt.close(fig)

    first = rows[0]; last = rows[-1]
    evaluation = config["evaluation"]
    transport_start, transport_end = map(float, evaluation["transport_interval_s"])
    transport_rows = [
        row for row in rows if transport_start <= row["time_s"] <= transport_end
    ]
    maximum_transport_lag = max(
        row["laser_vortex_distance_m"] for row in transport_rows
    )
    destination = np.asarray(evaluation["destination_m"], dtype=float)
    final_destination_distance = float(np.linalg.norm(
        np.asarray([last["x_m"], last["y_m"]]) - destination
    ))
    screened_success = (
        limiter_activations == 0
        and maximum_transport_lag <= float(evaluation["maximum_following_lag_m"])
        and final_destination_distance <= float(evaluation["delivery_radius_m"])
    )
    summary = {
        "physics_backend": "reduced_vortex",
        "output_directory": str(output),
        "steps": steps,
        "viscosity_N_s_per_m2": material.viscosity_N_s_per_m2,
        "thermal_force_coefficient_N_per_K": material.thermal_force_coefficient_N_per_K,
        "initial_position_m": [first["x_m"], first["y_m"]],
        "final_position_m": [last["x_m"], last["y_m"]],
        "net_displacement_m": float(np.hypot(last["x_m"]-first["x_m"], last["y_m"]-first["y_m"])),
        "final_laser_vortex_distance_m": last["laser_vortex_distance_m"],
        "maximum_transport_lag_m": maximum_transport_lag,
        "final_destination_distance_m": final_destination_distance,
        "maximum_unlimited_step_displacement_m": maximum_unlimited_displacement,
        "motion_limiter_activations": limiter_activations,
        "screened_manipulation_success": screened_success,
        "success_criteria": evaluation,
        "interpretation": "Screening result only; full TDGL and calibrated thermal/pinning inputs are required.",
    }
    with (output / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args(argv)
    summary = run(load_config(args.config))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
