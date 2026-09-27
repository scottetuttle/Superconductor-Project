"""Prepare a simple centered-vortex, stationary-laser AMR reference state."""

import argparse
import json
from dataclasses import fields as dataclass_fields
from pathlib import Path

import numpy as np

from shs.config.builder import build_simulation
from shs.config.simulation import LaserWaypointConfig
from shs.solvers.coupled_solver import build_tdgl_model, coupled_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.utils.output import reserve_output_directory
from tdgl_diagnostics import apply_uniform_magnetic_field, _seed_initial_state
from optical_feedback_tweezer import vortex_center


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path("tools/config/nested_mesh_simple_reference.json"))
    parser.add_argument("--dark-steps", type=int)
    parser.add_argument("--exposure-steps", type=int)
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    output = reserve_output_directory(config["output_directory"])
    sim = build_simulation(config["simulation_config"])
    sim.config.current = 0.0
    sim.config.electrical.drive_mode = "current"
    sim.config.thermal.model = "two_temperature"
    sim.config.pinning.enabled = False
    sim.config.laser.enabled = False
    apply_uniform_magnetic_field(sim, float(config["applied_field_T"]))
    _seed_initial_state(sim, {"vortices": [config["initial_vortex"]]})
    model = build_tdgl_model(sim)
    scales = tdgl_scales(sim, model)
    dt = sim.config.dt
    dark_steps = int(config["timing"]["dark_steps"] if args.dark_steps is None else args.dark_steps)
    exposure_steps = int(config["timing"]["exposure_steps"]
                         if args.exposure_steps is None else args.exposure_steps)
    if dark_steps < 0 or exposure_steps < 0:
        raise ValueError("Step counts cannot be negative.")
    for step in range(dark_steps):
        _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=step*dt)
        if not result.converged:
            raise RuntimeError(f"Dark preparation failed at step {step+1}.")
    center = vortex_center(sim, scales, np.asarray(config["seed_position_m"]))
    if center is None:
        raise RuntimeError("No unique centered vortex after dark preparation.")
    position = np.asarray(config["laser"]["position_m"], dtype=float)
    sim.config.laser.enabled = True
    sim.config.laser.absorbed_power_W = float(config["laser"]["absorbed_power_W"])
    sim.config.laser.sigma_m = float(config["laser"]["spot_sigma_m"])
    sim.config.laser.waypoints = [LaserWaypointConfig(0.0, *position),
                                  LaserWaypointConfig(max(dt, exposure_steps*dt), *position)]
    sim.config.laser.validate()
    for step in range(exposure_steps):
        _, result = coupled_step(sim, dt=dt, tdgl_model=model,
                                 time_s=(dark_steps+step)*dt)
        if not result.converged:
            raise RuntimeError(f"Exposure failed at step {step+1}.")
    center = vortex_center(sim, scales, center)
    if center is None:
        raise RuntimeError("Unique vortex lost during reference preparation.")
    arrays = {item.name: getattr(sim.fields, item.name) for item in dataclass_fields(sim.fields)
              if isinstance(getattr(sim.fields, item.name), np.ndarray)}
    checkpoint = f"checkpoint_{dark_steps+exposure_steps:08d}.npz"
    np.savez_compressed(output/checkpoint,
                        checkpoint_step=np.asarray(dark_steps+exposure_steps), **arrays)
    state = {"step": dark_steps+exposure_steps, "dt_s": dt,
             "last_center_m": center.tolist(), "checkpoint_file": checkpoint,
             "dark_steps": dark_steps, "exposure_steps": exposure_steps,
             "peak_electron_temperature_K": float(np.max(sim.fields.temperature))}
    (output/"state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")
    (output/"config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    print(output)
    return output


if __name__ == "__main__":
    main()
