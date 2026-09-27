"""Stateful zero -> positive -> negative TDGL Josephson protocol."""

from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from shs.config.builder import build_simulation
from shs.physics.josephson import junction_observables
from shs.solvers.coupled_solver import build_tdgl_model, coupled_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.utils.output import reserve_output_directory
from tdgl_diagnostics import apply_uniform_magnetic_field, _field_snapshot, _save_gif


def _current_parts(simulation):
    column = min(simulation.mesh.nx // 2, simulation.mesh.nx - 2)
    factor = simulation.mesh.dy * float(np.mean(simulation.material_map.thickness))
    return (float(np.trapezoid(simulation.fields.supercurrent_density_x[:, column]) * factor),
            float(np.trapezoid(simulation.fields.normal_current_density_x[:, column]) * factor))


def _stable(records, config):
    window = int(config["stability"]["window_steps"])
    if len(records) < max(window, int(config["stability"]["minimum_stage_steps"])): return False, "minimum duration"
    recent = records[-window:]
    phase = np.unwrap([r["junction_phase_difference_rad"] for r in recent])
    voltage = np.array([r["junction_voltage_V"] for r in recent])
    amplitude = np.array([r["junction_mean_amplitude"] for r in recent])
    phase_rate = np.diff(phase) / float(config["dt_s"])
    phase_stationary = np.max(np.abs(np.diff(phase))) <= float(config["stability"]["phase_step_tolerance_rad"])
    running_stable = len(phase_rate) > 2 and np.std(phase_rate) <= (
        float(config["stability"]["phase_rate_relative_tolerance"])
        * max(abs(np.mean(phase_rate)), 1.0))
    vscale = max(np.max(np.abs(voltage)), 1e-12)
    voltage_stable = np.ptp(voltage) <= float(config["stability"]["voltage_relative_tolerance"]) * vscale
    amplitude_stable = np.ptp(amplitude) <= float(config["stability"]["amplitude_tolerance"])
    return bool((phase_stationary or running_stable) and voltage_stable and amplitude_stable), ("stationary phase" if phase_stationary else "stable running phase")


def _plot(records, switches, path):
    t = np.array([r["time_s"] for r in records]) * 1e12
    fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True, constrained_layout=True)
    axes[0].step(t, np.array([r["commanded_current_A"] for r in records]) * 1e9, where="post")
    axes[0].set_ylabel("Command (nA)")
    axes[1].plot(t, [r["junction_phase_difference_rad"] for r in records]); axes[1].set_ylabel("Phase drop (rad)")
    axes[2].plot(t, np.array([r["junction_voltage_V"] for r in records]) * 1e6); axes[2].set_ylabel("Voltage (µV)")
    axes[3].plot(t, np.array([r["supercurrent_A"] for r in records]) * 1e9, label="super")
    axes[3].plot(t, np.array([r["normal_current_A"] for r in records]) * 1e9, label="normal")
    axes[3].set(xlabel="Time (ps)", ylabel="Current (nA)"); axes[3].legend()
    for axis in axes:
        axis.grid(alpha=.25)
        for switch in switches: axis.axvline(switch * 1e12, color="black", linestyle="--", alpha=.35)
    fig.savefig(path, dpi=160); plt.close(fig)


def run_case(config, magnitude, output):
    simulation = build_simulation(config["simulation_config"])
    apply_uniform_magnetic_field(simulation, float(config.get("applied_Bz_T", 0.0)), "symmetric")
    model = build_tdgl_model(simulation); scales = tdgl_scales(simulation, model)
    dt = float(config["dt_s"]); sample_every = int(config["output"]["sample_every"])
    stages = [("zero_relaxation", 0.0, int(config["zero_max_steps"])),
              ("positive_bias", magnitude, int(config["biased_max_steps"])),
              ("negative_bias", -magnitude, int(config["biased_max_steps"]))]
    records=[]; frames=[]; frame_rows=[]; switches=[]; time_s=0.0; global_step=0; stage_results=[]
    for stage_name, current, limit in stages:
        switches.append(time_s); stage_records=[]; reason="maximum steps"
        for stage_step in range(1, limit + 1):
            simulation.config.electrical.drive_mode="current"; simulation.config.current=current
            _, result=coupled_step(simulation, dt=dt, tdgl_model=model, time_s=time_s)
            if not result.converged: reason="coupling failure"; break
            time_s += dt; global_step += 1
            obs=junction_observables(simulation); super_i, normal_i=_current_parts(simulation)
            record={"step":global_step,"stage":stage_name,"stage_step":stage_step,"time_s":time_s,
                    "commanded_current_A":current,"supercurrent_A":super_i,"normal_current_A":normal_i,
                    "coupling_iterations":result.iterations,**obs}
            records.append(record); stage_records.append(record)
            if global_step % sample_every == 0 or stage_step == 1:
                row, frame=_field_snapshot(simulation,time_s,global_step,scales,1e-8,
                    {"applied_Bz_T":float(config.get("applied_Bz_T",0)),"calculate_self_field":True,
                     "observation_height_m":5e-9,"source_stride":int(config["output"].get("field_source_stride",2))},False,{"contour_inset_fraction":.2})
                frames.append(frame); frame_rows.append(row)
            stable, reason=_stable(stage_records,config)
            if stable: break
        stage_results.append({"stage":stage_name,"commanded_current_A":current,
                              "steps":len(stage_records),"stop_reason":reason,
                              "final":stage_records[-1] if stage_records else None})
    case=output/f"current_{magnitude*1e9:g}_nA"; case.mkdir(parents=True)
    with (case/"protocol.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=records[0].keys()); writer.writeheader(); writer.writerows(records)
    _plot(records,switches,case/"reversal_protocol.png")
    if config["output"]["gif_enabled"] and frames:
        _save_gif(frames,frame_rows,case/"order_parameter.gif",int(config["output"]["gif_fps"]),int(config["output"]["dpi"]))
    return {"current_magnitude_A":magnitude,"stages":stage_results,"case_directory":str(case)}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--config",type=Path,default=Path("tools/config/josephson_reversal_protocol.json")); args=parser.parse_args()
    config=json.loads(args.config.read_text(encoding="utf-8")); output=reserve_output_directory(config["output"]["directory"])
    summaries=[]
    for magnitude in map(float,config["current_magnitudes_A"]):
        print(f"Running ±{magnitude*1e9:g} nA",flush=True); summaries.append(run_case(config,magnitude,output))
    (output/"summary.json").write_text(json.dumps({"config":config,"cases":summaries},indent=2),encoding="utf-8")
    print(f"Protocol results written to {output}")


if __name__=="__main__": main()
