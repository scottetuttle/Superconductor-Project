"""Simulate current- or voltage-biased lumped SIS RCSJ junctions."""
import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from shs.physics.josephson_rcsj import (RCSJJunction, ambegaokar_baratoff_critical_current, extended_junction_interference,
    current_biased_rhs, junction_components)
from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM
from shs.utils.output import reserve_output_directory


def _junction(c):
    j=c["junction"]; ic=j.get("critical_current_A")
    if ic is None: ic=ambegaokar_baratoff_critical_current(c["temperature_K"],c["critical_temperature_K"],j["normal_resistance_ohm"])
    amplitude,offset=1.0,0.0
    extended=c.get("extended_junction")
    if extended:
        factor,_,_=extended_junction_interference(float(extended["width_m"]),float(extended.get("junction_x_m",0)),float(extended.get("magnetic_field_T",0)),float(extended.get("effective_magnetic_length_m",0)),extended.get("vortices",[]),int(extended.get("samples",401)))
        amplitude,offset=abs(factor),float(np.angle(factor))
    return RCSJJunction(float(ic),float(j["normal_resistance_ohm"]),float(j["capacitance_F"]),float(j.get("magnetic_flux_quanta",0))*SUPERCONDUCTING_FLUX_QUANTUM,amplitude,offset)


def _drive(t, drive):
    return float(drive["dc"])+float(drive.get("ac_amplitude",0))*np.sin(2*np.pi*float(drive.get("frequency_Hz",0))*t)


def simulate(c, dc):
    j=_junction(c); drive=dict(c["drive"]); drive["dc"]=dc; duration=float(c["run"]["duration_s"]); samples=int(c["run"]["samples"])
    times=np.linspace(0,duration,samples)
    configured_max_step=float(c["run"]["max_step_s"])
    frequency=abs(float(drive.get("frequency_Hz",0)))
    steps_per_cycle=int(c["run"].get("microwave_steps_per_cycle",80))
    max_step=min(configured_max_step,1/(frequency*steps_per_cycle)) if frequency and drive.get("ac_amplitude",0) else configured_max_step
    if drive["mode"]=="current":
        def rhs(t,y):
            first,second=current_biased_rhs(j,y[0],y[1],_drive(t,drive)); return [first,second]
        solution=solve_ivp(rhs,(0,duration),[float(c["initial_phase_rad"]),0.0],t_eval=times,
                           rtol=float(c["run"]["rtol"]),atol=float(c["run"]["atol"]),max_step=max_step)
        phase,rate=solution.y; bias=np.array([_drive(t,drive) for t in times])
        acceleration=np.gradient(rate,times)
    else:
        voltage=np.array([_drive(t,drive) for t in times]); rate=2*np.pi*voltage/SUPERCONDUCTING_FLUX_QUANTUM
        phase=float(c["initial_phase_rad"])+np.concatenate(([0],np.cumsum((rate[1:]+rate[:-1])*np.diff(times)/2)))
        acceleration=np.gradient(rate,times); bias=np.empty_like(times)
    rows=[]
    for i,t in enumerate(times):
        parts=junction_components(j,phase[i],rate[i],acceleration[i])
        if drive["mode"]=="voltage": bias[i]=parts["supercurrent_A"]+parts["resistive_current_A"]+parts["capacitive_current_A"]
        rows.append({"time_s":t,"bias_current_A":bias[i],**parts})
    return j,rows


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--config",type=Path,default=Path("tools/config/sis_junction_sweep.json")); args=parser.parse_args(); c=json.loads(args.config.read_text())
    out=reserve_output_directory(c["output_directory"]); cases=[]
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for dc in c["sweep_dc_values"]:
        j,rows=simulate(c,float(dc)); label=f"{float(dc)*1e6:g} µ{'A' if c['drive']['mode']=='current' else 'V'}"
        with (out/f"case_{float(dc):.6g}.csv").open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
        t=np.array([r["time_s"] for r in rows]); tail=rows[len(rows)//2:]
        axes[0,0].plot(t*1e9,[r["phase_rad"] for r in rows],label=label); axes[0,1].plot(t*1e9,np.array([r["voltage_V"] for r in rows])*1e6,label=label)
        axes[1,0].plot(t*1e9,np.array([r["supercurrent_A"] for r in rows])*1e6,label=label)
        cases.append({"dc":dc,"mean_voltage_V":float(np.mean([r["voltage_V"] for r in tail])),"mean_current_A":float(np.mean([r["bias_current_A"] for r in tail]))})
    axes[1,1].plot(np.array([x["mean_current_A"] for x in cases])*1e6,np.array([x["mean_voltage_V"] for x in cases])*1e6,"o-")
    axes[0,0].set(ylabel="Phase (rad)"); axes[0,1].set(ylabel="Voltage (µV)"); axes[1,0].set(ylabel="Josephson current (µA)"); axes[1,1].set(xlabel="Mean current (µA)",ylabel="Mean voltage (µV)",title="I–V")
    for a in axes.flat: a.grid(alpha=.25)
    axes[0,0].legend(fontsize=7); fig.savefig(out/"sis_diagnostics.png",dpi=160); plt.close(fig)
    (out/"summary.json").write_text(json.dumps({"junction":j.__dict__,"cases":cases},indent=2)); print(f"SIS results written to {out}")


if __name__=="__main__": main()
