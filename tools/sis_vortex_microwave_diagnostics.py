"""Sweep vortex layouts and microwave amplitudes for a finite-width SIS junction."""
import argparse, copy, csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from sis_junction_diagnostics import simulate, _junction
from shs.physics.josephson_rcsj import extended_junction_interference
from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM
from shs.utils.output import reserve_output_directory

def film_phase_map(config,vortices,global_phase=0.0,nx=240,ny=160):
    ext=config["extended_junction"]; width=float(ext["width_m"]); length=float(ext.get("electrode_length_m",width)); x=np.linspace(-length,length,nx); y=np.linspace(-width/2,width/2,ny); xx,yy=np.meshgrid(x,y); phase=np.where(xx>=0,global_phase,0.0)
    for v in vortices:
        contribution=int(v["charge"])*np.arctan2(yy-float(v["y_m"]),xx-float(v["x_m"])); side=(xx<0) if v.get("electrode","right")=="left" else (xx>=0); phase+=np.where(side,contribution,0.0)
    return x,y,np.angle(np.exp(1j*phase))

def save_film_image(config,scenario,path):
    x,y,phase=film_phase_map(config,scenario["vortices"]); fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True); im=ax.imshow(phase,origin="lower",extent=[x[0]*1e6,x[-1]*1e6,y[0]*1e6,y[-1]*1e6],cmap="twilight",vmin=-np.pi,vmax=np.pi,aspect="auto"); ax.axvline(0,color="white",linewidth=4,label="insulating barrier")
    for v in scenario["vortices"]: ax.scatter(float(v["x_m"])*1e6,float(v["y_m"])*1e6,marker="o" if int(v["charge"])>0 else "x",s=80,label=f"q={v['charge']} {v.get('electrode','right')}")
    ax.set(xlabel="x (µm)",ylabel="y (µm)",title=f"Electrode phase texture: {scenario['name']}"); ax.legend(fontsize=8); fig.colorbar(im,ax=ax,label="Phase (rad)"); fig.savefig(path,dpi=160); plt.close(fig)

def save_gif(config,scenario,rows,path):
    frequency=float(config["drive"]["frequency_Hz"]); period=1/frequency; cycles=float(config["output"].get("gif_cycles",2)); frame_count=int(config["output"].get("gif_frames_per_cycle",48)*cycles); end=rows[-1]["time_s"]; targets=np.linspace(max(rows[0]["time_s"],end-cycles*period),end,frame_count); times=np.array([r["time_s"] for r in rows]); selected=[rows[min(np.searchsorted(times,t),len(rows)-1)] for t in targets]; fig,(fa,ja)=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True); x,y,phase=film_phase_map(config,scenario["vortices"],selected[0]["phase_rad"]); image=fa.imshow(phase,origin="lower",extent=[x[0]*1e6,x[-1]*1e6,y[0]*1e6,y[-1]*1e6],cmap="twilight",vmin=-np.pi,vmax=np.pi,aspect="auto"); fa.axvline(0,color="white",linewidth=3)
    _,ay,theta=extended_junction_interference(float(config["extended_junction"]["width_m"]),0,float(config["extended_junction"]["magnetic_field_T"]),float(config["extended_junction"]["effective_magnetic_length_m"]),scenario["vortices"],int(config["extended_junction"]["samples"])); line,=ja.plot(ay*1e6,np.sin(theta+selected[0]["phase_rad"])); ja.set_ylim(-1.1,1.1); ja.set(xlabel="Junction y (µm)",ylabel="Normalized local Josephson current")
    def update(i):
        row=selected[i]; image.set_data(film_phase_map(config,scenario["vortices"],row["phase_rad"])[2]); line.set_ydata(np.sin(theta+row["phase_rad"])); fa.set_title(f"t={row['time_s']*1e9:.3f} ns, V={row['voltage_V']*1e6:.2f} µV"); return image,line
    animation.FuncAnimation(fig,update,frames=len(selected)).save(path,writer=animation.PillowWriter(fps=int(config["output"].get("gif_fps",12))),dpi=110); plt.close(fig)

def save_dynamic_diagnostics(config, scenario, rows, prefix):
    """Resolve drive, Josephson, plasma, and spatial-current dynamics."""
    tail = rows[int(len(rows) * float(config.get("steady_state_fraction", 0.5))):]
    time = np.asarray([row["time_s"] for row in tail])
    phase = np.unwrap(np.asarray([row["phase_rad"] for row in tail]))
    phase_rate = np.asarray([row["phase_rate_rad_per_s"] for row in tail])
    voltage = np.asarray([row["voltage_V"] for row in tail])
    supercurrent = np.asarray([row["supercurrent_A"] for row in tail])
    resistive = np.asarray([row["resistive_current_A"] for row in tail])
    capacitive = np.asarray([row["capacitive_current_A"] for row in tail])
    drive_frequency = float(config["drive"]["frequency_Hz"])
    dt = float(np.mean(np.diff(time)))
    frequencies = np.fft.rfftfreq(len(time), dt)
    window = np.hanning(len(time))
    _, junction_y, vortex_theta = extended_junction_interference(
        float(config["extended_junction"]["width_m"]),
        float(config["extended_junction"].get("junction_x_m", 0)),
        float(config["extended_junction"]["magnetic_field_T"]),
        float(config["extended_junction"]["effective_magnetic_length_m"]),
        scenario["vortices"], int(config["extended_junction"]["samples"]),
    )
    center_local_current = np.sin(vortex_theta[len(vortex_theta)//2] + phase)

    def normalized_spectrum(values):
        amplitude = np.abs(np.fft.rfft((values-np.mean(values))*window))
        return amplitude/max(float(np.max(amplitude)), np.finfo(float).tiny)

    component_spectra = {
        "voltage": normalized_spectrum(voltage),
        "phase_rate": normalized_spectrum(phase_rate),
        "josephson_current": normalized_spectrum(supercurrent),
        "resistive_current": normalized_spectrum(resistive),
        "capacitive_current": normalized_spectrum(capacitive),
        "center_local_josephson_current": normalized_spectrum(center_local_current),
    }
    spectrum = component_spectra["voltage"]
    spectral_rows = [
        {"frequency_Hz": frequencies[index], **{
            f"normalized_{name}_amplitude": values[index]
            for name, values in component_spectra.items()
        }} for index in range(len(frequencies))
    ]
    _write_rows(prefix.with_name(prefix.name + "_transient.csv"), tail)
    _write_rows(prefix.with_name(prefix.name + "_spectrum.csv"), spectral_rows)

    junction = _junction({**config, "extended_junction": {
        **config["extended_junction"], "vortices": scenario["vortices"]
    }})
    plasma_frequency = np.sqrt(
        2*np.pi*abs(junction.effective_critical_current_A)
        / (SUPERCONDUCTING_FLUX_QUANTUM*junction.capacitance_F)
    )/(2*np.pi) if junction.capacitance_F else np.nan
    beta_c = (2*np.pi*abs(junction.effective_critical_current_A)
              * junction.normal_resistance_ohm**2 * junction.capacitance_F
              / SUPERCONDUCTING_FLUX_QUANTUM)
    josephson_frequency = float(np.mean(phase_rate)/(2*np.pi))
    def dominant(values):
        peaks = np.argsort(values[1:])[-12:][::-1] + 1
        return [float(frequencies[index]) for index in peaks[:8]]

    frequency_summary = {
        "drive_frequency_Hz": drive_frequency,
        "mean_josephson_frequency_Hz": josephson_frequency,
        "josephson_to_drive_ratio": josephson_frequency/drive_frequency,
        "small_signal_plasma_frequency_Hz": float(plasma_frequency),
        "mccumber_parameter": float(beta_c),
        "dominant_frequencies_Hz": {
            name: dominant(values) for name, values in component_spectra.items()
        },
        "phase_rate_resistive_current_correlation": float(np.corrcoef(
            phase_rate, resistive
        )[0, 1]),
    }
    prefix.with_name(prefix.name + "_frequency_summary.json").write_text(
        json.dumps(frequency_summary, indent=2), encoding="utf-8"
    )

    plot_cycles = float(config["output"].get("dynamic_plot_cycles", 5.0))
    plot_start = time[-1] - plot_cycles/drive_frequency
    plot_mask = time >= plot_start
    plot_time = time[plot_mask]
    plot_phase = phase[plot_mask]
    plot_phase_rate = phase_rate[plot_mask]
    plot_supercurrent = supercurrent[plot_mask]
    plot_resistive = resistive[plot_mask]
    plot_capacitive = capacitive[plot_mask]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    relative_time = (plot_time-plot_time[0])*1e9
    axes[0, 0].plot(relative_time, plot_phase/(2*np.pi), label="unwrapped")
    axes[0, 0].plot(relative_time, np.mod(plot_phase, 2*np.pi)/(2*np.pi), alpha=.7,
                    label="wrapped")
    axes[0, 0].set(xlabel="Steady-state time (ns)", ylabel="Phase / 2π",
                   title="Global junction phase")
    axes[0, 0].legend()
    axes[0, 1].plot(relative_time, plot_phase_rate/(2*np.pi*drive_frequency))
    axes[0, 1].axhline(josephson_frequency/drive_frequency, color="black",
                       linestyle="--", label="time average")
    axes[0, 1].set(xlabel="Steady-state time (ns)", ylabel="Instantaneous fJ / fdrive",
                   title="Nonuniform phase rotation")
    axes[0, 1].legend()
    axes[1, 0].plot(relative_time, plot_supercurrent*1e6, label="Josephson")
    axes[1, 0].plot(relative_time, plot_resistive*1e6, label="resistive")
    axes[1, 0].plot(relative_time, plot_capacitive*1e6, label="capacitive")
    axes[1, 0].set(xlabel="Steady-state time (ns)", ylabel="Current (µA)",
                   title="RCSJ current balance")
    axes[1, 0].legend()
    positive = frequencies > 0
    axes[1, 1].semilogy(frequencies[positive]/1e9,
                        np.maximum(spectrum[positive], 1e-12))
    axes[1, 1].axvline(drive_frequency/1e9, color="black", linestyle="--",
                       label="microwave")
    if np.isfinite(plasma_frequency):
        axes[1, 1].axvline(plasma_frequency/1e9, color="tab:red", linestyle=":",
                           label="small-signal plasma")
    axes[1, 1].set(xlim=(0, min(frequencies[-1], 8*drive_frequency)/1e9),
                   xlabel="Frequency (GHz)", ylabel="Normalized voltage FFT",
                   title="Voltage spectrum")
    axes[1, 1].legend()
    for axis in axes.flat: axis.grid(alpha=.25)
    fig.savefig(prefix.with_name(prefix.name + "_dynamics.png"), dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True,
                             constrained_layout=True)
    groups = [
        (("phase_rate", "voltage", "resistive_current"), "Direct Josephson relation"),
        (("josephson_current", "center_local_josephson_current"),
         "Josephson-current harmonics"),
        (("capacitive_current",), "Capacitive response"),
    ]
    positive = frequencies > 0
    for axis, (names, title) in zip(axes, groups):
        for name in names:
            axis.semilogy(frequencies[positive]/1e9,
                          np.maximum(component_spectra[name][positive], 1e-12),
                          label=name.replace("_", " "))
        axis.set(ylabel="Normalized FFT", title=title,
                 xlim=(0, min(frequencies[-1], 8*drive_frequency)/1e9))
        axis.grid(alpha=.25); axis.legend()
    axes[-1].set_xlabel("Frequency (GHz)")
    fig.savefig(prefix.with_name(prefix.name + "_component_spectra.png"), dpi=170)
    plt.close(fig)

    stride = max(1, len(plot_time)//1200)
    local_current = np.sin(vortex_theta[:, None] + plot_phase[None, ::stride])
    fig, axis = plt.subplots(figsize=(10, 5), constrained_layout=True)
    image = axis.imshow(local_current, origin="lower", aspect="auto", cmap="RdBu_r",
                        vmin=-1, vmax=1, extent=[relative_time[0], relative_time[-1],
                        junction_y[0]*1e6, junction_y[-1]*1e6])
    axis.set(xlabel="Steady-state time (ns)", ylabel="Junction coordinate y (µm)",
             title="Normalized local Josephson current: sin[θvortex(y) + φ(t)]")
    fig.colorbar(image, ax=axis, label="Local current / local Ic")
    fig.savefig(prefix.with_name(prefix.name + "_local_current_kymograph.png"), dpi=170)
    plt.close(fig)
    return frequency_summary

def _write_rows(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--config",type=Path,default=Path("tools/config/sis_vortex_shapiro.json")); args=parser.parse_args(); c=json.loads(args.config.read_text()); out=reserve_output_directory(c["output_directory"]); quantum=SUPERCONDUCTING_FLUX_QUANTUM*float(c["drive"]["frequency_Hz"]); amplitudes=list(map(float,c["microwave_amplitudes_A"])); biases=list(map(float,c["sweep_dc_values"])); summaries=[]
    fig,axes=plt.subplots(len(c["vortex_scenarios"]),len(amplitudes),figsize=(4*len(amplitudes),3.2*len(c["vortex_scenarios"])),squeeze=False,constrained_layout=True)
    for si,scenario in enumerate(c["vortex_scenarios"]):
        save_film_image(c,scenario,out/f"{scenario['name']}_film_phase.png"); ss={"name":scenario["name"],"microwave_sweeps":[]}
        for ai,amp in enumerate(amplitudes):
            sc=copy.deepcopy(c); sc["extended_junction"]["vortices"]=scenario["vortices"]; sc["drive"]["ac_amplitude"]=amp; junction=_junction(sc); points=[]; gif_rows=None
            for bias in biases:
                _,rows=simulate(sc,bias); tail=rows[int(len(rows)*float(c.get("steady_state_fraction",0.5))):]; voltages=np.array([r["voltage_V"] for r in tail]); voltage=float(np.mean(voltages)); normalized=voltage/quantum; nearest=int(np.rint(normalized)); phase=np.unwrap([r["phase_rad"] for r in tail]); duration=tail[-1]["time_s"]-tail[0]["time_s"]; phase_lock=(phase[-1]-phase[0])/(2*np.pi*float(c["drive"]["frequency_Hz"])*duration); points.append({"scenario":scenario["name"],"microwave_amplitude_A":amp,"bias_current_A":bias,"mean_voltage_V":voltage,"voltage_std_V":float(np.std(voltages)),"normalized_voltage":normalized,"phase_lock_ratio":phase_lock,"nearest_step":nearest,"locking_error":max(abs(normalized-nearest),abs(phase_lock-nearest)),"effective_Ic_A":junction.effective_critical_current_A,"phase_offset_rad":junction.phase_offset_rad})
                if np.isclose(bias,float(c["output"].get("gif_bias_current_A",0))) and np.isclose(amp,float(c["output"].get("gif_microwave_amplitude_A",-1))): gif_rows=rows
            tag=f"{scenario['name']}_mw_{amp:.3g}"; f=(out/f"{tag}_sweep.csv").open("w",newline="",encoding="utf-8"); w=csv.DictWriter(f,fieldnames=points[0]); w.writeheader(); w.writerows(points); f.close(); axes[si,ai].plot(np.array(biases)*1e6,[p["normalized_voltage"] for p in points],"o-"); axes[si,ai].set(title=f"{scenario['name']}\nAC={amp*1e6:g} µA",xlabel="DC bias (µA)",ylabel="V/(Φ0 f)"); axes[si,ai].grid(alpha=.25)
            for n in range(-5,6): axes[si,ai].axhline(n,color="gray",linewidth=.4,linestyle="--")
            frequency_summary = None
            if gif_rows:
                prefix = out/tag
                frequency_summary = save_dynamic_diagnostics(sc,scenario,gif_rows,prefix)
                if c["output"].get("gif_enabled",False):
                    save_gif(sc,scenario,gif_rows,out/f"{tag}.gif")
            if len(c["vortex_scenarios"])==1 and len(amplitudes)==1:
                current=np.array(biases)*1e6; norm=np.array([p["normalized_voltage"] for p in points]); error=np.array([p["locking_error"] for p in points]); voltage=np.array([p["mean_voltage_V"] for p in points]); derivative=np.gradient(voltage,np.array(biases)); detail,da=plt.subplots(2,2,figsize=(11,8),constrained_layout=True); da[0,0].plot(current,voltage*1e6,"o-"); da[0,0].set(xlabel="DC bias (µA)",ylabel="Mean V (µV)",title="Fine I–V")
                da[0,1].plot(current,norm,"o-"); da[0,1].set(xlabel="DC bias (µA)",ylabel="V/(Φ0 f)",title="Normalized Shapiro plateaus"); da[1,0].plot(current,derivative,"o-"); da[1,0].set(xlabel="DC bias (µA)",ylabel="dV/dI (Ω)",title="Differential resistance"); da[1,1].semilogy(current,np.maximum(error,1e-12),"o-"); da[1,1].set(xlabel="DC bias (µA)",ylabel="Locking error",title="Voltage and phase-rate lock error")
                for a in da.flat: a.grid(alpha=.25)
                detail.savefig(out/"focused_shapiro_diagnostics.png",dpi=170); plt.close(detail)
            ss["microwave_sweeps"].append({"amplitude_A":amp,"effective_Ic_A":junction.effective_critical_current_A,"phase_offset_rad":junction.phase_offset_rad,"selected_bias_frequency_diagnostics":frequency_summary,"points":points})
        summaries.append(ss)
    fig.savefig(out/"shapiro_sweep_matrix.png",dpi=160); plt.close(fig); (out/"summary.json").write_text(json.dumps({"configuration":c,"voltage_step_V":quantum,"scenarios":summaries},indent=2),encoding="utf-8"); print(f"Vortex/microwave sweep written to {out}")

if __name__=="__main__": main()
