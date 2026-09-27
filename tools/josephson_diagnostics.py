"""Run TDGL weak-link sweeps and create Josephson researcher plots."""

import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from tdgl_diagnostics import main as run_tdgl_diagnostics


def plot_junction_results(output):
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    records = []
    for case in summary["cases"]:
        path = output / case["label"] / "diagnostics.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        final = rows[-1]
        records.append((float(final["transport_current_A"]),
                        float(final["junction_voltage_V"]),
                        float(final["junction_phase_difference_rad"]),
                        float(final["junction_mean_amplitude"]), rows, case["label"]))
    records.sort(key=lambda item: item[0])
    current = np.array([r[0] for r in records])
    voltage = np.array([r[1] for r in records])
    phase = np.array([r[2] for r in records])
    amplitude = np.array([r[3] for r in records])
    resistance = np.gradient(voltage, current) if len(current) > 1 and np.all(np.diff(current)) else np.full_like(current, np.nan)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes[0, 0].plot(current * 1e9, voltage * 1e6, "o-")
    axes[0, 0].set(xlabel="Current (nA)", ylabel="Junction voltage (µV)", title="I–V characteristic")
    axes[0, 1].plot(current * 1e9, resistance, "o-")
    axes[0, 1].set(xlabel="Current (nA)", ylabel="dV/dI (Ω)", title="Differential resistance")
    axes[1, 0].plot(current * 1e9, phase, "o-")
    axes[1, 0].set(xlabel="Current (nA)", ylabel="Gauge-invariant phase drop (rad)", title="Junction phase response")
    axes[1, 1].plot(current * 1e9, amplitude, "o-")
    axes[1, 1].set(xlabel="Current (nA)", ylabel="Mean |ψ| in weak link", title="Weak-link suppression")
    for axis in axes.flat: axis.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(output / "josephson_sweep.png", dpi=160); plt.close(fig)
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for _, _, _, _, rows, label in records:
        time = np.array([float(r["time_s"]) for r in rows]) * 1e12
        axes[0].plot(time, [float(r["junction_voltage_V"]) * 1e6 for r in rows], label=label)
        axes[1].plot(time, [float(r["junction_phase_difference_rad"]) for r in rows], label=label)
        axes[2].plot(time, [float(r["cross_section_supercurrent_A"]) * 1e9 for r in rows], label=f"super: {label}")
        axes[2].plot(time, [float(r["cross_section_normal_current_A"]) * 1e9 for r in rows], "--", label=f"normal: {label}")
    axes[0].set(ylabel="Voltage (µV)", title="Junction transient response")
    axes[1].set(ylabel="Phase drop (rad)")
    axes[2].set(xlabel="Time (ps)", ylabel="Section current (nA)")
    axes[0].legend(fontsize=7)
    axes[2].legend(fontsize=6, ncol=2)
    for axis in axes: axis.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(output / "josephson_transients.png", dpi=160); plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("tools/config/josephson_current_sweep.json"))
    args = parser.parse_args()
    output = run_tdgl_diagnostics(["--config", str(args.config)])
    plot_junction_results(output)
    print(f"Josephson plots written to {output}")


if __name__ == "__main__": main()
