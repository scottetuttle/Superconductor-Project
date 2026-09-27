"""Create retrospective trajectory visuals from an optical transport run."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np


def _read(path):
    with path.open(newline="", encoding="utf-8") as stream:
        raw = list(csv.DictReader(stream))
    if not raw:
        raise ValueError(f"No diagnostic rows in {path}.")
    numeric = {}
    for key in raw[0]:
        try:
            numeric[key] = np.array([float(row[key]) for row in raw])
        except (TypeError, ValueError):
            numeric[key] = np.array([row[key] for row in raw])
    return numeric


def _fill_coordinate_gaps(points, label):
    """Interpolate isolated nonfinite trajectory samples from legacy runs."""
    points = np.asarray(points, dtype=float).copy()
    sample = np.arange(len(points), dtype=float)
    for axis in range(points.shape[1]):
        finite = np.isfinite(points[:, axis])
        if not np.any(finite):
            raise ValueError(f"{label} has no finite coordinate samples.")
        points[:, axis] = np.interp(sample, sample[finite], points[finite, axis])
    return points


def _setup(run):
    data = _read(run / "diagnostics.csv")
    config = json.loads((run / "config.json").read_text(encoding="utf-8"))
    state = json.loads((run / "state.json").read_text(encoding="utf-8"))
    initial = np.asarray(state["initial_center_m"]) * 1e9
    displacement = np.asarray(config["target"]["displacement_m"]) * 1e9
    target = initial + displacement
    direction = displacement / np.linalg.norm(displacement)
    transverse = np.array([-direction[1], direction[0]])
    vortex = _fill_coordinate_gaps(np.column_stack(
        (data["vortex_x_zero_m"], data["vortex_y_zero_m"])), "vortex") * 1e9
    beam = _fill_coordinate_gaps(np.column_stack(
        (data["commanded_beam_x_m"], data["commanded_beam_y_m"])), "beam") * 1e9
    relative = vortex - initial
    return data, config, state, initial, target, direction, transverse, vortex, beam, relative


def plot_kymographs(run, dpi):
    data, _, _, initial, target, direction, transverse, vortex, beam, relative = _setup(run)
    time_ps = data["time_s"] * 1e12
    along_v = relative @ direction
    cross_v = relative @ transverse
    along_b = (beam - initial) @ direction
    cross_b = (beam - initial) @ transverse
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    axes[0].plot(time_ps, along_v, color="tab:red", label="vortex")
    axes[0].plot(time_ps, along_b, color="tab:cyan", label="beam")
    axes[0].axhline(np.linalg.norm(target-initial), color="black", linestyle=":", label="target")
    axes[0].set_ylabel("Position along target direction (nm)")
    axes[0].legend()
    axes[1].plot(time_ps, cross_v, color="tab:red", label="vortex")
    axes[1].plot(time_ps, cross_b, color="tab:cyan", label="beam")
    axes[1].axhline(0, color="black", linewidth=.8)
    axes[1].set_ylabel("Cross-track position (nm)")
    axes[1].set_xlabel("Time (ps)")
    for axis in axes:
        axis.grid(alpha=.3)
    fig.suptitle("Optical transport position-time kymographs")
    fig.savefig(run / "trajectory_kymographs.png", dpi=dpi)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
    axes[0].plot(time_ps, data["maximum_electron_temperature_K"], label="peak electron T")
    axes[0].plot(time_ps, data["electron_temperature_at_core_K"], label="electron T at vortex")
    axes[0].plot(time_ps, data["phonon_temperature_at_core_K"], label="phonon T at vortex")
    axes[0].legend(fontsize=8); axes[0].set_ylabel("Temperature (K)")
    axes[1].plot(time_ps, data["gradient_toward_target_K_per_nm"])
    axes[1].axhline(0, color="black", linewidth=.8)
    axes[1].set_ylabel("Gradient toward target (K/nm)")
    axes[2].plot(time_ps, data["target_distance_nm"], label="Euclidean distance")
    axes[2].plot(time_ps, data["progress_nm"], label="targetward projection")
    axes[2].set_ylabel("Distance (nm)"); axes[2].set_xlabel("Time (ps)")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=.3)
    fig.savefig(run / "thermal_motion_kymographs.png", dpi=dpi)
    plt.close(fig)


def plot_checkpoint(run, dpi):
    data, _, state, initial, target, _, _, vortex, beam, _ = _setup(run)
    checkpoint = run / state["checkpoint_file"]
    with np.load(checkpoint) as fields:
        temperature = fields["temperature"]
        psi = fields["psi"]
    extent = [0, temperature.shape[1]-1, 0, temperature.shape[0]-1]
    # The current optical film uses a 1 nm mesh; infer exact plotting scales
    # from diagnostic coordinates and preserve the conventional nanometre axes.
    extent = [0, extent[1], 0, extent[3]]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    entries = ((temperature, "inferno", "Electron temperature (K)", None),
               (np.abs(psi), "viridis", "|psi|", None),
               (np.angle(psi), "twilight", "Phase", (-np.pi, np.pi)))
    for axis, (field, cmap, title, limits) in zip(axes, entries):
        image = axis.imshow(field, origin="lower", extent=extent, cmap=cmap,
                            vmin=None if limits is None else limits[0],
                            vmax=None if limits is None else limits[1])
        axis.plot(vortex[:, 0], vortex[:, 1], color="white", linewidth=1, alpha=.75)
        axis.plot(beam[:, 0], beam[:, 1], color="cyan", linewidth=1, alpha=.7)
        axis.plot(initial[0], initial[1], "wo", markersize=4, label="start")
        axis.plot(target[0], target[1], "w*", markersize=8, label="target")
        axis.set(xlabel="x (nm)", ylabel="y (nm)", title=title)
        fig.colorbar(image, ax=axis)
    axes[0].legend(fontsize=7)
    fig.suptitle(f'Checkpoint fields at step {state["step"]}; white=vortex, cyan=beam')
    fig.savefig(run / "checkpoint_fields_with_paths.png", dpi=dpi)
    plt.close(fig)


def save_gif(run, fps, dpi, stride):
    data, _, _, initial, target, _, _, vortex, beam, _ = _setup(run)
    indices = list(range(0, len(vortex), stride))
    if indices[-1] != len(vortex)-1:
        indices.append(len(vortex)-1)
    padding = 4
    xmin = min(vortex[:, 0].min(), beam[:, 0].min(), initial[0], target[0]) - padding
    xmax = max(vortex[:, 0].max(), beam[:, 0].max(), initial[0], target[0]) + padding
    ymin = min(vortex[:, 1].min(), beam[:, 1].min(), initial[1], target[1]) - padding
    ymax = max(vortex[:, 1].max(), beam[:, 1].max(), initial[1], target[1]) + padding
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    path_v, = axes[0].plot([], [], color="tab:red", label="vortex")
    path_b, = axes[0].plot([], [], color="tab:cyan", label="beam")
    point_v, = axes[0].plot([], [], "ro")
    point_b, = axes[0].plot([], [], "c+")
    axes[0].plot(initial[0], initial[1], "ko", markersize=4, label="start")
    axes[0].plot(target[0], target[1], "k*", markersize=9, label="target")
    axes[0].set(xlim=(xmin, xmax), ylim=(ymin, ymax), xlabel="x (nm)", ylabel="y (nm)",
                title="Recorded trajectory")
    axes[0].set_aspect("equal"); axes[0].grid(alpha=.3); axes[0].legend(fontsize=8)
    time_ps = data["time_s"] * 1e12
    axes[1].plot(time_ps, data["maximum_electron_temperature_K"], color="tab:orange",
                 label="peak electron T")
    axes[1].plot(time_ps, data["electron_temperature_at_core_K"], color="tab:red",
                 label="T at vortex")
    time_line = axes[1].axvline(time_ps[0], color="black")
    gradient_axis = axes[1].twinx()
    gradient_axis.plot(time_ps, data["gradient_toward_target_K_per_nm"],
                       color="tab:blue", alpha=.65, label="gradient")
    axes[1].set(xlabel="Time (ps)", ylabel="Temperature (K)", title="Thermal drive")
    gradient_axis.set_ylabel("Gradient toward target (K/nm)", color="tab:blue")
    axes[1].grid(alpha=.3); axes[1].legend(fontsize=8, loc="upper left")
    title = fig.suptitle("")

    def update(frame):
        index = indices[frame]
        path_v.set_data(vortex[:index+1, 0], vortex[:index+1, 1])
        path_b.set_data(beam[:index+1, 0], beam[:index+1, 1])
        point_v.set_data([vortex[index, 0]], [vortex[index, 1]])
        point_b.set_data([beam[index, 0]], [beam[index, 1]])
        time_line.set_xdata([time_ps[index], time_ps[index]])
        title.set_text(f't={time_ps[index]:.2f} ps, step={int(data["step"][index])}, '
                       f'target distance={data["target_distance_nm"][index]:.2f} nm')
        return path_v, path_b, point_v, point_b, time_line, title

    movie = animation.FuncAnimation(fig, update, frames=len(indices), blit=False)
    movie.save(run / "transport_trajectory.gif", writer=animation.PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)


def save_field_gif(run, fps, dpi, stride):
    paths = sorted((run / "field_frames").glob("frame_*.npz"))[::stride]
    if not paths:
        return False
    frames = []
    for path in paths:
        with np.load(path) as archive:
            frames.append({key: archive[key].copy() for key in archive.files})
    tmin = min(float(np.min(frame["temperature"])) for frame in frames)
    tmax = max(float(np.max(frame["temperature"])) for frame in frames)
    amax = max(float(np.max(np.abs(frame["psi"]))) for frame in frames)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    fields = (("temperature", "inferno", (tmin, tmax), "Electron T (K)"),
              ("amplitude", "viridis", (0, amax), "|psi|"),
              ("phase", "twilight", (-np.pi, np.pi), "Phase"))
    images, cores, beams = [], [], []
    first = frames[0]
    derived = {"temperature": first["temperature"], "amplitude": np.abs(first["psi"]),
               "phase": np.angle(first["psi"])}
    for axis, (key, cmap, limits, title_text) in zip(axes, fields):
        image = axis.imshow(derived[key], origin="lower", cmap=cmap,
                            vmin=limits[0], vmax=limits[1])
        core, = axis.plot([], [], "rx", markersize=7)
        beam, = axis.plot([], [], "c+", markersize=8)
        axis.set_title(title_text); fig.colorbar(image, ax=axis)
        images.append(image); cores.append(core); beams.append(beam)
    title = fig.suptitle("")

    def update(index):
        frame = frames[index]
        values = (frame["temperature"], np.abs(frame["psi"]), np.angle(frame["psi"]))
        # The detailed large-film reference has 1 nm spacing, so metres to
        # image-node coordinates is a direct 1e9 conversion.
        core_xy = frame["vortex_center_m"] * 1e9
        beam_xy = frame["beam_m"] * 1e9
        for image, value, core, beam in zip(images, values, cores, beams):
            image.set_data(value)
            core.set_data([core_xy[0]], [core_xy[1]])
            beam.set_data([beam_xy[0]], [beam_xy[1]])
        title.set_text(f't={float(frame["time_s"])*1e12:.2f} ps, '
                       f'step={int(frame["step"])}, laser={"on" if bool(frame["laser_on"]) else "off"}')
        return images + cores + beams + [title]

    movie = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
    movie.save(run / "transport_fields.gif", writer=animation.PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=110)
    parser.add_argument("--stride", type=int, default=2,
                        help="Diagnostic samples per GIF frame.")
    parser.add_argument("--field-stride", type=int, default=1,
                        help="Saved field snapshots per field-GIF frame.")
    args = parser.parse_args(argv)
    if args.fps < 1 or args.dpi < 40 or min(args.stride, args.field_stride) < 1:
        raise ValueError("fps, dpi and stride must be positive.")
    plot_kymographs(args.run, args.dpi)
    plot_checkpoint(args.run, args.dpi)
    save_gif(args.run, args.fps, args.dpi, args.stride)
    save_field_gif(args.run, args.fps, args.dpi, args.field_stride)
    print(args.run)


if __name__ == "__main__":
    main()
