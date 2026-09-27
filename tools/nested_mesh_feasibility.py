"""Evaluate an isolated coarse-exterior/fine-interior representation."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shs.config.builder import build_simulation
from shs.experimental.adaptive_mesh import (
    composite_coarse_exterior,
    gauge_covariant_reconstruction,
)
from shs.solvers.coupled_solver import build_tdgl_model
from shs.solvers.tdgl_solver import tdgl_scales


def _rms(value, mask=None):
    array = np.asarray(value)
    if mask is not None:
        array = array[mask]
    return float(np.sqrt(np.mean(np.abs(array) ** 2)))


def _link_phases(psi, ax, ay, dx_dimensionless, dy_dimensionless):
    ux = np.exp(-1j * ax[:, :-1] * dx_dimensionless)
    uy = np.exp(-1j * ay[:-1, :] * dy_dimensionless)
    px = np.angle(np.conjugate(psi[:, :-1]) * ux * psi[:, 1:])
    py = np.angle(np.conjugate(psi[:-1, :]) * uy * psi[1:, :])
    return px, py


def _angular_error(first, second):
    return np.angle(np.exp(1j * (first - second)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path("tools/config/nested_mesh_feasibility.json"))
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    run = Path(config["source_run_directory"])
    state = json.loads((run / "state.json").read_text(encoding="utf-8"))
    checkpoint = run / state["checkpoint_file"]
    simulation = build_simulation(config["simulation_config"])
    model = build_tdgl_model(simulation)
    scales = tdgl_scales(simulation, model)
    with np.load(checkpoint) as archive:
        temperature = archive["temperature"]
        phonons = archive["phonon_temperature"]
        psi = archive["psi"]
        ax = archive["vector_potential_x"]
        ay = archive["vector_potential_y"]
    if temperature.shape != simulation.fields.temperature.shape:
        raise ValueError("Checkpoint and configured mesh shapes differ.")
    center_m = np.asarray(state["last_center_m"])
    center_cell = center_m / [simulation.mesh.dx, simulation.mesh.dy]
    radius_cells = float(config["fine_patch_radius_m"]) / min(simulation.mesh.dx,
                                                               simulation.mesh.dy)
    transition_cells = float(config["transition_width_m"]) / min(simulation.mesh.dx,
                                                                  simulation.mesh.dy)
    factor = int(config["coarse_factor"])
    reconstructed = {}
    fine_weight = None
    for name, field in (("temperature", temperature), ("phonon_temperature", phonons),
                        ("psi", psi), ("vector_potential_x", ax),
                        ("vector_potential_y", ay)):
        reconstructed[name], weight, _ = composite_coarse_exterior(
            field, factor=factor, fine_center_cell=center_cell,
            fine_radius_cells=radius_cells, transition_cells=transition_cells)
        fine_weight = weight if fine_weight is None else fine_weight
    exterior = fine_weight < 1e-12
    transition = (fine_weight > 0) & (fine_weight < 1)
    axd = scales.vector_potential_to_dimensionless(ax)
    ayd = scales.vector_potential_to_dimensionless(ay)
    axr = scales.vector_potential_to_dimensionless(reconstructed["vector_potential_x"])
    ayr = scales.vector_potential_to_dimensionless(reconstructed["vector_potential_y"])
    covariant_coarse = gauge_covariant_reconstruction(
        psi, axd, ayd, factor=factor,
        dx_dimensionless=simulation.mesh.dx/scales.xi,
        dy_dimensionless=simulation.mesh.dy/scales.xi)
    covariant_psi = fine_weight * psi + (1-fine_weight) * covariant_coarse
    links = _link_phases(psi, axd, ayd, simulation.mesh.dx/scales.xi,
                         simulation.mesh.dy/scales.xi)
    reconstructed_links = _link_phases(
        reconstructed["psi"], axr, ayr, simulation.mesh.dx/scales.xi,
        simulation.mesh.dy/scales.xi)
    covariant_links = _link_phases(
        covariant_psi, axd, ayd, simulation.mesh.dx/scales.xi,
        simulation.mesh.dy/scales.xi)
    phase_error_x = _angular_error(reconstructed_links[0], links[0])
    phase_error_y = _angular_error(reconstructed_links[1], links[1])
    covariant_phase_error_x = _angular_error(covariant_links[0], links[0])
    covariant_phase_error_y = _angular_error(covariant_links[1], links[1])
    coarse_nodes = ((simulation.mesh.nx - 1) // factor + 1) * (
        (simulation.mesh.ny - 1) // factor + 1)
    fine_nodes = int(np.count_nonzero(fine_weight > 0))
    overlap = int(np.count_nonzero((fine_weight[::factor, ::factor] > 0)))
    estimated_nodes = coarse_nodes + fine_nodes - overlap
    amplitude = np.abs(psi)
    amplitude_error = np.abs(np.abs(reconstructed["psi"]) - amplitude)
    temperature_error = reconstructed["temperature"] - temperature
    summary = {
        "source_run_directory": str(run), "checkpoint_step": state["step"],
        "coarse_factor": factor, "fine_patch_center_m": center_m.tolist(),
        "fine_patch_radius_m": config["fine_patch_radius_m"],
        "transition_width_m": config["transition_width_m"],
        "uniform_node_count": int(temperature.size),
        "estimated_composite_node_count": estimated_nodes,
        "estimated_node_fraction": estimated_nodes / temperature.size,
        "estimated_node_reduction_percent": 100 * (1 - estimated_nodes / temperature.size),
        "temperature_exterior_rms_error_K": _rms(temperature_error, exterior),
        "temperature_max_abs_error_K": float(np.max(np.abs(temperature_error))),
        "order_amplitude_exterior_rms_error": _rms(amplitude_error, exterior),
        "order_amplitude_max_abs_error": float(np.max(amplitude_error)),
        "covariant_x_link_phase_rms_error_rad": _rms(phase_error_x),
        "covariant_y_link_phase_rms_error_rad": _rms(phase_error_y),
        "parallel_transport_x_link_phase_rms_error_rad": _rms(covariant_phase_error_x),
        "parallel_transport_y_link_phase_rms_error_rad": _rms(covariant_phase_error_y),
        "parallel_transport_order_amplitude_exterior_rms_error": _rms(
            np.abs(covariant_psi)-amplitude, exterior),
        "transition_temperature_rms_error_K": _rms(temperature_error, transition)
        if np.any(transition) else 0.0,
        "interpretation": "This reconstructs a coarse exterior back onto the fine grid; it tests representation error but does not yet evolve a true nested mesh or demonstrate runtime speedup."
    }
    output = run / config["output_subdirectory"]
    output.mkdir(exist_ok=True)
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    plots = ((temperature, "inferno", "Reference electron T"),
             (reconstructed["temperature"], "inferno", "Composite electron T"),
             (temperature_error, "coolwarm", "T error (K)"),
             (amplitude, "viridis", "Reference |psi|"),
             (np.abs(reconstructed["psi"]), "viridis", "Composite |psi|"),
             (amplitude_error, "magma", "|psi| absolute error"))
    for axis, (field, cmap, title) in zip(axes.flat, plots):
        image = axis.imshow(field, origin="lower", cmap=cmap)
        axis.contour(fine_weight, levels=[0.01, .99], colors=["white", "cyan"],
                     linewidths=.7)
        axis.set_title(title)
        fig.colorbar(image, ax=axis)
    fig.savefig(output / "coarse_exterior_comparison.png", dpi=int(config["output_dpi"]))
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    comparison = ((amplitude, "Reference |psi|"),
                  (np.abs(covariant_psi), "Parallel-transport |psi|"),
                  (np.abs(covariant_psi)-amplitude, "Parallel-transport amplitude error"))
    for axis, (field, title) in zip(axes, comparison):
        image = axis.imshow(field, origin="lower",
                            cmap="coolwarm" if "error" in title else "viridis")
        axis.contour(fine_weight, levels=[0.01, .99], colors=["white", "cyan"], linewidths=.7)
        axis.set_title(title); fig.colorbar(image, ax=axis)
    fig.savefig(output / "gauge_covariant_transfer.png", dpi=int(config["output_dpi"]))
    plt.close(fig)
    print(output)
    return output


if __name__ == "__main__":
    main()
