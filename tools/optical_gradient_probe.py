"""Diagnose stationary optical thermal gradients and the response of one vortex."""
import argparse
import csv
import json
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from shs.config.builder import build_simulation
from shs.config.simulation import LaserWaypointConfig
from shs.physics.pinning import pinning_suppression
from shs.solvers.coupled_solver import build_tdgl_model, coupled_step
from shs.solvers.tdgl_solver import tdgl_scales
from shs.utils.output import reserve_output_directory
from tdgl_diagnostics import apply_uniform_magnetic_field, detect_vortices, _seed_initial_state


def subcell_core(amplitude, x_cell, y_cell, radius=3):
    """Quadratic local |psi|^2 minimum; NaNs mean the fit cannot locate a core."""
    ny, nx = amplitude.shape
    x0, y0 = int(np.floor(x_cell)), int(np.floor(y_cell))
    xs = np.arange(max(0, x0-radius), min(nx, x0+radius+2))
    ys = np.arange(max(0, y0-radius), min(ny, y0+radius+2))
    xx, yy = np.meshgrid(xs-x_cell, ys-y_cell)
    if xx.size < 9:
        return float("nan"), float("nan")
    design = np.column_stack((np.ones(xx.size), xx.ravel(), yy.ravel(),
                              xx.ravel()**2, xx.ravel()*yy.ravel(), yy.ravel()**2))
    coefficients = np.linalg.lstsq(design, amplitude[np.ix_(ys, xs)].ravel()**2,
                                   rcond=None)[0]
    hessian = np.array([[2*coefficients[3], coefficients[4]],
                        [coefficients[4], 2*coefficients[5]]])
    if np.min(np.linalg.eigvalsh(hessian)) <= 0:
        return float("nan"), float("nan")
    center = -np.linalg.solve(hessian, coefficients[1:3])
    if np.max(np.abs(center)) > radius+1:
        return float("nan"), float("nan")
    return float(x_cell+center[0]), float(y_cell+center[1])


def complex_zero_core(psi, x_cell, y_cell):
    """Locate the bilinear complex-order-parameter zero in a vortex plaquette."""
    x0, y0 = int(np.floor(x_cell)), int(np.floor(y_cell))
    a = psi[y0, x0]
    b = psi[y0, x0+1]-a
    c = psi[y0+1, x0]-a
    d = psi[y0+1, x0+1]-a-b-c
    u = v = .5
    for _ in range(12):
        value = a+b*u+c*v+d*u*v
        jacobian = np.array([[np.real(b+d*v), np.real(c+d*u)],
                             [np.imag(b+d*v), np.imag(c+d*u)]])
        if abs(np.linalg.det(jacobian)) < 1e-14:
            return float("nan"), float("nan")
        delta = np.linalg.solve(jacobian, [np.real(value), np.imag(value)])
        u -= delta[0]
        v -= delta[1]
        if np.max(np.abs(delta)) < 1e-10:
            break
    if not (-1e-7 <= u <= 1+1e-7 and -1e-7 <= v <= 1+1e-7):
        return float("nan"), float("nan")
    return float(x0+u), float(y0+v)


def _sample(array, x_cell, y_cell):
    x0, y0 = int(np.floor(x_cell)), int(np.floor(y_cell))
    x0 = np.clip(x0, 0, array.shape[1]-2)
    y0 = np.clip(y0, 0, array.shape[0]-2)
    fx, fy = x_cell-x0, y_cell-y0
    return float((1-fx)*(1-fy)*array[y0, x0] + fx*(1-fy)*array[y0, x0+1]
                 + (1-fx)*fy*array[y0+1, x0] + fx*fy*array[y0+1, x0+1])


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def _snapshot(sim, model, scales, seed_m, time_s, step, phase, mode):
    f, mesh = sim.fields, sim.mesh
    vortices = detect_vortices(
        f.psi, scales.vector_potential_to_dimensionless(f.vector_potential_x),
        scales.vector_potential_to_dimensionless(f.vector_potential_y),
        mesh.dx/scales.xi, mesh.dy/scales.xi)
    winding = np.argwhere(vortices.winding == 1)
    seed_cell = np.asarray([seed_m[0]/mesh.dx, seed_m[1]/mesh.dy])
    if len(winding):
        y, x = min(winding, key=lambda p: np.hypot(p[1]+.5-seed_cell[0],
                                                   p[0]+.5-seed_cell[1]))
        plaquette_x, plaquette_y = x+.5, y+.5
        core_x, core_y = subcell_core(np.abs(f.psi), plaquette_x, plaquette_y)
        zero_x, zero_y = complex_zero_core(f.psi, plaquette_x, plaquette_y)
    else:
        plaquette_x = plaquette_y = core_x = core_y = zero_x = zero_y = float("nan")
    sample_x = zero_x if np.isfinite(zero_x) else seed_cell[0]
    sample_y = zero_y if np.isfinite(zero_y) else seed_cell[1]
    grad_ty, grad_tx = np.gradient(f.temperature, mesh.dy, mesh.dx)
    epsilon = model.temperature_coefficient(f.temperature/sim.material_map.Tc)
    epsilon -= pinning_suppression(sim)
    grad_ey, grad_ex = np.gradient(epsilon, mesh.dy, mesh.dx)
    phonons = f.phonon_temperature if f.phonon_temperature is not None else f.temperature
    peak_y, peak_x = np.unravel_index(np.argmax(f.temperature), f.temperature.shape)
    row = {
        "mode": mode, "phase": phase, "step": step, "time_s": time_s,
        "positive_vortices": vortices.positive_count,
        "negative_vortices": vortices.negative_count,
        "vortex_x_plaquette_m": plaquette_x*mesh.dx,
        "vortex_y_plaquette_m": plaquette_y*mesh.dy,
        "vortex_x_subcell_m": core_x*mesh.dx,
        "vortex_y_subcell_m": core_y*mesh.dy,
        "subcell_fit_valid": bool(np.isfinite(core_x)),
        "vortex_x_zero_m": zero_x*mesh.dx,
        "vortex_y_zero_m": zero_y*mesh.dy,
        "complex_zero_valid": bool(np.isfinite(zero_x)),
        "electron_temperature_at_core_K": _sample(f.temperature, sample_x, sample_y),
        "phonon_temperature_at_core_K": _sample(phonons, sample_x, sample_y),
        "electron_gradient_x_at_core_K_per_m": _sample(grad_tx, sample_x, sample_y),
        "electron_gradient_y_at_core_K_per_m": _sample(grad_ty, sample_x, sample_y),
        "gl_coefficient_gradient_x_at_core_per_m": _sample(grad_ex, sample_x, sample_y),
        "gl_coefficient_gradient_y_at_core_per_m": _sample(grad_ey, sample_x, sample_y),
        "maximum_electron_temperature_K": float(np.max(f.temperature)),
        "electron_temperature_peak_x_m": float(peak_x*mesh.dx),
        "electron_temperature_peak_y_m": float(peak_y*mesh.dy),
        "maximum_phonon_temperature_K": float(np.max(phonons)),
        "mean_order_amplitude": float(np.mean(np.abs(f.psi))),
        "laser_x_m": f.laser_position_x_m,
        "laser_y_m": f.laser_position_y_m,
    }
    return row


def _plot_results(rows, final_fields, output, config, seed_m, mesh):
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True,
                             constrained_layout=True)
    for mode in ("dark", "plus", "minus"):
        subset = [row for row in rows if row["mode"] == mode]
        time = np.array([row["time_s"] for row in subset])*1e12
        axes[0].plot(time, [row["electron_gradient_y_at_core_K_per_m"]*1e-9
                            for row in subset], label=mode)
        axes[1].plot(time, [(row["vortex_y_zero_m"]-seed_m[1])*1e9
                            for row in subset], label=mode)
        axes[2].plot(time, [row["electron_temperature_at_core_K"]
                            for row in subset], label=mode)
    axes[0].set_ylabel("Electron dT/dy at core (K/nm)")
    axes[1].set_ylabel("Complex-zero core y - seed (nm)")
    axes[2].set_ylabel("Electron T at core (K)")
    axes[2].set_xlabel("Time after dark preparation (ps)")
    for axis in axes:
        axis.grid(alpha=.3)
        axis.legend()
    fig.savefig(output/"gradient_and_response.png", dpi=config["output"]["dpi"])
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True)
    for column, mode in enumerate(("dark", "plus", "minus")):
        fields = final_fields[mode]
        for axis, array, label in ((axes[0, column], fields.temperature, "Electron T (K)"),
                                   (axes[1, column], np.abs(fields.psi), "|psi|")):
            image = axis.imshow(array, origin="lower")
            axis.set_title(f"{mode}: {label}")
            axis.plot(seed_m[0]/mesh.dx, seed_m[1]/mesh.dy,
                      "wx", markersize=6)
            fig.colorbar(image, ax=axis)
    fig.savefig(output/"final_film_maps.png", dpi=config["output"]["dpi"])
    plt.close(fig)


def response_summary(final, critical_temperature, threshold_nm):
    """Quantify beam-direction response after subtracting common dark drift."""
    dark = final["dark"]
    plus = final["plus"]
    minus = final["minus"]
    zero_plus = (plus["vortex_y_zero_m"]-dark["vortex_y_zero_m"])*1e9
    zero_minus = (minus["vortex_y_zero_m"]-dark["vortex_y_zero_m"])*1e9
    fit_plus = (plus["vortex_y_subcell_m"]-dark["vortex_y_subcell_m"])*1e9
    fit_minus = (minus["vortex_y_subcell_m"]-dark["vortex_y_subcell_m"])*1e9
    directional = bool(
        zero_plus >= threshold_nm and zero_minus <= -threshold_nm
        and fit_plus > 0 and fit_minus < 0
        and plus["electron_gradient_y_at_core_K_per_m"] > 0
        and minus["electron_gradient_y_at_core_K_per_m"] < 0
        and max(plus["maximum_electron_temperature_K"],
                minus["maximum_electron_temperature_K"]) < critical_temperature
        and all(final[mode]["positive_vortices"] == 1 and
                final[mode]["negative_vortices"] == 0 for mode in final)
    )
    return {"complex_zero_plus_minus_dark_nm": zero_plus,
            "complex_zero_minus_minus_dark_nm": zero_minus,
            "amplitude_fit_plus_minus_dark_nm": fit_plus,
            "amplitude_fit_minus_minus_dark_nm": fit_minus,
            "directional_response_in_model": directional,
            "critical_temperature_K": critical_temperature,
            "minimum_response_nm": threshold_nm}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path("tools/config/optical_gradient_probe.json"))
    parser.add_argument("--dark-steps", type=int, help="Short smoke-test override.")
    parser.add_argument("--exposure-steps", type=int, help="Short smoke-test override.")
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    output = reserve_output_directory(config["output_directory"])
    sim = build_simulation(config["simulation_config"])
    sim.config.current = 0.0
    sim.config.electrical.drive_mode = "current"
    sim.config.laser.enabled = False
    sim.config.thermal.model = "two_temperature"
    sim.config.thermal.validate()
    apply_uniform_magnetic_field(sim, config["applied_field_T"])
    _seed_initial_state(sim, {"vortices": [config["initial_vortex"]]})
    seed_m = np.asarray(config["seed_position_m"])
    dt = sim.config.dt
    model = build_tdgl_model(sim)
    scales = tdgl_scales(sim, model)
    dark_steps = config["timing"]["dark_steps"] if args.dark_steps is None else args.dark_steps
    exposure_steps = (config["timing"]["exposure_steps"] if args.exposure_steps is None
                      else args.exposure_steps)
    if dark_steps < 0 or exposure_steps < 1:
        raise ValueError("Dark steps must be nonnegative; exposure steps must be positive.")
    dark_change = []
    for step in range(dark_steps):
        before = sim.fields.psi.copy()
        _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=step*dt)
        if not result.converged:
            raise RuntimeError(f"Dark preparation failed at step {step+1}.")
        dark_change.append(float(np.max(np.abs(sim.fields.psi-before))))
    prepared = deepcopy(sim.fields)
    rows, final_fields = [], {}
    sample_every = int(config["timing"]["sample_every"])
    for mode, offset in (("dark", None), ("plus", 1), ("minus", -1)):
        sim.fields = deepcopy(prepared)
        sim.config.laser.enabled = offset is not None
        if offset is not None:
            target = seed_m + np.asarray(config["beam_offset_m"])*offset
            sim.config.laser.absorbed_power_W = config["absorbed_power_W"]
            sim.config.laser.sigma_m = config["spot_sigma_m"]
            sim.config.laser.waypoints = [LaserWaypointConfig(0.0, *target),
                                          LaserWaypointConfig(exposure_steps*dt, *target)]
            sim.config.laser.validate()
        for step in range(exposure_steps+1):
            if step:
                _, result = coupled_step(sim, dt=dt, tdgl_model=model,
                                         time_s=(step-1)*dt)
                if not result.converged:
                    raise RuntimeError(f"{mode} failed at step {step}.")
            if step % sample_every == 0 or step == exposure_steps:
                rows.append(_snapshot(sim, model, scales, seed_m, step*dt,
                                      step, "exposure", mode))
            if step and step % 50 == 0:
                print(f"{mode}: {step}/{exposure_steps}", flush=True)
        final_fields[mode] = deepcopy(sim.fields)
    _write_csv(output/"force_probe.csv", rows)
    _plot_results(rows, final_fields, output, config, seed_m, sim.mesh)
    summary = {
        "config": config, "dark_preparation_last_psi_change":
        dark_change[-1] if dark_change else None,
        "dark_preparation_max_psi_change": max(dark_change) if dark_change else None,
        "final": {mode: [row for row in rows if row["mode"] == mode][-1]
                  for mode in final_fields},
        "interpretation": "Temperature and GL-coefficient gradients are local drive proxies, not a calibrated force decomposition. The bilinear complex zero and quadratic |psi| minimum are independent subcell position estimates; compare both with plaquette winding and dark control.",
    }
    summary["response"] = response_summary(
        summary["final"], float(np.min(sim.material_map.Tc)),
        float(config["minimum_directional_response_nm"]))
    (output/"summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(output)
    return output


if __name__ == "__main__":
    main()
