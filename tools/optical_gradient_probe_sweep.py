"""Sweep stationary optical hotspot offset, power, and width for one vortex.

Every illuminated branch starts from the same dark-prepared state.  The dark
control is evolved once, then each (+offset, -offset) pair is run independently
for every Cartesian combination of beam offset, absorbed power, and spot width.
"""

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
from optical_gradient_probe import _plot_results as _plot_case_results


def subcell_core(amplitude, x_cell, y_cell, radius=3):
    """Quadratic local |psi|^2 minimum; NaNs mean the fit cannot locate a core."""
    ny, nx = amplitude.shape
    x0, y0 = int(np.floor(x_cell)), int(np.floor(y_cell))
    xs = np.arange(max(0, x0 - radius), min(nx, x0 + radius + 2))
    ys = np.arange(max(0, y0 - radius), min(ny, y0 + radius + 2))
    xx, yy = np.meshgrid(xs - x_cell, ys - y_cell)
    if xx.size < 9:
        return float("nan"), float("nan")
    design = np.column_stack(
        (
            np.ones(xx.size),
            xx.ravel(),
            yy.ravel(),
            xx.ravel() ** 2,
            xx.ravel() * yy.ravel(),
            yy.ravel() ** 2,
        )
    )
    coefficients = np.linalg.lstsq(
        design, amplitude[np.ix_(ys, xs)].ravel() ** 2, rcond=None
    )[0]
    hessian = np.array(
        [
            [2 * coefficients[3], coefficients[4]],
            [coefficients[4], 2 * coefficients[5]],
        ]
    )
    if np.min(np.linalg.eigvalsh(hessian)) <= 0:
        return float("nan"), float("nan")
    center = -np.linalg.solve(hessian, coefficients[1:3])
    if np.max(np.abs(center)) > radius + 1:
        return float("nan"), float("nan")
    return float(x_cell + center[0]), float(y_cell + center[1])


def complex_zero_core(psi, x_cell, y_cell):
    """Locate the bilinear complex-order-parameter zero in a vortex plaquette."""
    x0, y0 = int(np.floor(x_cell)), int(np.floor(y_cell))
    a = psi[y0, x0]
    b = psi[y0, x0 + 1] - a
    c = psi[y0 + 1, x0] - a
    d = psi[y0 + 1, x0 + 1] - a - b - c
    u = v = 0.5
    for _ in range(12):
        value = a + b * u + c * v + d * u * v
        jacobian = np.array(
            [
                [np.real(b + d * v), np.real(c + d * u)],
                [np.imag(b + d * v), np.imag(c + d * u)],
            ]
        )
        if abs(np.linalg.det(jacobian)) < 1e-14:
            return float("nan"), float("nan")
        delta = np.linalg.solve(jacobian, [np.real(value), np.imag(value)])
        u -= delta[0]
        v -= delta[1]
        if np.max(np.abs(delta)) < 1e-10:
            break
    if not (-1e-7 <= u <= 1 + 1e-7 and -1e-7 <= v <= 1 + 1e-7):
        return float("nan"), float("nan")
    return float(x0 + u), float(y0 + v)


def _sample(array, x_cell, y_cell):
    x0, y0 = int(np.floor(x_cell)), int(np.floor(y_cell))
    x0 = np.clip(x0, 0, array.shape[1] - 2)
    y0 = np.clip(y0, 0, array.shape[0] - 2)
    fx, fy = x_cell - x0, y_cell - y0
    return float(
        (1 - fx) * (1 - fy) * array[y0, x0]
        + fx * (1 - fy) * array[y0, x0 + 1]
        + (1 - fx) * fy * array[y0 + 1, x0]
        + fx * fy * array[y0 + 1, x0 + 1]
    )


def _write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    temporary = path.with_suffix(path.suffix+".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _snapshot(
    sim,
    model,
    scales,
    seed_m,
    time_s,
    step,
    phase,
    mode,
    *,
    offset_vector_m=None,
    absorbed_power_W=None,
    spot_sigma_m=None,
):
    f, mesh = sim.fields, sim.mesh
    vortices = detect_vortices(
        f.psi,
        scales.vector_potential_to_dimensionless(f.vector_potential_x),
        scales.vector_potential_to_dimensionless(f.vector_potential_y),
        mesh.dx / scales.xi,
        mesh.dy / scales.xi,
    )
    winding = np.argwhere(vortices.winding == 1)
    seed_cell = np.asarray([seed_m[0] / mesh.dx, seed_m[1] / mesh.dy])
    if len(winding):
        y, x = min(
            winding,
            key=lambda p: np.hypot(
                p[1] + 0.5 - seed_cell[0], p[0] + 0.5 - seed_cell[1]
            ),
        )
        plaquette_x, plaquette_y = x + 0.5, y + 0.5
        core_x, core_y = subcell_core(np.abs(f.psi), plaquette_x, plaquette_y)
        zero_x, zero_y = complex_zero_core(f.psi, plaquette_x, plaquette_y)
    else:
        plaquette_x = plaquette_y = core_x = core_y = zero_x = zero_y = float("nan")

    sample_x = zero_x if np.isfinite(zero_x) else seed_cell[0]
    sample_y = zero_y if np.isfinite(zero_y) else seed_cell[1]
    grad_ty, grad_tx = np.gradient(f.temperature, mesh.dy, mesh.dx)
    epsilon = model.temperature_coefficient(f.temperature / sim.material_map.Tc)
    epsilon -= pinning_suppression(sim)
    grad_ey, grad_ex = np.gradient(epsilon, mesh.dy, mesh.dx)
    phonons = f.phonon_temperature if f.phonon_temperature is not None else f.temperature
    peak_y, peak_x = np.unravel_index(np.argmax(f.temperature), f.temperature.shape)

    if offset_vector_m is None:
        offset_x_m = offset_y_m = offset_distance_m = float("nan")
    else:
        offset_vector_m = np.asarray(offset_vector_m, dtype=float)
        offset_x_m = float(offset_vector_m[0])
        offset_y_m = float(offset_vector_m[1])
        offset_distance_m = float(np.linalg.norm(offset_vector_m))

    return {
        "mode": mode,
        "phase": phase,
        "step": step,
        "time_s": time_s,
        "offset_x_m": offset_x_m,
        "offset_y_m": offset_y_m,
        "offset_distance_m": offset_distance_m,
        "absorbed_power_W": absorbed_power_W,
        "spot_sigma_m": spot_sigma_m,
        "positive_vortices": vortices.positive_count,
        "negative_vortices": vortices.negative_count,
        "vortex_x_plaquette_m": plaquette_x * mesh.dx,
        "vortex_y_plaquette_m": plaquette_y * mesh.dy,
        "vortex_x_subcell_m": core_x * mesh.dx,
        "vortex_y_subcell_m": core_y * mesh.dy,
        "subcell_fit_valid": bool(np.isfinite(core_x)),
        "vortex_x_zero_m": zero_x * mesh.dx,
        "vortex_y_zero_m": zero_y * mesh.dy,
        "complex_zero_valid": bool(np.isfinite(zero_x)),
        "electron_temperature_at_core_K": _sample(f.temperature, sample_x, sample_y),
        "phonon_temperature_at_core_K": _sample(phonons, sample_x, sample_y),
        "electron_gradient_x_at_core_K_per_m": _sample(grad_tx, sample_x, sample_y),
        "electron_gradient_y_at_core_K_per_m": _sample(grad_ty, sample_x, sample_y),
        "gl_coefficient_gradient_x_at_core_per_m": _sample(grad_ex, sample_x, sample_y),
        "gl_coefficient_gradient_y_at_core_per_m": _sample(grad_ey, sample_x, sample_y),
        "maximum_electron_temperature_K": float(np.max(f.temperature)),
        "electron_temperature_peak_x_m": float(peak_x * mesh.dx),
        "electron_temperature_peak_y_m": float(peak_y * mesh.dy),
        "maximum_phonon_temperature_K": float(np.max(phonons)),
        "mean_order_amplitude": float(np.mean(np.abs(f.psi))),
        "laser_x_m": f.laser_position_x_m,
        "laser_y_m": f.laser_position_y_m,
    }


def _as_sweep(config):
    """Read sweep settings while remaining compatible with the old single-case config."""
    sweeps = config.get("sweeps", {})

    offsets = sweeps.get("beam_offsets_m")
    if offsets is None:
        offsets = [config["beam_offset_m"]]

    powers = sweeps.get("absorbed_powers_W")
    if powers is None:
        powers = [config["absorbed_power_W"]]

    sigmas = sweeps.get("spot_sigmas_m")
    if sigmas is None:
        sigmas = [config["spot_sigma_m"]]

    offsets = [np.asarray(value, dtype=float) for value in offsets]
    powers = [float(value) for value in powers]
    sigmas = [float(value) for value in sigmas]

    if not offsets or not powers or not sigmas:
        raise ValueError("Sweep offset, power, and spot-width lists must be nonempty.")

    for offset in offsets:
        if offset.shape != (2,):
            raise ValueError(
                "Each beam offset must be [dx_m, dy_m]. For a sweep use "
                '"sweeps": {"beam_offsets_m": [[0, 1e-8], [0, 2e-8], ...]}.'
            )
        if not np.all(np.isfinite(offset)) or np.linalg.norm(offset) <= 0:
            raise ValueError("Each beam offset must be finite and nonzero.")
    if any(power <= 0 or not np.isfinite(power) for power in powers):
        raise ValueError("All absorbed powers must be finite and positive.")
    if any(sigma <= 0 or not np.isfinite(sigma) for sigma in sigmas):
        raise ValueError("All spot sigmas must be finite and positive.")

    return offsets, powers, sigmas


def _core_position(row, prefix="zero"):
    if prefix == "zero":
        return np.array([row["vortex_x_zero_m"], row["vortex_y_zero_m"]], dtype=float)
    if prefix == "subcell":
        return np.array([row["vortex_x_subcell_m"], row["vortex_y_subcell_m"]], dtype=float)
    raise ValueError(prefix)


def response_summary(final, offset_vector_m, critical_temperature, threshold_nm):
    """Quantify response projected along the actual beam direction."""
    dark = final["dark"]
    plus = final["plus"]
    minus = final["minus"]

    offset_vector_m = np.asarray(offset_vector_m, dtype=float)
    offset_distance_m = float(np.linalg.norm(offset_vector_m))
    unit = offset_vector_m / offset_distance_m

    dark_zero = _core_position(dark, "zero")
    plus_zero = _core_position(plus, "zero")
    minus_zero = _core_position(minus, "zero")
    dark_fit = _core_position(dark, "subcell")
    plus_fit = _core_position(plus, "subcell")
    minus_fit = _core_position(minus, "subcell")

    # Project motion onto each branch's own beam direction. Positive means the
    # vortex moved toward that branch's laser location relative to dark drift.
    plus_zero_nm = float(np.dot(plus_zero - dark_zero, unit) * 1e9)
    minus_zero_nm = float(np.dot(minus_zero - dark_zero, -unit) * 1e9)
    plus_fit_nm = float(np.dot(plus_fit - dark_fit, unit) * 1e9)
    minus_fit_nm = float(np.dot(minus_fit - dark_fit, -unit) * 1e9)

    # Symmetric response is the mean of the two opposed controls. It suppresses
    # common drift and is the primary y-value used by the sweep plots.
    mean_directional_nm = 0.5 * (plus_zero_nm + minus_zero_nm)
    asymmetry_nm = 0.5 * (plus_zero_nm - minus_zero_nm)

    grad_plus = np.array(
        [
            plus["electron_gradient_x_at_core_K_per_m"],
            plus["electron_gradient_y_at_core_K_per_m"],
        ]
    )
    grad_minus = np.array(
        [
            minus["electron_gradient_x_at_core_K_per_m"],
            minus["electron_gradient_y_at_core_K_per_m"],
        ]
    )
    plus_gradient_toward_beam = float(np.dot(grad_plus, unit))
    minus_gradient_toward_beam = float(np.dot(grad_minus, -unit))

    valid_topology = all(
        final[mode]["positive_vortices"] == 1 and final[mode]["negative_vortices"] == 0
        for mode in final
    )
    below_tc = max(
        plus["maximum_electron_temperature_K"],
        minus["maximum_electron_temperature_K"],
    ) < critical_temperature

    directional = bool(
        plus_zero_nm >= threshold_nm
        and minus_zero_nm >= threshold_nm
        and plus_fit_nm > 0
        and minus_fit_nm > 0
        and plus_gradient_toward_beam > 0
        and minus_gradient_toward_beam > 0
        and below_tc
        and valid_topology
    )

    return {
        "offset_x_m": float(offset_vector_m[0]),
        "offset_y_m": float(offset_vector_m[1]),
        "offset_distance_m": offset_distance_m,
        "offset_distance_nm": offset_distance_m * 1e9,
        "plus_motion_toward_beam_nm": plus_zero_nm,
        "minus_motion_toward_beam_nm": minus_zero_nm,
        "mean_motion_toward_beam_nm": mean_directional_nm,
        "directional_asymmetry_nm": asymmetry_nm,
        "plus_amplitude_fit_toward_beam_nm": plus_fit_nm,
        "minus_amplitude_fit_toward_beam_nm": minus_fit_nm,
        "plus_temperature_gradient_toward_beam_K_per_m": plus_gradient_toward_beam,
        "minus_temperature_gradient_toward_beam_K_per_m": minus_gradient_toward_beam,
        "directional_response_in_model": directional,
        "critical_temperature_K": critical_temperature,
        "minimum_response_nm": threshold_nm,
        "topology_preserved": valid_topology,
        "maximum_electron_temperature_K": max(
            plus["maximum_electron_temperature_K"],
            minus["maximum_electron_temperature_K"],
        ),
    }


def _run_branch(
    sim,
    prepared,
    model,
    scales,
    seed_m,
    dt,
    exposure_steps,
    sample_every,
    mode,
    sign,
    offset_vector_m,
    absorbed_power_W,
    spot_sigma_m,
):
    sim.fields = deepcopy(prepared)
    sim.config.laser.enabled = sign is not None

    if sign is not None:
        target = seed_m + offset_vector_m * sign
        sim.config.laser.absorbed_power_W = absorbed_power_W
        sim.config.laser.sigma_m = spot_sigma_m
        sim.config.laser.waypoints = [
            LaserWaypointConfig(0.0, *target),
            LaserWaypointConfig(exposure_steps * dt, *target),
        ]
        sim.config.laser.validate()

    rows = []
    critical_temperature = float(np.min(sim.material_map.Tc))
    for step in range(exposure_steps + 1):
        if step:
            _, result = coupled_step(
                sim, dt=dt, tdgl_model=model, time_s=(step - 1) * dt
            )
            if not result.converged:
                raise RuntimeError(f"{mode} failed at step {step}.")
        overheated = bool(np.max(sim.fields.temperature) >= critical_temperature)
        if step % sample_every == 0 or step == exposure_steps or overheated:
            rows.append(
                _snapshot(
                    sim,
                    model,
                    scales,
                    seed_m,
                    step * dt,
                    step,
                    "exposure",
                    mode,
                    offset_vector_m=offset_vector_m if sign is not None else None,
                    absorbed_power_W=absorbed_power_W if sign is not None else None,
                    spot_sigma_m=spot_sigma_m if sign is not None else None,
                )
            )
        if overheated:
            print(f"{mode}: stopped at step {step} after electron temperature reached Tc",
                  flush=True)
            break
    return rows, deepcopy(sim.fields)


def _plot_sweep(summary_rows, output, dpi):
    sigmas = sorted({row["spot_sigma_m"] for row in summary_rows})
    powers = sorted({row["absorbed_power_W"] for row in summary_rows})

    for sigma in sigmas:
        fig, axis = plt.subplots(figsize=(9, 6), constrained_layout=True)
        for power in powers:
            subset = [
                row
                for row in summary_rows
                if row["spot_sigma_m"] == sigma and row["absorbed_power_W"] == power
            ]
            subset.sort(key=lambda row: row["offset_distance_nm"])
            if not subset:
                continue
            x = [row["offset_distance_nm"] for row in subset]
            y = [row["mean_motion_toward_beam_nm"] for row in subset]
            axis.plot(x, y, marker="o", label=f"{power * 1e9:g} nW")

        axis.axhline(0.0, linewidth=1)
        axis.set_xlabel("Beam-vortex offset (nm)")
        axis.set_ylabel("Mean vortex motion toward beam (nm)")
        axis.set_title(f"Optical-gradient response, sigma = {sigma * 1e9:g} nm")
        axis.grid(alpha=0.3)
        axis.legend(title="Absorbed power")
        filename = f"distance_vs_offset_sigma_{sigma * 1e9:g}".replace(".", "p") + "_nm.png"
        fig.savefig(output / filename, dpi=dpi)
        plt.close(fig)

    # If several widths are swept, also make one comparison plot at each power.
    if len(sigmas) > 1:
        for power in powers:
            fig, axis = plt.subplots(figsize=(9, 6), constrained_layout=True)
            for sigma in sigmas:
                subset = [
                    row
                    for row in summary_rows
                    if row["spot_sigma_m"] == sigma and row["absorbed_power_W"] == power
                ]
                subset.sort(key=lambda row: row["offset_distance_nm"])
                if not subset:
                    continue
                axis.plot(
                    [row["offset_distance_nm"] for row in subset],
                    [row["mean_motion_toward_beam_nm"] for row in subset],
                    marker="o",
                    label=f"sigma={sigma * 1e9:g} nm",
                )
            axis.axhline(0.0, linewidth=1)
            axis.set_xlabel("Beam-vortex offset (nm)")
            axis.set_ylabel("Mean vortex motion toward beam (nm)")
            axis.set_title(f"Spot-width comparison, absorbed power = {power * 1e9:g} nW")
            axis.grid(alpha=0.3)
            axis.legend(title="Beam width")
            filename = f"distance_vs_offset_power_{power * 1e9:g}".replace(".", "p") + "_nW.png"
            fig.savefig(output / filename, dpi=dpi)
            plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("tools/config/optical_gradient_probe.json"),
    )
    parser.add_argument("--dark-steps", type=int, help="Short smoke-test override.")
    parser.add_argument("--exposure-steps", type=int, help="Short smoke-test override.")
    parser.add_argument("--max-cases", type=int,
                        help="Stop after this many total sweep combinations.")
    parser.add_argument("--resume", type=Path,
                        help="Continue a prior sweep directory with the same config and step counts.")
    args = parser.parse_args(argv)

    config = json.loads(args.config.read_text(encoding="utf-8"))
    output = (reserve_output_directory(config["output_directory"])
              if args.resume is None else args.resume)
    offsets, powers, sigmas = _as_sweep(config)

    sim = build_simulation(config["simulation_config"])
    sim.config.current = 0.0
    sim.config.electrical.drive_mode = "current"
    sim.config.laser.enabled = False
    sim.config.thermal.model = "two_temperature"
    sim.config.thermal.validate()
    apply_uniform_magnetic_field(sim, config["applied_field_T"])
    _seed_initial_state(sim, {"vortices": [config["initial_vortex"]]})

    seed_m = np.asarray(config["seed_position_m"], dtype=float)
    dt = sim.config.dt
    model = build_tdgl_model(sim)
    scales = tdgl_scales(sim, model)

    dark_steps = (
        config["timing"]["dark_steps"] if args.dark_steps is None else args.dark_steps
    )
    exposure_steps = (
        config["timing"]["exposure_steps"]
        if args.exposure_steps is None
        else args.exposure_steps
    )
    sample_every = int(config["timing"]["sample_every"])
    if dark_steps < 0 or exposure_steps < 1:
        raise ValueError("Dark steps must be nonnegative; exposure steps must be positive.")
    if sample_every < 1:
        raise ValueError("timing.sample_every must be positive.")

    # One shared dark preparation. Every subsequent branch clones this exact state.
    dark_change = []
    for step in range(dark_steps):
        before = sim.fields.psi.copy()
        _, result = coupled_step(sim, dt=dt, tdgl_model=model, time_s=step * dt)
        if not result.converged:
            raise RuntimeError(f"Dark preparation failed at step {step + 1}.")
        dark_change.append(float(np.max(np.abs(sim.fields.psi - before))))
    prepared = deepcopy(sim.fields)

    # The dark exposure control is independent of beam parameters, so run it once.
    dark_rows, dark_fields = _run_branch(
        sim,
        prepared,
        model,
        scales,
        seed_m,
        dt,
        exposure_steps,
        sample_every,
        "dark",
        None,
        offsets[0],
        powers[0],
        sigmas[0],
    )
    dark_final = dark_rows[-1]

    if args.resume is None:
        all_rows = list(dark_rows)
        sweep_rows = []
    else:
        prior = json.loads((output/"summary.json").read_text(encoding="utf-8"))
        if ({key: value for key, value in prior["config"].items()
             if key != "_documentation"} !=
            {key: value for key, value in config.items() if key != "_documentation"}
                or prior["effective_dark_steps"] != dark_steps
                or prior["effective_exposure_steps"] != exposure_steps):
            raise ValueError("Resume settings differ from the saved sweep.")
        sweep_rows = prior["sweep_results"]
        with (output/"force_probe.csv").open(newline="", encoding="utf-8") as stream:
            all_rows = list(csv.DictReader(stream))
        completed_keys = {
            (row["offset_x_m"], row["offset_y_m"], row["absorbed_power_W"],
             row["spot_sigma_m"]) for row in sweep_rows}
        all_rows = [row for row in all_rows if row["mode"] == "dark" or
                    (float(row["offset_x_m"]), float(row["offset_y_m"]),
                     float(row["absorbed_power_W"]), float(row["spot_sigma_m"]))
                    in completed_keys]
    best_case = None
    best_score = max((row["mean_motion_toward_beam_nm"] for row in sweep_rows
                      if row.get("topology_preserved_throughout")
                      and row.get("peak_electron_temperature_throughout_K", float("inf"))
                      < float(np.min(sim.material_map.Tc))), default=-float("inf"))
    total_cases = len(offsets) * len(powers) * len(sigmas)
    requested_cases = min(total_cases, args.max_cases) if args.max_cases is not None else total_cases
    if requested_cases < 0:
        raise ValueError("--max-cases must be nonnegative.")
    print(f"Planned: {requested_cases} cases, about "
          f"{dark_steps + exposure_steps*(1+2*requested_cases)} coupled steps "
          "before early thermal cutoffs.", flush=True)
    case_number = 0
    threshold_nm = float(config.get("minimum_directional_response_nm", 0.0))
    critical_temperature = float(np.min(sim.material_map.Tc))

    for sigma in sigmas:
        for power in powers:
            for offset_vector in offsets:
                case_number += 1
                if args.max_cases is not None and case_number > args.max_cases:
                    break
                if case_number <= len(sweep_rows):
                    continue
                print(
                    f"case {case_number}/{total_cases}: "
                    f"offset={np.linalg.norm(offset_vector) * 1e9:g} nm, "
                    f"power={power * 1e9:g} nW, sigma={sigma * 1e9:g} nm",
                    flush=True,
                )

                plus_rows, plus_fields = _run_branch(
                    sim,
                    prepared,
                    model,
                    scales,
                    seed_m,
                    dt,
                    exposure_steps,
                    sample_every,
                    "plus",
                    +1,
                    offset_vector,
                    power,
                    sigma,
                )
                minus_rows, minus_fields = _run_branch(
                    sim,
                    prepared,
                    model,
                    scales,
                    seed_m,
                    dt,
                    exposure_steps,
                    sample_every,
                    "minus",
                    -1,
                    offset_vector,
                    power,
                    sigma,
                )
                all_rows.extend(plus_rows)
                all_rows.extend(minus_rows)

                final = {
                    "dark": dark_final,
                    "plus": plus_rows[-1],
                    "minus": minus_rows[-1],
                }
                response = response_summary(
                    final,
                    offset_vector,
                    critical_temperature,
                    threshold_nm,
                )
                response["absorbed_power_W"] = power
                response["absorbed_power_nW"] = power * 1e9
                response["spot_sigma_m"] = sigma
                response["spot_sigma_nm"] = sigma * 1e9
                response["plus_completed_steps"] = plus_rows[-1]["step"]
                response["minus_completed_steps"] = minus_rows[-1]["step"]
                response["fully_evolved"] = bool(
                    plus_rows[-1]["step"] == exposure_steps
                    and minus_rows[-1]["step"] == exposure_steps)
                branch_rows = plus_rows + minus_rows
                response["topology_preserved_throughout"] = all(
                    row["positive_vortices"] == 1 and row["negative_vortices"] == 0
                    for row in branch_rows)
                response["peak_electron_temperature_throughout_K"] = max(
                    row["maximum_electron_temperature_K"] for row in branch_rows)
                response["directional_response_in_model"] = bool(
                    response["directional_response_in_model"]
                    and response["fully_evolved"]
                    and response["topology_preserved_throughout"]
                    and response["peak_electron_temperature_throughout_K"] < critical_temperature)
                sweep_rows.append(response)

                if (response["topology_preserved_throughout"]
                        and response["peak_electron_temperature_throughout_K"] < critical_temperature
                        and np.isfinite(response["mean_motion_toward_beam_nm"])
                        and response["mean_motion_toward_beam_nm"] > best_score):
                    best_case = {"response": response, "rows": dark_rows+plus_rows+minus_rows,
                                 "fields": {"dark": dark_fields, "plus": plus_fields,
                                            "minus": minus_fields}}
                    best_score = response["mean_motion_toward_beam_nm"]

                _write_csv(output / "force_probe.csv", all_rows)
                _write_csv(output / "sweep_summary.csv", sweep_rows)
                (output / "summary.json").write_text(json.dumps({
                    "config": config,
                    "effective_dark_steps": dark_steps,
                    "effective_exposure_steps": exposure_steps,
                    "dark_preparation_last_psi_change": dark_change[-1] if dark_change else None,
                    "dark_final": dark_final,
                    "planned_cases": total_cases,
                    "completed_cases": len(sweep_rows),
                    "sweep_results": sweep_rows,
                }, indent=2), encoding="utf-8")

    if sweep_rows:
        _plot_sweep(sweep_rows, output, int(config["output"]["dpi"]))
    if best_case is not None:
        best_output = output/"best_case"
        best_output.mkdir(exist_ok=True)
        _plot_case_results(best_case["rows"], best_case["fields"], best_output,
                           config, seed_m, sim.mesh)
        _write_csv(best_output/"force_probe.csv", best_case["rows"])
        (best_output/"summary.json").write_text(
            json.dumps(best_case["response"], indent=2), encoding="utf-8")

    summary = {
        "config": config,
        "effective_dark_steps": dark_steps,
        "effective_exposure_steps": exposure_steps,
        "dark_preparation_last_psi_change": dark_change[-1] if dark_change else None,
        "dark_preparation_max_psi_change": max(dark_change) if dark_change else None,
        "dark_final": dark_final,
        "sweep_results": sweep_rows,
        "planned_cases": total_cases,
        "completed_cases": len(sweep_rows),
        "best_current_case": (max((row for row in sweep_rows
                                   if row.get("topology_preserved_throughout")
                                   and row.get("peak_electron_temperature_throughout_K", float("inf"))
                                   < critical_temperature),
                                  key=lambda row: row["mean_motion_toward_beam_nm"],
                                  default=None)),
        "interpretation": (
            "The primary sweep metric is mean_motion_toward_beam_nm: the mean of "
            "the +offset and -offset complex-zero displacements after subtracting "
            "the same dark-control drift and projecting each branch onto its own "
            "beam direction. Positive values mean motion toward the hotspot."
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(output)
    return output


if __name__ == "__main__":
    main()
