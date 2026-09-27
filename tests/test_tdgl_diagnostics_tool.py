import runpy

import numpy as np

from shs.config.builder import build_simulation


DIAGNOSTICS = runpy.run_path("tools/tdgl_diagnostics.py", run_name="tdgl_diagnostics_test")


def test_vortex_detector_finds_signed_phase_winding():
    yy, xx = np.indices((8, 8))
    positive = np.exp(1j * np.arctan2(yy - 3.5, xx - 3.5))
    negative = np.conjugate(positive)

    positive_result = DIAGNOSTICS["detect_vortices"](positive)
    negative_result = DIAGNOSTICS["detect_vortices"](negative)

    assert (positive_result.positive_count, positive_result.negative_count) == (1, 0)
    assert (negative_result.positive_count, negative_result.negative_count) == (0, 1)


def test_vortex_detector_is_invariant_under_discrete_gauge_transform():
    yy, xx = np.indices((8, 8))
    psi = np.exp(1j * np.arctan2(yy - 3.5, xx - 3.5))
    chi = 0.03 * xx**2 + 0.07 * yy
    ax = np.zeros_like(xx, dtype=float)
    ay = np.zeros_like(yy, dtype=float)
    transformed_ax = ax.copy()
    transformed_ay = ay.copy()
    transformed_ax[:, :-1] += chi[:, 1:] - chi[:, :-1]
    transformed_ay[:-1, :] += chi[1:, :] - chi[:-1, :]

    original = DIAGNOSTICS["detect_vortices"](psi, ax, ay)
    transformed = DIAGNOSTICS["detect_vortices"](
        psi * np.exp(1j * chi), transformed_ax, transformed_ay
    )

    assert np.array_equal(original.winding, transformed.winding)


def test_diagnostic_sweeps_expand_cartesian_product():
    cases = DIAGNOSTICS["build_cases"](
        {"tdgl.gamma": [0.0, 2.0], "electrical.voltage_left": [1e-6, 2e-6]}
    )

    assert len(cases) == 4
    assert cases[0][1] == {"tdgl.gamma": 0.0, "electrical.voltage_left": 1e-6}


def test_reference_transport_passes_declared_tdgl_regime():
    config = DIAGNOSTICS["load_config"]("tools/config/tdgl_diagnostics.json")
    simulation = build_simulation(config["simulation_config"])

    assessment = DIAGNOSTICS["assess_model_regime"](
        simulation, config["validity"]
    )

    assert assessment["passes"]


def test_diagnostic_current_sweep_has_one_authoritative_control():
    config = DIAGNOSTICS["load_config"]("tools/config/tdgl_diagnostics.json")
    assert "transport" not in config
    assert config["sweeps"]["diagnostics.target_current_A"] == [1e-8]
    assert "sweeps.diagnostics.target_current_A" in config["_documentation"]["fields"]


def test_relative_phase_change_removes_global_rotation_but_keeps_structure():
    amplitude = np.ones((4, 5))
    reference = {"phase": np.zeros((4, 5)), "amplitude": amplitude}
    global_only = {"phase": np.full((4, 5), 0.4), "amplitude": amplitude}
    structured = {"phase": 0.4 + np.linspace(-0.1, 0.1, 5)[None, :]
                  + np.zeros((4, 1)), "amplitude": amplitude}

    removed = DIAGNOSTICS["relative_phase_change"](global_only, reference)
    retained = DIAGNOSTICS["relative_phase_change"](structured, reference)

    assert np.allclose(removed, 0.0)
    assert np.isclose(np.nanmax(retained) - np.nanmin(retained), 0.2)
