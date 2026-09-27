import importlib.util
import sys
from pathlib import Path

import numpy as np


TOOLS = Path(__file__).resolve().parents[1]/"tools"
sys.path.insert(0, str(TOOLS))
SPEC = importlib.util.spec_from_file_location(
    "optical_gradient_probe_sweep", TOOLS/"optical_gradient_probe_sweep.py")
SWEEP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SWEEP)


def _row(x, gx):
    return {
        "vortex_x_zero_m": x, "vortex_y_zero_m": 50e-9,
        "vortex_x_subcell_m": x, "vortex_y_subcell_m": 50e-9,
        "electron_gradient_x_at_core_K_per_m": gx,
        "electron_gradient_y_at_core_K_per_m": 0.0,
        "positive_vortices": 1, "negative_vortices": 0,
        "maximum_electron_temperature_K": 14.5,
    }


def test_sweep_projects_both_mirrored_offsets_toward_beam():
    final = {"dark": _row(50e-9, 0),
             "plus": _row(50.6e-9, 1e7),
             "minus": _row(49.4e-9, -1e7)}
    result = SWEEP.response_summary(final, [10e-9, 0], 15.5, 0.1)
    assert result["directional_response_in_model"]
    assert np.isclose(result["mean_motion_toward_beam_nm"], .6)


def test_sweep_rejects_empty_axes():
    try:
        SWEEP._as_sweep({"sweeps": {"beam_offsets_m": [],
                                     "absorbed_powers_W": [1e-7],
                                     "spot_sigmas_m": [1e-8]}})
    except ValueError:
        pass
    else:
        raise AssertionError("empty sweep should fail")

