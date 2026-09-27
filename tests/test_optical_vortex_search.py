import importlib.util
import sys
from pathlib import Path

import numpy as np


TOOLS = Path(__file__).resolve().parents[1]/"tools"
sys.path.insert(0, str(TOOLS))
SPEC = importlib.util.spec_from_file_location("optical_vortex_search", TOOLS/"optical_vortex_search.py")
SEARCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SEARCH)


THRESHOLDS = {
    "minimum_beam_displacement_m": 5e-9,
    "minimum_vortex_displacement_m": 2e-9,
    "minimum_position_correlation": .75,
    "maximum_rms_lag_m": 8e-9,
    "maximum_initial_offset_m": 4e-9,
    "minimum_control_difference_m": 1e-9,
}


def test_follower_scores_only_persistent_directional_motion():
    start = np.array([40e-9, 25e-9])
    laser = [start + np.array([0, x]) for x in np.linspace(0, 8e-9, 9)]
    follower = [point + np.array([0, -1e-9]) for point in laser]
    drift = [start + np.array([x, 0]) for x in np.linspace(0, 8e-9, 9)]
    assert SEARCH.score_motion(follower, laser, start, [0, 1], THRESHOLDS)["qualified"]
    assert not SEARCH.score_motion(drift, laser, start, [0, 1], THRESHOLDS)["qualified"]
    assert not SEARCH.score_motion(follower[:-1], laser, start, [0, 1], THRESHOLDS)["qualified"]
    assert SEARCH.score_motion(follower[:2], laser[:2], start, [0, 1], THRESHOLDS)["reason"] == "insufficient transport samples"


def test_controls_reject_common_drift_and_accept_reversed_following():
    forward = {"qualified": True, "converged": True,
               "displacement_along_path_m": 6e-9}
    controls = {
        "no_laser": {"converged": True, "displacement_along_path_m": 0.5e-9},
        "stationary": {"converged": True, "displacement_along_path_m": 1e-9},
        "reverse": {"converged": True, "qualified": True},
    }
    assert SEARCH.controls_confirm(forward, controls, THRESHOLDS)
    controls["stationary"]["displacement_along_path_m"] = 5.5e-9
    assert not SEARCH.controls_confirm(forward, controls, THRESHOLDS)
    controls["stationary"]["displacement_along_path_m"] = 1e-9
    controls["reverse"]["qualified"] = False
    assert not SEARCH.controls_confirm(forward, controls, THRESHOLDS)
