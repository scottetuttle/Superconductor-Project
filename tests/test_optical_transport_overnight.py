import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))
SPEC = importlib.util.spec_from_file_location(
    "optical_transport_overnight", TOOLS / "optical_transport_overnight.py")
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
VIS_SPEC = importlib.util.spec_from_file_location(
    "visualize_optical_transport", TOOLS / "visualize_optical_transport.py")
VISUALIZER = importlib.util.module_from_spec(VIS_SPEC)
VIS_SPEC.loader.exec_module(VISUALIZER)


def test_thermal_hysteresis_and_dwell():
    switch = lambda on, temp, last, step: RUNNER.thermal_switch(
        on, temp, off_K=15.0, on_K=14.6, dwell=10,
        last_change=last, step=step)
    assert switch(True, 15.1, -10, 1) == (False, 1)
    assert switch(False, 14.5, 1, 5) == (False, 1)
    assert switch(False, 14.8, 1, 11) == (False, 1)
    assert switch(False, 14.5, 1, 11) == (True, 11)


def test_stall_requires_multiple_completed_windows():
    assert not RUNNER.stalled([0, 0.01, 0.02], 200, 100, 0.05, 3)
    assert RUNNER.stalled([0, 0.01, 0.02, 0.03], 300, 100, 0.05, 3)
    assert not RUNNER.stalled([0, 0.01, 0.2, 0.21], 300, 100, 0.05, 3)


def test_optional_microwave_contact_current():
    assert RUNNER.drive_current(
        0.0, enabled=False, dc_A=2, ac_A=3, frequency_Hz=5) == 0
    assert RUNNER.drive_current(
        .25, enabled=True, dc_A=2, ac_A=3, frequency_Hz=1) == pytest.approx(5)
    with pytest.raises(ValueError, match="positive frequency"):
        RUNNER.drive_current(
            0.0, enabled=True, dc_A=0, ac_A=1, frequency_Hz=0)


def test_visualizer_interpolates_isolated_legacy_coordinate_gap():
    points = np.array([[1.0, 4.0], [np.nan, np.nan], [3.0, 8.0]])
    filled = VISUALIZER._fill_coordinate_gaps(points, "vortex")
    assert np.allclose(filled, [[1, 4], [2, 6], [3, 8]])


def test_visualizer_rejects_trajectory_without_any_finite_coordinate():
    with pytest.raises(ValueError, match="no finite coordinate"):
        VISUALIZER._fill_coordinate_gaps(np.full((2, 2), np.nan), "vortex")


def test_resume_ignores_documentation_and_observational_controls():
    old = {"physics": {"field": 1}, "_documentation": {"purpose": "old"}}
    new = {"physics": {"field": 1}, "_documentation": {"purpose": "new"},
           "live_preview": {"enabled": True}, "profiling": {"enabled": True}}
    assert (RUNNER._resume_compatible_config(old)
            == RUNNER._resume_compatible_config(new))


def test_resume_still_detects_a_physics_change():
    first = {"physics": {"field": 1}, "live_preview": {"enabled": False}}
    second = {"physics": {"field": 2}, "live_preview": {"enabled": True}}
    assert (RUNNER._resume_compatible_config(first)
            != RUNNER._resume_compatible_config(second))
