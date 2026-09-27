import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))
SPEC = importlib.util.spec_from_file_location(
    "optical_pinning_response_sweep", TOOLS / "optical_pinning_response_sweep.py")
SWEEP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SWEEP)


def test_dark_subtracted_velocity_and_gradient_direction(monkeypatch):
    mesh = SimpleNamespace(nx=10, ny=10, dx=1e-9, dy=1e-9)
    sim = SimpleNamespace(mesh=mesh, material_map=SimpleNamespace(Tc=np.full((10, 10), 15.5)))
    monkeypatch.setattr(SWEEP, "pinning_suppression", lambda _: np.zeros((10, 10)))

    def row(step, y, gradient):
        return {"step": step, "time_s": step * 1e-12,
                "vortex_x_zero_m": 5e-9, "vortex_y_zero_m": y,
                "electron_gradient_x_at_core_K_per_m": 0.0,
                "electron_gradient_y_at_core_K_per_m": gradient,
                "maximum_electron_temperature_K": 14.5,
                "positive_vortices": 1, "negative_vortices": 0}

    dark = [row(0, 5e-9, 0), row(1, 5.1e-9, 0)]
    illuminated = [row(0, 5e-9, 0), row(1, 5.3e-9, 2e7)]
    result = SWEEP._trajectory(illuminated, dark, np.array([0.0, 1.0]), sim)
    assert np.isclose(result[-1]["dark_subtracted_displacement_toward_beam_nm"], 0.2)
    assert np.isclose(result[-1]["dark_subtracted_velocity_toward_beam_m_per_s"], 200)
    assert np.isclose(result[-1]["thermal_gradient_toward_beam_K_per_nm"], 0.02)
