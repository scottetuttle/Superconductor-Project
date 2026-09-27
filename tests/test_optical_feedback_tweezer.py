import importlib.util
import sys
from pathlib import Path

import numpy as np


TOOLS = Path(__file__).resolve().parents[1]/"tools"
sys.path.insert(0, str(TOOLS))
SPEC = importlib.util.spec_from_file_location(
    "optical_feedback_tweezer", TOOLS/"optical_feedback_tweezer.py")
FEEDBACK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FEEDBACK)


def test_beam_command_respects_speed_limit_and_target():
    current = np.array([50e-9, 60e-9])
    target = np.array([50e-9, 70e-9])
    next_beam = FEEDBACK.bounded_command(current, target, 1e-9)
    assert np.allclose((next_beam-current)/1e-9, [0, 1], atol=1e-12)
    assert np.allclose(FEEDBACK.bounded_command(next_beam, target, 20e-9), target)


def test_vortex_tracking_can_preserve_identity_when_an_extra_vortex_appears(monkeypatch):
    class Fields:
        psi = np.ones((4, 5), dtype=complex)
        vector_potential_x = np.zeros((4, 5))
        vector_potential_y = np.zeros((4, 5))

    class Mesh:
        dx = 1e-9
        dy = 1e-9
        x = np.arange(5)*dx
        y = np.arange(4)*dy

    class Scales:
        xi = 1e-9

        @staticmethod
        def vector_potential_to_dimensionless(value):
            return value

    class Detection:
        winding = np.array([[1, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
        positive_count = 2
        negative_count = 0

    sim = type("Simulation", (), {"fields": Fields(), "mesh": Mesh()})()
    monkeypatch.setattr(FEEDBACK, "detect_vortices", lambda *args: Detection())
    monkeypatch.setattr(FEEDBACK, "complex_zero_core", lambda psi, x, y: (x, y))

    reference = np.array([3.4e-9, .5e-9])
    assert FEEDBACK.vortex_center(sim, Scales(), reference) is None
    tracked = FEEDBACK.vortex_center(
        sim, Scales(), reference, require_unique=False, maximum_jump_m=1e-9)
    assert np.allclose(tracked, [3.5e-9, .5e-9])


def test_vortex_tracking_rejects_an_identity_jump(monkeypatch):
    class Fields:
        psi = np.ones((3, 3), dtype=complex)
        vector_potential_x = np.zeros((3, 3))
        vector_potential_y = np.zeros((3, 3))

    class Mesh:
        dx = dy = 1e-9
        x = y = np.arange(3)*1e-9

    class Scales:
        xi = 1e-9

        @staticmethod
        def vector_potential_to_dimensionless(value):
            return value

    class Detection:
        winding = np.array([[0, 0], [0, 1]])
        positive_count = 1
        negative_count = 0

    sim = type("Simulation", (), {"fields": Fields(), "mesh": Mesh()})()
    monkeypatch.setattr(FEEDBACK, "detect_vortices", lambda *args: Detection())
    monkeypatch.setattr(FEEDBACK, "complex_zero_core", lambda psi, x, y: (x, y))
    assert FEEDBACK.vortex_center(
        sim, Scales(), np.array([0.0, 0.0]), require_unique=False,
        maximum_jump_m=.5e-9) is None


def test_core_locator_falls_back_when_complex_zero_is_ill_conditioned(monkeypatch):
    psi = np.ones((7, 7), dtype=complex)
    monkeypatch.setattr(FEEDBACK, "complex_zero_core",
                        lambda *args: (float("nan"), float("nan")))
    monkeypatch.setattr(FEEDBACK, "subcell_core",
                        lambda *args, **kwargs: (3.2, 3.7))
    core, method = FEEDBACK.robust_core_position(psi, 3, 3)
    assert np.allclose(core, [3.2, 3.7])
    assert method == "amplitude_minimum"
