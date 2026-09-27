import importlib.util
import sys
from pathlib import Path

import numpy as np


TOOLS = Path(__file__).resolve().parents[1]/"tools"
sys.path.insert(0, str(TOOLS))
SPEC = importlib.util.spec_from_file_location("optical_gradient_probe",
                                              TOOLS/"optical_gradient_probe.py")
PROBE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PROBE)


def test_subcell_core_recovers_quadratic_minimum():
    yy, xx = np.indices((30, 30))
    amplitude = np.sqrt(.1 + (xx-14.3)**2 + 2*(yy-12.7)**2)
    x, y = PROBE.subcell_core(amplitude, 14.5, 12.5)
    assert np.allclose([x, y], [14.3, 12.7], atol=1e-10)


def test_sample_interpolates_linear_field():
    yy, xx = np.indices((10, 10))
    assert np.isclose(PROBE._sample(2*xx+3*yy, 4.25, 2.5), 16.0)


def test_complex_zero_tracks_true_core_under_amplitude_gradient():
    yy, xx = np.indices((30, 30))
    psi = (xx-14.3 + 1j*(yy-12.7))*(1+.02*yy)
    x, y = PROBE.complex_zero_core(psi, 14.5, 12.5)
    assert np.allclose([x, y], [14.3, 12.7], atol=1e-2)
