import numpy as np
import pytest

from shs.config.builder import build_simulation
from shs.config.simulation import PinningSiteConfig
from shs.physics.pinning import pinning_suppression
from shs.solvers.tdgl_solver import tdgl_step
from shs.tdgl.model import TDGLModel
from shs.tdgl.parameters import TDGLParameters


def _configure_site(simulation, *, strength=0.4, sigma=2e-9):
    simulation.config.pinning.enabled = True
    simulation.config.pinning.maximum_suppression = 0.5
    simulation.config.pinning.sites = [PinningSiteConfig(
        x_m=2.5e-8, y_m=2.5e-8, sigma_m=sigma, strength=strength,
    )]
    simulation.config.pinning.validate()


def test_gaussian_pinning_map_is_centered_bounded_and_static():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    _configure_site(simulation, strength=0.8)
    suppression = pinning_suppression(simulation)
    center = (simulation.mesh.ny // 2, simulation.mesh.nx // 2)
    assert np.max(suppression) == pytest.approx(0.5)
    assert suppression[center] == pytest.approx(0.5)
    assert np.all(suppression[~simulation.region_map.active_mask] == 0.0)
    assert pinning_suppression(simulation) is suppression


def test_pinning_defect_locally_suppresses_condensate_dynamics():
    simulation = build_simulation("configs/simulations/nbn_tdgl_transport.json")
    _configure_site(simulation)
    simulation.fields.temperature.fill(14.0)
    simulation.fields.psi.fill(0.2 + 0j)
    model = TDGLModel(TDGLParameters(
        normalization="pytdgl", temperature_model="tc_over_t_minus_one",
        max_normalized_timestep=0.01,
    ))
    tdgl_step(simulation, 1e-14, model)
    amplitude = np.abs(simulation.fields.psi)
    center = (simulation.mesh.ny // 2, simulation.mesh.nx // 2)
    assert amplitude[center] < amplitude[center[0], center[1] + 15]
