import numpy as np

from shs.config.builder import build_simulation
from shs.solvers import run_coupled_simulation
from shs.tdgl import TDGLBoundarySide, TDGLBoundaryType


def test_small_step_tdgl_example_runs_one_coupled_step():
    simulation = build_simulation(
        "configs/simulations/nbn_tdgl_small_step.json"
    )

    assert simulation.config.tdgl.include_scalar_potential
    assert simulation.tdgl_boundaries.get(
        TDGLBoundarySide.LEFT
    ).type == TDGLBoundaryType.NORMAL_CONTACT

    result = run_coupled_simulation(simulation, steps=1)

    assert result.converged
    assert result.elapsed_time == simulation.config.dt
    assert np.all(np.isfinite(simulation.fields.psi))
    assert np.all(simulation.fields.psi[:, 0] == 0.0)
    assert np.all(simulation.fields.psi[:, -1] == 0.0)
