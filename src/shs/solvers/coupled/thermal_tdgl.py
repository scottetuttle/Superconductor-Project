"""
Thermal-TDGL coupled solver.

Coordinates the evolution of:

    superconducting order parameter
    thermal field

The individual physics models remain responsible
for their own equations.
"""

from shs.config.simulation_state import Simulation

from shs.tdgl.model import TDGLModel
from shs.physics.thermal import ThermalModel

from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.thermal_solver import thermal_step


def thermal_tdgl_step(
    simulation: Simulation,
    dt: float,
    tdgl_model: TDGLModel,
    thermal_model: ThermalModel,
):
    """
    Advance the coupled TDGL-thermal system by one timestep.

    Sequence:

        1. Advance TDGL using the current temperature.
        2. Advance the thermal field using the current heat source.

    Parameters
    ----------
    simulation:
        Complete simulation state.

    dt:
        Coupled timestep.

    tdgl_model:
        TDGL physics model.

    thermal_model:
        Thermal physics model.

    Returns
    -------
    Fields
        Updated simulation fields.
    """

    # TDGL evolution
    tdgl_step(
        simulation,
        dt=dt,
        tdgl_model=tdgl_model,
    )

    # Thermal evolution
    thermal_step(
        simulation,
        dt=dt,
        thermal_model=thermal_model,
    )

    return simulation.fields