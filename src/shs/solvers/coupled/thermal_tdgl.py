"""
Thermal-TDGL coupled solver.

Coordinates the evolution of the superconducting order
parameter and thermal field.
"""

from shs.config.simulation_state import Simulation

from shs.tdgl.model import TDGLModel
from shs.physics.thermal import ThermalModel
from shs.physics.electrical import ElectricalModel

from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.thermal_solver import thermal_step

from shs.solvers.electrical_solver import electrical_step

import numpy as np


def thermal_tdgl_step(
    simulation,
    dt,
    tdgl_model,
    thermal_model,
    electrical_model,
):

    # dt is physical seconds

    tdgl_step(
        simulation,
        dt=dt,
        tdgl_model=tdgl_model,
    )

    fields = simulation.fields

    superconducting_fraction = (
        np.abs(fields.psi) ** 2
    )

    electrical_step(
        fields=fields,
        mesh=simulation.mesh,
        material_map=simulation.material_map,
        contact_map=simulation.contact_map,
        superconducting_current_x=(
            fields.supercurrent_density_x
        ),
        superconducting_current_y=(
            fields.supercurrent_density_y
        ),
        superconducting_fraction=(
            superconducting_fraction
        ),
    )

    thermal_step(
        simulation,
        dt=dt,
        thermal_model=thermal_model,
    )

    return fields