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

import numpy as np


def thermal_tdgl_step(
    simulation: Simulation,
    dt: float,
    tdgl_model: TDGLModel,
    thermal_model: ThermalModel,
    electrical_model: ElectricalModel,
):
    """
    Advance the coupled TDGL-thermal system by one timestep.

    Coupling:

        T → TDGL → psi
        psi → normal fraction → J_n
        J_n · E → Joule heating
        Joule heating → T
    """

    # --------------------------------------------------
    # 1. Advance TDGL
    # --------------------------------------------------

    tdgl_step(
        simulation,
        dt=dt,
        tdgl_model=tdgl_model,
    )

    fields = simulation.fields
    material_map = simulation.material_map

    # --------------------------------------------------
    # 2. Calculate superconducting fraction
    # --------------------------------------------------

    superconducting_fraction = (
        np.abs(fields.psi) ** 2
    )

    # --------------------------------------------------
    # 3. Calculate normal current
    # --------------------------------------------------

    normal_current_x, normal_current_y = (
        electrical_model.normal_current(
            fields.electric_field_x,
            fields.electric_field_y,
            material_map.normal_resistivity,
            superconducting_fraction,
        )
    )

    # --------------------------------------------------
    # 4. Calculate total current
    # --------------------------------------------------

    total_current_x, total_current_y = (
        electrical_model.total_current(
            fields.supercurrent_density_x,
            fields.supercurrent_density_y,
            normal_current_x,
            normal_current_y,
        )
    )

    fields.current_density_x = total_current_x
    fields.current_density_y = total_current_y

    # --------------------------------------------------
    # 5. Calculate dissipative heating
    # --------------------------------------------------

    fields.heat_source = (
        electrical_model.joule_heating(
            normal_current_x,
            normal_current_y,
            fields.electric_field_x,
            fields.electric_field_y,
        )
    )

    # --------------------------------------------------
    # 6. Advance thermal field
    # --------------------------------------------------

    thermal_step(
        simulation,
        dt=dt,
        thermal_model=thermal_model,
    )

    return fields