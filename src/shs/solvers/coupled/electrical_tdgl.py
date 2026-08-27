"""
Electro-TDGL coupled solver.

Coordinates the interaction between:

    TDGL order parameter
        ↓
    superconducting current
        ↓
    electrical transport
        ↓
    electric field
        ↓
    dissipative normal current
"""

from shs.config.simulation_state import Simulation

from shs.tdgl.model import TDGLModel
from shs.physics.electrical import ElectricalModel

from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.electrical_solver import electrical_step


def electrical_tdgl_step(
    simulation: Simulation,
    dt: float,
    tdgl_model: TDGLModel,
    electrical_model: ElectricalModel,
    voltage_left: float = 1.0,
    voltage_right: float = 0.0,
):
    """
    Advance the coupled electrical-TDGL system by one timestep.

    Coupling:

        T, A → TDGL → psi
                     ↓
                   Js
                     ↓
              electrical solve
                     ↓
                    E
                     ↓
                   Jn
    """

    fields = simulation.fields

    # --------------------------------------------------
    # 1. Advance TDGL
    # --------------------------------------------------

    tdgl_step(
        simulation,
        dt=dt,
        tdgl_model=tdgl_model,
    )

    # --------------------------------------------------
    # 2. Obtain superconducting current
    # --------------------------------------------------

    superconducting_current_x = (
        fields.supercurrent_density_x
    )

    superconducting_current_y = (
        fields.supercurrent_density_y
    )

    # --------------------------------------------------
    # 3. Solve electrical transport
    # --------------------------------------------------

    electrical_step(
        fields=fields,
        mesh=simulation.mesh,
        material_map=simulation.material_map,
        contact_map=simulation.contact_map,
        voltage_left=voltage_left,
        voltage_right=voltage_right,
        superconducting_current_x=superconducting_current_x,
        superconducting_current_y=superconducting_current_y,
    )

    return fields