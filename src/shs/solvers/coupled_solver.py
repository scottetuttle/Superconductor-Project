"""
Coupled electrothermal solver.

Coordinates:

Electrical transport
        |
        v
Joule heating
        |
        v
Thermal evolution

The coupled solver does not contain physics equations.
It coordinates existing physics modules.
"""


from dataclasses import dataclass


from shs.config.simulation_state import Simulation

from shs.physics.thermal import ThermalModel

from shs.solvers.electrical_solver import electrical_step

from shs.solvers.thermal_solver import thermal_step



@dataclass
class CoupledSolverResult:
    """
    Stores coupled simulation results.
    """

    simulation: Simulation

    steps: int



def coupled_step(
    simulation: Simulation,
    dt: float,
    thermal_model: ThermalModel,
    voltage_left: float = 1.0,
    voltage_right: float = 0.0,
):
    """
    Perform one coupled electrothermal timestep.

    Sequence:

    1. Solve electrical transport
    2. Generate Joule heating
    3. Advance thermal state

    """



    #
    # Electrical solve
    #

    electrical_step(
        fields=simulation.fields,
        mesh=simulation.mesh,
        material_map=simulation.material_map,
        contact_map=simulation.contact_map,
        voltage_left=voltage_left,
        voltage_right=voltage_right,
    )



    #
    # Thermal solve
    #

    thermal_step(
        simulation=simulation,
        dt=dt,
        thermal_model=thermal_model,
    )



    return simulation



def run_coupled_simulation(
    simulation: Simulation,
    thermal_model: ThermalModel,
    steps: int,
    dt: float,
    voltage_left: float = 1.0,
    voltage_right: float = 0.0,
):
    """
    Run multiple electrothermal timesteps.
    """


    for _ in range(steps):

        coupled_step(
            simulation,
            dt,
            thermal_model,
            voltage_left,
            voltage_right,
        )


    return CoupledSolverResult(
        simulation=simulation,
        steps=steps,
    )