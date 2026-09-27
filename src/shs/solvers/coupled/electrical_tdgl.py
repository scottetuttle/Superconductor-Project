"""Sequential electrical/TDGL integration compatibility interface."""
from copy import copy, deepcopy
from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.electrical_solver import electrical_simulation_step


def electrical_tdgl_step(simulation, dt, tdgl_model, electrical_model=None,
                         voltage_left=None, voltage_right=None):
    work = copy(simulation)
    work.fields = deepcopy(simulation.fields)
    tdgl_step(work, dt, tdgl_model)
    electrical_simulation_step(work, voltage_left, voltage_right)
    simulation.fields.__dict__.update(vars(work.fields))
    return simulation.fields
