"""Sequential thermal/electrical/TDGL integration compatibility interface."""
from copy import copy, deepcopy
from shs.solvers.tdgl_solver import tdgl_step
from shs.solvers.electrical_solver import electrical_simulation_step
from shs.solvers.thermal_solver import thermal_step


def thermal_tdgl_step(simulation, dt, tdgl_model, thermal_model,
                      electrical_model=None):
    work = copy(simulation)
    work.fields = deepcopy(simulation.fields)
    tdgl_step(work, dt, tdgl_model)
    electrical_simulation_step(work)
    thermal_step(work, dt, thermal_model)
    simulation.fields.__dict__.update(vars(work.fields))
    return simulation.fields
