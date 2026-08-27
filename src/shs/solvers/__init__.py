"""
SHS Solver Package

Contains numerical solvers for SHS physics modules.
"""

from .thermal_solver import thermal_step

from .electrical_solver import electrical_step

from .coupled_solver import (
    coupled_step,
    run_coupled_simulation,
    CoupledSolverResult,
)

from .tdgl_solver import tdgl_step

from .coupled import thermal_tdgl_step