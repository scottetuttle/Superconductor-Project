"""
Numerical infrastructure for SHS.

Contains reusable mathematical operators
and numerical solvers used by physics modules.
"""

from .operators import (
    gradient,
    divergence,
    laplacian,
)

from .iterative import (
    gauss_seidel,
    red_black_sor,
    SolverResult,
)


__all__ = [
    "gradient",
    "divergence",
    "laplacian",
    "gauss_seidel",
    "SolverResult",
    "red_black_sor",
]