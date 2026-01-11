"""Solver backend registry for ALPACA."""

from .base import ModelInputs, SolverResult, SolverBackend
from .factory import create_solver_backend, AVAILABLE_SOLVERS

__all__ = [
    "ModelInputs",
    "SolverResult",
    "SolverBackend",
    "create_solver_backend",
    "AVAILABLE_SOLVERS",
]
