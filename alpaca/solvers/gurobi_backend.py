"""Gurobi-backed implementation of the ALPACA solver interface."""

from __future__ import annotations

from typing import Dict, Optional

# get Gurobi backend from the original implementation
from alpaca.ALPACA_model_class import Model as GurobiModel

from .base import ModelInputs, SolverBackend, SolverResult


class GurobiBackend(SolverBackend):
    name = "gurobi"

    def solve(self, allowed_complexity: int) -> SolverResult:
        model_kwargs = {**self.config, "allowed_tree_complexity": allowed_complexity}
        model = GurobiModel(
            segment=self.inputs.segment,
            ci_table=self.inputs.ci_table,
            fractional_copy_number_table=self.inputs.fractional_copy_number_table,
            tree=self.inputs.tree,
            clone_proportions=self.inputs.clone_proportions,
            **model_kwargs,
        )
        model.model.optimize()
        model.get_output()
        runtime = getattr(model.model, "Runtime", None)
        objective_values: Dict[str, float] = {}
        if getattr(model, "add_CI_objective", False):
            objective_values["CI"] = float(model.Z.getValue())
        if getattr(model, "add_D_objective", False):
            objective_values["D"] = float(model.D.getValue())
        solver_specific = {
            "gurobi_model": model.model,
            "gurobi_solution": model.solution,
        }
        return SolverResult(
            solution=model.solution.copy(),
            runtime=runtime,
            backend_name=self.name,
            gap_status=getattr(model, "gap_status", {}),
            raw_model=model.model,
            wrapper=model,
            objective_values=objective_values,
            solver_specific=solver_specific,
        )

    def supports_feature(self, feature_name: str) -> bool:
        if feature_name in {"iis", "objective_gap"}:
            return True
        return super().supports_feature(feature_name)
