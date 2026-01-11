#!/bin/bash

# compare 4 implementations using a realistic test segment. NB cbc and HiGHS do not work here.

tumour_id=X0001_average
input_tumour_directory="tests/compare_solvers/input/${tumour_id}"


# run alpaca with default backend:
output_directory="tests/compare_solvers/output_gurobi_backend_all_solutions/${tumour_id}"

alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --output_all_solutions \
    --debug \
    --solver_logs "${output_directory}"

# run alpaca with gurobi backend mediated via pyomo:
output_directory="tests/compare_solvers/output_pyomo_backend_gurobi_solver_all_solutions/${tumour_id}"
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --solver pyomo \
    --pyomo_solver gurobi \
    --output_all_solutions \
    --debug

# use glpk via pyomo

output_directory="tests/compare_solvers/output_pyomo_backend_glpk_solver_all_solutions/${tumour_id}"
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --solver pyomo \
    --pyomo_solver glpk \
    --output_all_solutions \
    --debug 

# use scip via pyomo

output_directory="tests/compare_solvers/output_pyomo_backend_scip_solver_all_solutions/${tumour_id}"
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --solver pyomo \
    --pyomo_solver scip \
    --debug \
    --output_all_solutions \
    --solver_logs "${output_directory}"


