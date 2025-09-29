#!/bin/bash
input_tumour_directory="tests/solution_provided/input/LTX0000-Tumour1"
output_directory="tests/solution_provided/output"

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --overwrite_output 1 \
    --debug \
    --debug_solution_file tests/solution_provided/input/LTX0000-Tumour1/tetraploid_solution.csv \
    --complexity 100 \
    --gurobi_logs "${output_directory}"
