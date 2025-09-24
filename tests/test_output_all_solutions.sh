#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="examples/example_cohort/input/${tumour_id}"
output_directory="tests/all_solutions/output/${tumour_id}"
mkdir -p $output_directory

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --overwrite_output 1 \
    --output_all_solutions
  