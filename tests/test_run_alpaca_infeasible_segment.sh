#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="tests/infeasible_segment/input/${tumour_id}"
output_directory="tests/infeasible_segment/output/${tumour_id}"

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --overwrite_output 0 \
    --debug \
    --simulate_infeasibility 1_762496_28527452