#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="tests/x_chromosome/input/${tumour_id}"
output_directory="tests/x_chromosome/output/${tumour_id}"

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --overwrite_output 0 \
    --debug