#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="tests/find_monoclonal_samples/input/${tumour_id}"
output_directory="tests/find_monoclonal_samples/output/${tumour_id}"

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --overwrite_output 0 \
    --debug