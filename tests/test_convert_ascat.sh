#!/bin/bash
tumour_id=LTX0000-Tumour1
SCRIPT_DIR="alpaca/scripts/submodules/alpaca_input_formatting/convert_ascat_output"
output_dir="tests/test_ascat_conversion"

ascat_rds_path="dev/ASCAT_test_data/patient_1/ASCAT_objects.Rdata"

Rscript "${SCRIPT_DIR}/unpack_RDS_ascat_helper.R" \
--rdata $ascat_rds_path \
--output_dir $output_dir \
--out_prefix $tumour_id

ascat_segments_path="$output_dir/${tumour_id}_segments.tsv" 
ascat_snps_path="$output_dir/${tumour_id}_snps.tsv"
ascat_purity_ploidy_path="$output_dir/${tumour_id}_purity_ploidy.tsv"

python3 "${SCRIPT_DIR}/convert_ascat_output.py" \
    --tumour_id $tumour_id \
    --output_dir $output_dir \
    --segments $ascat_segments_path \
    --snps $ascat_snps_path \
    --purity_ploidy $ascat_purity_ploidy_path \
    --n_boot 10 \
    --gamma 1
