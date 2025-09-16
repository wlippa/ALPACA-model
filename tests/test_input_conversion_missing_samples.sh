#!/bin/bash
tumour_id=LTX0000-Tumour1
output_dir="test_output/missing_samples"
CONIPHER_tree_object="missing_samples/input/${tumour_id}/${tumour_id}.tree.RDS"
refphase_rData="missing_samples/input/${tumour_id}/${tumour_id}-refphase-results.RData"
CONIPHER_tree_index=1
SCRIPT_DIR="../alpaca/scripts/submodules/alpaca_input_formatting"

bash "${SCRIPT_DIR}/input_conversion.sh" \
    --refphase_rData $refphase_rData \
    --tumour_id $tumour_id \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --CONIPHER_tree_index $CONIPHER_tree_index \
    --output_dir $output_dir

python - <<EOF
import pandas as pd
output_dir="${output_dir}"
cp_table = pd.read_csv(f"{output_dir}/cp_table.csv", index_col="clone")
ci_table = pd.read_csv(f"{output_dir}/ci_table.csv")
cp_table_samples = set(cp_table.columns)
ci_table_samples = set(ci_table['sample'])
assert cp_table_samples == ci_table_samples, f"Samples in cp_table ({cp_table_samples}) do not match samples in ci_table ({ci_table_samples})"
EOF