#!/usr/bin/env bash

# Fail fast: exit on error, unset variables are errors, propagate failures in pipelines
set -euo pipefail
IFS=$'\n\t'

# Print a helpful message and exit when any command fails
trap 'code=$?; cmd="${BASH_COMMAND:-unknown}"; echo "ERROR: command failed: \"$cmd\" (exit code $code)" >&2; exit $code' ERR

# Ensure script is run with bash (not sh)
if [ -z "${BASH_VERSION:-}" ]; then
    echo "This script requires bash. Run with 'bash $0' or make executable and run './$0'" >&2
    exit 2
fi
# default arguments

usage() {
    echo "Usage: $0 --tumour_id TUMOUR_ID --refphase_rData RDATA_PATH --CONIPHER_tree_object TREE_OBJECT_PATH --output_dir OUTPUT_DIR"
    echo
    echo "Arguments:"
    echo "  --tumour_id              Tumour ID (required)"
    echo "  --refphase_rData         Path to refphase .RData file (required)"
    echo "  --CONIPHER_tree_object   Path to CONIPHER tree object .RDS file (required)"
    echo "  --CONIPHER_tree_index    Selected CONIPHER tree index (optional, default: 1)"
    echo "  --output_dir             Output directory (required)"
    echo "  --help                   Display this help message"
    exit 1
}
# default arguments
CONIPHER_tree_index=1

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --tumour_id) tumour_id="$2"; shift ;;
        --refphase_rData) refphase_rData="$2"; shift ;;
        --CONIPHER_tree_object) CONIPHER_tree_object="$2"; shift ;;
        --CONIPHER_tree_index) CONIPHER_tree_index="$2"; shift ;;
        --output_dir) output_dir="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if required arguments are provided
if [ -z "$tumour_id" ] || [ -z "$refphase_rData" ] || [ -z "$CONIPHER_tree_object" ] || [ -z "$output_dir" ]; then
    echo "Error: All arguments are required"
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


echo "Tumour ID: ${tumour_id}"
echo "Extract TSV files from RData"
echo "make output directory: $output_dir"
mkdir -p $output_dir

# exit if refphase_rData does not exist:
if [ ! -f $refphase_rData ]; then
    echo "Error: refphase_rData (${refphase_rData}) file does not exist" >&2
    exit 1
fi

# Check if Rscript is available
if ! command -v Rscript &> /dev/null; then
    echo "Error: Rscript is not available in the environment." >&2
    echo "Please install R and ensure Rscript is in your PATH." >&2
    exit 1
fi
echo "===================================="
echo "Extracting data from REFPHASE output"
Rscript "${SCRIPT_DIR}/convert_refphase_output/extract_rephase_data.R" \
    --refphase_rData $refphase_rData \
    --output_dir $output_dir
if [ $? -ne 0 ]; then
    echo "Rscript extract_rephase_data.R failed" >&2
    exit 1
fi

refphase_segments_path="${output_dir}/phased_segs.tsv"
refphase_snps_path="${output_dir}/phased_snps.tsv"
refphase_purity_ploidy_path="${output_dir}/purity_ploidy.tsv"
echo "===================================="
echo "Extracting data from CONIPHER output"
Rscript "${SCRIPT_DIR}/convert_conipher_output/convert_conipher_output.R" \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --CONIPHER_tree_index $CONIPHER_tree_index \
    --output_dir $output_dir
if [ $? -ne 0 ]; then
    echo "Rscript convert_conipher_output.R failed" >&2
    exit 1
fi
echo "===================================="
echo "Converting REFPHASE output"
python3 "${SCRIPT_DIR}/convert_refphase_output/convert_refphase.py" \
    --tumour_id $tumour_id \
    --output_dir $output_dir \
    --refphase_segments $refphase_segments_path \
    --refphase_snps $refphase_snps_path \
    --refphase_purity_ploidy $refphase_purity_ploidy_path \
    --conipher_cp_table "${output_dir}/cp_table.csv" \

if [ $? -ne 0 ]; then
    echo "Python convert_refphase.py failed" >&2
    exit 1
fi
