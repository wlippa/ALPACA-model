import os
import subprocess
import sys
from pathlib import Path
from importlib.resources import files
from alpaca.analysis import get_cn_change_to_ancestor
from alpaca.analysis import calculate_ccd
import argparse
from datetime import datetime
import logging
from alpaca.utils import create_logger, save_dataframe_to_csv, SUPPORTED_GENOME_BUILDS
from alpaca.plotting import (
    export_plot_outputs,
    _DEFAULT_HEATMAP_PALETTE,
    _SUPPORTED_HEATMAP_CHOICES,
    _DEFAULT_GENOME_BUILD,
)


def input_conversion():
    """
    Wrapper function to execute input_conversion.sh from submodule
    """
    print("Running input_conversion - it may take a few minutes")
    try:
        # check if all file in sys.argv exist:
        for i, arg in enumerate(sys.argv[1:], start=1):
            if "=" not in arg and ("/" in arg or "\\" in arg or "." in arg):
                # check if argument is a number
                try:
                    float(arg)
                    continue  # it's a number, skip existence check
                except ValueError:
                    pass  # not a number, proceed to check existence
                exists = os.path.exists(arg)
                print(
                    f"Argument {i} ({arg}): {'Exists' if exists else 'DOES NOT EXIST'}"
                )
        # Locate the input_conversion.sh script
        script_path = str(
            files("alpaca").joinpath(
                "scripts/submodules/alpaca_input_formatting/input_conversion.sh"
            )
        )
        print(script_path)
        # Locate the submodules directory
        submodules_path = str(files("alpaca").joinpath("scripts/submodules"))
        # Set environment variable to help script locate its dependencies
        env = os.environ.copy()
        env["SUBMODULES_PATH"] = submodules_path
        cmd_args = sys.argv[1:]
        print(f"Command line arguments: {cmd_args}")
        cmd = [script_path] + sys.argv[2:]
        print(f"Running: {' '.join(cmd)}")
        try:
            result = subprocess.run(
                cmd,
                check=True,
                text=True,
                env=env,
                stdout=sys.stdout,
                stderr=sys.stderr,
            )
            print(f"Return code: {result.returncode}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing input_conversion (return code {e.returncode}).\nCommand: {' '.join(cmd)}", file=sys.stderr)
            return e.returncode
        print("")
        print("")
        print("---------------------------------------")
        print("Input conversion completed successfully.")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Error executing input_conversion: {e.stderr}", file=sys.stderr)
        return e.returncode
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        return 1


def run_get_cn_change_to_ancestor():
    """CLI wrapper for get_cn_change_to_ancestor"""
    logger = create_logger(name="ancestor-delta", log_dir="logs")
    parser = argparse.ArgumentParser(
        description="Compute copy number changes to ancestor and save to CSV."
    )
    parser.add_argument("command", choices=["ancestor-delta"], help="Command to run")
    parser.add_argument("--tree_path", help="Path to the tree file", required=True)
    parser.add_argument(
        "--tumour_df_path",
        help="Path to the tumour dataframe file (CSV format)",
        required=True,
    )
    parser.add_argument(
        "--output_directory",
        help="Directory to save the output CSV file",
        required=True,
    )

    args = parser.parse_args()
    # Validate input files exist
    if not os.path.isfile(args.tree_path):
        logger.error(f"Tree file not found: {args.tree_path}")
        exit(1)

    if not os.path.isfile(args.tumour_df_path):
        logger.error(f"Tumour dataframe file not found: {args.tumour_df_path}")
        exit(1)

    try:
        logger.info("Starting analysis...")
        cn_change_to_ancestor_df = get_cn_change_to_ancestor(
            args.tree_path, args.tumour_df_path
        )
        save_dataframe_to_csv(
            df=cn_change_to_ancestor_df,
            output_dir=args.output_directory,
            output_filename="cn_change_to_ancestor.csv",
        )
        logger.info(
            f"Analysis completed successfully. Output saved to: {args.output_directory}"
        )

    except Exception as e:
        logger.exception(f"An error occurred during analysis: {e}")
        exit(1)


def run_calculate_ccd():
    """CLI wrapper for calculate_ccd"""
    logger = create_logger(name="ccd_analysis", log_dir="logs")
    parser = argparse.ArgumentParser(
        description="Compute clone copy number diversity and save results."
    )
    parser.add_argument("command", choices=["ccd"], help="Command to run")
    parser.add_argument(
        "--alpaca_output_path",
        help="Path to the results dataframe file (CSV format), either the entire cohort or a single tumour",
        required=True,
    )
    parser.add_argument(
        "--output_directory", help="Path to save the output CSV file", required=True
    )

    args = parser.parse_args()
    # Validate input files exist
    if not os.path.isfile(args.alpaca_output_path):
        logger.error(f"Tumour dataframe file not found: {args.alpaca_output_path}")
        exit(1)
    # check if first row of the file contains columns 'tumour_id' and 'pred_CN_A' and 'pred_CN_B':
    with open(args.alpaca_output_path, "r") as f:
        header = f.readline().strip().split(",")
        required_columns = ["tumour_id", "clone", "segment", "pred_CN_A", "pred_CN_B"]
        missing_columns = [col for col in required_columns if col not in header]
        if missing_columns:
            logger.error(
                f"Tumour dataframe file does not contain required columns: {missing_columns}"
            )
            exit(1)
    try:
        logger.info("Starting CCD analysis...")
        ccd_scores_df = calculate_ccd(args.alpaca_output_path)

        # Ensure output directory exists
        output_dir = args.output_directory
        os.makedirs(output_dir, exist_ok=True)
        output_name = f"{output_dir}/clone_copy_number_diversity_scores.csv"
        ccd_scores_df.to_csv(output_name, index=False)
        logger.info(f"Analysis completed successfully. Output saved to: {output_name}")

    except Exception as e:
        logger.exception(f"An error occurred during analysis: {e}")
        exit(1)


def run_plot_tumour():
    """CLI wrapper for generating plots/notebooks once ALPACA outputs exist."""

    logger = create_logger(name="plotting", log_dir="logs")
    parser = argparse.ArgumentParser(
        description="Generate ALPACA heatmap plots (PDF or notebook) from existing outputs."
    )
    parser.add_argument(
        "command",
        choices=["plot-tumour"],
        help="Command to run",
    )
    parser.add_argument(
        "--input_directory",
        required=True,
        help="Tumour input directory used during the original ALPACA run.",
    )
    parser.add_argument(
        "--output_directory",
        required=True,
        help="Existing ALPACA output directory containing ALPACA_output*.csv files.",
    )
    parser.add_argument(
        "--alpaca-output-path",
        dest="alpaca_output_path",
        default=None,
        help="Optional explicit path to the ALPACA_output CSV (defaults to autodetect in the output directory).",
    )
    parser.add_argument(
        "--chr-table",
        dest="chr_table",
        default=None,
        help="Optional path to a chromosome length table (CSV/TSV).",
    )
    parser.add_argument(
        "--mutation-table",
        dest="mutation_table",
        default=None,
        help="Optional path to a driver mutation table (CSV/TSV).",
    )
    parser.add_argument(
        "--plot-output-mode",
        dest="plot_output_mode",
        choices=["pdf", "notebook", "none"],
        default="notebook",
        help="Choose 'pdf' for static figures, 'notebook' (default), or 'none' to skip generation.",
    )
    parser.add_argument(
        "--notebook-name",
        dest="notebook_name",
        default=None,
        help="Optional filename for the generated notebook (defaults to <tumour_id>_plots.ipynb).",
    )
    parser.add_argument(
        "--heatmap-palette",
        dest="heatmap_palette",
        default=_DEFAULT_HEATMAP_PALETTE,
        choices=_SUPPORTED_HEATMAP_CHOICES,
        help="Colour palette for copy-number gains (>=2 states).",
    )
    parser.add_argument(
        "--genome-build",
        dest="genome_build",
        default=_DEFAULT_GENOME_BUILD,
        choices=SUPPORTED_GENOME_BUILDS,
        help="Reference genome build for chromosome lengths (default: hg19).",
    )

    args = parser.parse_args()

    def _resolve_dir(path_str, label):
        path = Path(path_str).expanduser().resolve()
        if not path.exists() or not path.is_dir():
            logger.error(f"{label} not found or is not a directory: {path}")
            sys.exit(1)
        return path

    def _resolve_file(path_str, label):
        if path_str is None:
            return None
        path = Path(path_str).expanduser().resolve()
        if not path.exists():
            logger.error(f"{label} not found: {path}")
            sys.exit(1)
        return path

    input_dir = _resolve_dir(args.input_directory, "Input directory")
    output_dir = _resolve_dir(args.output_directory, "Output directory")
    alpaca_output_path = _resolve_file(
        args.alpaca_output_path, "ALPACA output CSV"
    )
    chr_table_path = _resolve_file(args.chr_table, "Chromosome table")
    mutation_table_path = _resolve_file(args.mutation_table, "Mutation table")

    logger.info(
        "Generating %s plot artefacts for %s",
        args.plot_output_mode,
        output_dir,
    )

    try:
        generated = export_plot_outputs(
            mode=args.plot_output_mode,
            tumour_input_dir=input_dir,
            tumour_output_dir=output_dir,
            alpaca_output_path=alpaca_output_path,
            chr_table_override=chr_table_path,
            mutation_table_override=mutation_table_path,
            notebook_name=args.notebook_name,
            heatmap_palette=args.heatmap_palette,
            genome_build=args.genome_build,
        )
    except Exception as exc:
        logger.exception(f"Failed to generate plots: {exc}")
        sys.exit(1)

    if not generated:
        logger.info(
            "Plot generation complete. No artefacts created (mode=%s)",
            args.plot_output_mode,
        )
    else:
        for artifact in generated:
            logger.info(f"Created {artifact}")
        logger.info("Plot generation complete.")
