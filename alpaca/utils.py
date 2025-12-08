import os
import pandas as pd
from pathlib import Path
import json
import importlib
import logging
from datetime import datetime
from typing import Optional
import glob
import sys
import csv


def show_version():
    try:
        version = importlib.metadata.version("alpaca")
        print(f"alpaca {version}")
    except importlib.metadata.PackageNotFoundError:
        print("alpaca version unknown (not installed)")


def show_help():
    print(
        "ALPACA = ALlele-specific Phylogenetic Analysis of clone Copy-number Alterations"
    )
    print_logo()
    print("")
    print("Usage:")
    print("  alpaca [command]")
    print("")
    print("Commands:")
    print("  version              Show version")
    print("  help                 Show this help")
    print("  run                  Run ALPACA")
    print("  input-conversion     Run input conversion")
    print("  ccd                  Calculate clone copy number diversity")
    print("")


def create_logger(name: str, log_dir: Optional[str] = None) -> logging.Logger:
    """
    Create a named logger with both console and file handlers.

    Args:
        name: Name for the logger
        log_dir: Optional directory for log files (defaults to current directory)

    Returns:
        Configured logger instance
    """
    # create logger
    logger = logging.getLogger(name)
    # check for active handlers
    if logger.handlers:
        return logger

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler
    log_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{name}_log_{log_time}.log"

    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, log_filename)
    else:
        log_path = log_filename

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info(f"Logger '{name}' initialized. Log file: {log_path}")

    return logger


def split_to_segments(tumour_dir: str) -> list[str]:
    segments_dir_path = f"{tumour_dir}/segments"
    df_path = f"{tumour_dir}/ALPACA_input_table.csv"
    os.makedirs(segments_dir_path, exist_ok=True)
    df = pd.read_csv(df_path)
    tumour_id = df["tumour_id"].iloc[0]
    assert len(
        df["tumour_id"].unique()
    ), "Found multiple tumour ids. In tummour mode only one tumour_id is allowed per input csv"
    segments = []
    for segment, segment_df in df.groupby("segment"):
        segment_df_path = (
            f"{segments_dir_path}/ALPACA_input_table_{tumour_id}_{segment}.csv"
        )
        segment_df.to_csv(segment_df_path, index=False)
        segments.append(segment_df_path)
    return segments


def concatenate_output(output_dir: str) -> str:
    logger = logging.getLogger("ALPACA")
    # keep only segment files in output files list
    output_files = [
        f
        for f in os.listdir(output_dir)
        if f.endswith(".csv") and (("optimal" in f) or ("all" in f))
    ]
    dfs = [pd.read_csv(f"{output_dir}/{f}") for f in output_files]
    concatenated_df = pd.concat(dfs)
    tumour_id = concatenated_df["tumour_id"].iloc[0]
    output_name = f"{output_dir}/ALPACA_output_{tumour_id}.csv"
    concatenated_df.to_csv(output_name, index=False)
    if os.path.exists(output_name):
        logger.info(f"Combined output saved to {output_name}")
        # remove segment files
        for f in output_files:
            if f != os.path.basename(output_name):
                os.remove(f"{output_dir}/{f}")
    else:
        logger.error(f"Failed to save combined output to {output_name}")
        raise FileNotFoundError(f"Output file not found: {output_name}")
    return output_name


def set_run_mode(config: dict) -> tuple[dict, str]:
    run_mode = config["preprocessing_config"]["mode"]
    if run_mode == "tumour":
        print("Running in tumour mode")
        # create segment files:
        config["preprocessing_config"]["input_files"] = [
            Path(x).name
            for x in split_to_segments(
                config["preprocessing_config"]["input_tumour_directory"]
            )
        ]
        config["preprocessing_config"]["input_data_directory"] = Path(
            config["preprocessing_config"]["input_tumour_directory"]
        ).parent
    return config, run_mode


def read_tree_json(json_path: str) -> list[list[str]]:
    with open(json_path, "r") as f:
        tree = json.load(f)
    return tree


def find_path_edges(branch, tree_edges):
    branch_edges = []
    for edge in tree_edges:
        if (edge[0] in branch) and (edge[1] in branch):
            branch_edges.append(edge)
    return set(branch_edges)


def get_tree_edges(tree):
    all_edges = list()
    for branch in tree:
        if len(branch) == 2:
            all_edges.append(tuple(branch))
        else:
            for i in range(len(branch) - 1):
                all_edges.append((branch[i], branch[i + 1]))
    unique_edges = set(all_edges)
    return unique_edges


def find_parent(tree, clone_name):
    for branch in tree:
        if branch[0] == clone_name:
            return "diploid"
        if clone_name in branch:
            clone_index = branch.index(clone_name)
            return branch[clone_index - 1]


def flat_list(target_list):
    if isinstance(target_list[0], list):
        return [item for sublist in target_list for item in sublist]
    else:
        return target_list


def get_length_from_name(segment):
    e = int(segment.split("_")[-1])
    s = int(segment.split("_")[-2])
    return e - s


def print_logo():
    print(
        """
     _____ __    _____ _____ _____ _____
    |  _  |  |  |  _  |  _  |     |  _  |
    |     |  |__|   __|     |   --|     |
    |__|__|_____|__|  |__|__|_____|__|__|
    /\\⌒⌒⌒/\\
    (⦿   ⦿)
    ( 'Y' )
     (   )
     (   )
     (   )
     (~ ~~~~~~~~~~)
     ( ~ ~~   ~~  )
     ( ~  ~ ~  ~  )
     (~  ~~~~~   ~)
     │ │      │ │
     │ │      │ │
    """
    )


def save_dataframe_to_csv(df: pd.DataFrame, output_dir: str, output_filename: str):
    """
    Save a DataFrame to a CSV file.

    Args:
        df: DataFrame to save
        output_dir: Path to directory where the CSV file will be saved
        output_filename: Name of the output CSV file
    """
    output_path = os.path.join(output_dir, output_filename)
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    df.to_csv(output_path, index=False)


def process_ci_reports(dirpath: str, delete: bool = False, outpath: Optional[str] = None) -> pd.DataFrame:
    pattern = os.path.join(dirpath, "*_ci_report.json")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        outpath = os.path.join(dirpath, "ci_modified_report.csv")

    header = [
        "tumour_id",
        "segment",
        "affected_sample",
        "affected_allele",
        "min_ci",
        "timestamp",
        "source_file",
    ]

    # write header even if no files
    with open(outpath, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        writer.writerow(header)
        for fp in files:
            try:
                with open(fp) as fh:
                    data = json.load(fh)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue

            tumour = data.get("tumour_id")
            segment = data.get("segment")
            min_ci = data.get("min_ci")
            ts = data.get("timestamp")
            affected_samples = data.get("affected_samples") or []
            affected_alleles = data.get("affected_alleles") or []

            # produce a row per sample x allele combination
            if not affected_samples and not affected_alleles:
                writer.writerow(
                    [tumour, segment, "", "", min_ci, ts, os.path.basename(fp)]
                )
                continue

            for s in affected_samples:
                if affected_alleles:
                    for a in affected_alleles:
                        writer.writerow(
                            [tumour, segment, s, a, min_ci, ts, os.path.basename(fp)]
                        )
                else:
                    writer.writerow(
                        [tumour, segment, s, "", min_ci, ts, os.path.basename(fp)]
                    )

            if delete and bool(files):
                try:
                    os.remove(fp)
                except Exception as e:
                    print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)


def process_monoclonal_reports(dirpath: str, delete: bool = False, outpath: Optional[str] = None) -> pd.DataFrame:
    pattern = os.path.join(dirpath, "*_monoclonal_samples_report.csv")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        combined_df = os.path.join(dirpath, "monoclonal_samples_report.csv")
    # each file is a csv with header, concatenate them or make an empty dataframe
    if not files:
        combined_df = pd.DataFrame(columns=['tumour_id', 'sample', 'segment', 'cpnA', 'cpnB', 'distance_to_integer_A', 'distance_to_integer_B'])
    else:
        dfs = []
        for fp in files:
            try:
                df = pd.read_csv(fp)
                dfs.append(df)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv(outpath, index=False)
    if delete and bool(files):
        try:
            os.remove(fp)
        except Exception as e:
            print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)


def process_elbow_increase_reports(dirpath: str, delete: bool = False, outpath: Optional[str] = None) -> pd.DataFrame:
    pattern = os.path.join(dirpath, "*_elbow_increase_report.csv")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        combined_df = os.path.join(dirpath, "elbow_increase_report.csv")
    # each file is a csv with header, concatenate them or make an empty dataframe
    if not files:
        combined_df = pd.DataFrame(columns=['complexity', 'D_score', 'CI_score', 'allowed_complexity', 'issue', 'tumour_id', 'segment'])
    else:
        dfs = []
        for fp in files:
            try:
                df = pd.read_csv(fp)
                dfs.append(df)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv(outpath, index=False)
    if delete and bool(files):
        try:
            os.remove(fp)
        except Exception as e:
            print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)


def process_run_summary_reports(dirpath: str, delete: bool = False, outpath: Optional[str] = None) -> pd.DataFrame:
    """Combine segment-level run_summary reports into a single file.
    
    Only includes segments with non-zero gap (non-optimal solutions).
    Reports gap value and reason (time_limit, gap_tolerance, or other).
    """
    pattern = os.path.join(dirpath, "*_run_summary.csv")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        outpath = os.path.join(dirpath, "run_summary.csv")
    
    # Concatenate all segment reports
    if not files:
        combined_df = pd.DataFrame(columns=[
            'tumour_id', 'segment', 'max_gap', 'gap_reason', 
            'runtime', 'optimal_complexity', 'strict_gap_enabled'
        ])
    else:
        dfs = []
        for fp in files:
            try:
                df = pd.read_csv(fp)
                dfs.append(df)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
        else:
            combined_df = pd.DataFrame(columns=[
                'tumour_id', 'segment', 'max_gap', 'gap_reason', 
                'runtime', 'optimal_complexity', 'strict_gap_enabled'
            ])
    
    combined_df.to_csv(outpath, index=False)
    
    if delete and files:
        for fp in files:
            try:
                os.remove(fp)
            except Exception as e:
                print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)