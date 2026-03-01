    #!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd


def _read_sample_set(path: Path) -> Set[str]:
    df = pd.read_csv(path)
    if "sample" not in df.columns:
        raise ValueError(f"File {path} is missing required 'sample' column.")
    return set(df["sample"].dropna().astype(str).unique())


def _subset_sample_table(path: Path, keep_samples: Set[str]) -> int:
    df = pd.read_csv(path)
    if "sample" not in df.columns:
        raise ValueError(f"File {path} is missing required 'sample' column.")
    df = df[df["sample"].astype(str).isin(keep_samples)].copy()
    df.to_csv(path, index=False)
    return len(df)


def _read_cp_samples(path: Path) -> List[str]:
    cp_table = pd.read_csv(path, index_col="clone")
    return [str(col) for col in cp_table.columns]


def _subset_cp_table(path: Path, keep_samples: Set[str]) -> int:
    cp_table = pd.read_csv(path, index_col="clone")
    keep_cols = [col for col in cp_table.columns if str(col) in keep_samples]
    cp_table = cp_table.loc[:, keep_cols].copy()
    cp_table.to_csv(path)
    return len(keep_cols)


def _collect_segment_files(segments_dir: Path) -> List[Path]:
    if not segments_dir.exists() or not segments_dir.is_dir():
        return []
    return sorted([p for p in segments_dir.iterdir() if p.is_file() and p.suffix == ".csv"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Ensure final ALPACA conversion outputs have a consistent sample set."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory containing ALPACA_input_table.csv, ci_table.csv, cp_table.csv and optional segments/.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).expanduser().resolve()

    alpaca_input_path = output_dir / "ALPACA_input_table.csv"
    ci_table_path = output_dir / "ci_table.csv"
    cp_table_path = output_dir / "cp_table.csv"
    segments_dir = output_dir / "segments"
    report_path = output_dir / "sample_filter_report.txt"

    for required in [alpaca_input_path, ci_table_path, cp_table_path]:
        if not required.exists():
            raise FileNotFoundError(f"Required output file does not exist: {required}")

    segment_files = _collect_segment_files(segments_dir)

    sample_sets: Dict[str, Set[str]] = {
        str(alpaca_input_path): _read_sample_set(alpaca_input_path),
        str(ci_table_path): _read_sample_set(ci_table_path),
        str(cp_table_path): set(_read_cp_samples(cp_table_path)),
    }
    for seg_path in segment_files:
        sample_sets[str(seg_path)] = _read_sample_set(seg_path)

    all_samples = sorted(set().union(*sample_sets.values())) if sample_sets else []
    common_samples = set.intersection(*sample_sets.values()) if sample_sets else set()

    if not common_samples:
        raise ValueError(
            "No shared samples remain across output files after intersection. "
            f"Observed sample sets: { {k: sorted(v) for k, v in sample_sets.items()} }"
        )

    removed_samples = sorted(set(all_samples) - common_samples)
    kept_samples_sorted = sorted(common_samples)

    _subset_sample_table(alpaca_input_path, common_samples)
    _subset_sample_table(ci_table_path, common_samples)
    _subset_cp_table(cp_table_path, common_samples)
    for seg_path in segment_files:
        _subset_sample_table(seg_path, common_samples)

    report_lines = [
        "Sample consistency check",
        "========================",
        f"files_checked: {len(sample_sets)}",
        f"kept_samples_count: {len(kept_samples_sorted)}",
        f"kept_samples: {', '.join(kept_samples_sorted)}",
        f"removed_samples_count: {len(removed_samples)}",
        (
            "removed_samples: none"
            if not removed_samples
            else f"removed_samples: {', '.join(removed_samples)}"
        ),
    ]
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"Wrote sample consistency report to {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
