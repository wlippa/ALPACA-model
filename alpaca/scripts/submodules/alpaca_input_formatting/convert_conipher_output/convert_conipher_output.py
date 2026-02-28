#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from typing import Any, List

import numpy as np
import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
INPUT_FORMATTING_DIR = THIS_DIR.parent
if str(INPUT_FORMATTING_DIR) not in sys.path:
    sys.path.insert(0, str(INPUT_FORMATTING_DIR))

from r_object_io import get_field, normalize_cluster_id, read_rds, to_dataframe


def _to_mapping(obj: Any):
    if isinstance(obj, dict):
        return obj
    if isinstance(obj, SimpleNamespace):
        return vars(obj)
    if hasattr(obj, "__dict__"):
        try:
            return vars(obj)
        except TypeError:
            return None
    return None


def _as_vector(obj: Any) -> list[Any]:
    if obj is None:
        return []
    if isinstance(obj, np.ndarray):
        return obj.reshape(-1).tolist()
    if isinstance(obj, pd.Series):
        return obj.tolist()
    if isinstance(obj, (list, tuple)):
        return list(obj)
    mapping = _to_mapping(obj)
    if mapping is not None:
        if "data" in mapping:
            return _as_vector(mapping["data"])
        if ".Data" in mapping:
            return _as_vector(mapping[".Data"])
        if "values" in mapping and "lengths" in mapping:
            out = []
            values = _as_vector(mapping["values"])
            lengths = _as_vector(mapping["lengths"])
            for value, length in zip(values, lengths):
                try:
                    n = int(length)
                except (TypeError, ValueError):
                    n = 0
                if n > 0:
                    out.extend([value] * n)
            return out
    return [obj]


def _edges_from_mapping(mapping: dict) -> pd.DataFrame | None:
    # Direct parent/child columns
    lower_to_key = {str(k).lower(): k for k in mapping.keys()}
    if "parent" in lower_to_key and "child" in lower_to_key:
        parent = _as_vector(mapping[lower_to_key["parent"]])
        child = _as_vector(mapping[lower_to_key["child"]])
        n = min(len(parent), len(child))
        return pd.DataFrame({"Parent": parent[:n], "Child": child[:n]})

    # listData-style container (R data.frame / S4 slots)
    list_data = mapping.get("listData")
    if list_data is not None:
        list_mapping = _to_mapping(list_data)
        if isinstance(list_mapping, dict):
            list_lower = {str(k).lower(): k for k in list_mapping.keys()}
            if "parent" in list_lower and "child" in list_lower:
                parent = _as_vector(list_mapping[list_lower["parent"]])
                child = _as_vector(list_mapping[list_lower["child"]])
                n = min(len(parent), len(child))
                return pd.DataFrame({"Parent": parent[:n], "Child": child[:n]})

    # matrix-style flattened .Data with dim [n,2]
    raw = mapping.get(".Data", mapping.get("data"))
    dim = mapping.get("dim")
    if raw is not None and dim is not None:
        values = _as_vector(raw)
        dims = _as_vector(dim)
        if len(dims) >= 2:
            try:
                nrow = int(dims[0])
                ncol = int(dims[1])
            except (TypeError, ValueError):
                nrow = 0
                ncol = 0
            if nrow > 0 and ncol >= 2 and len(values) >= nrow * ncol:
                # R matrices are column-major: first nrow values are first column.
                parent = values[:nrow]
                child = values[nrow : 2 * nrow]
                return pd.DataFrame({"Parent": parent, "Child": child})

    return None


def _coerce_edges_dataframe(tree_graph: Any) -> pd.DataFrame:
    # 0) raw ndarray handling before generic DataFrame conversion
    if isinstance(tree_graph, np.ndarray):
        arr = tree_graph
        if arr.ndim == 2 and arr.shape[1] >= 2:
            return pd.DataFrame({"Parent": arr[:, 0], "Child": arr[:, 1]})
        if arr.ndim == 2 and arr.shape[1] == 1:
            flat = arr[:, 0]
            # Each row stores a pair-like object [parent, child].
            if all(
                isinstance(x, (list, tuple, np.ndarray, pd.Series)) and len(x) >= 2
                for x in flat
            ):
                return pd.DataFrame(
                    {
                        "Parent": [x[0] for x in flat],
                        "Child": [x[1] for x in flat],
                    }
                )
            # Column-major stacked representation: first half Parent, second half Child.
            if len(flat) % 2 == 0 and len(flat) > 0:
                half = len(flat) // 2
                return pd.DataFrame({"Parent": flat[:half], "Child": flat[half:]})
        if arr.ndim == 1 and len(arr) % 2 == 0 and len(arr) > 0:
            half = len(arr) // 2
            return pd.DataFrame({"Parent": arr[:half], "Child": arr[half:]})

    # 1) standard route
    try:
        df = to_dataframe(tree_graph, name="tree_graph")
    except Exception:
        df = pd.DataFrame()
    if not df.empty and df.shape[1] >= 2:
        out = df.iloc[:, :2].copy()
        out.columns = ["Parent", "Child"]
        return out

    # 2) one-column of pair-like rows
    if not df.empty and df.shape[1] == 1 and len(df) > 0:
        col = df.iloc[:, 0]
        if col.apply(
            lambda x: isinstance(x, (list, tuple, np.ndarray, pd.Series)) and len(x) >= 2
        ).all():
            return pd.DataFrame(
                {"Parent": col.apply(lambda x: x[0]), "Child": col.apply(lambda x: x[1])}
            )
        # Stacked single-column representation: first half Parent, second half Child.
        if len(col) % 2 == 0:
            half = len(col) // 2
            return pd.DataFrame(
                {
                    "Parent": col.iloc[:half].to_list(),
                    "Child": col.iloc[half:].to_list(),
                }
            )

    # 3) mapping-based fallbacks for S4/matrix/listData representations
    mapping = _to_mapping(tree_graph)
    if mapping is not None:
        parsed = _edges_from_mapping(dict(mapping))
        if parsed is not None and parsed.shape[1] >= 2:
            return parsed

    mapping = _to_mapping(tree_graph)
    keys = sorted(mapping.keys()) if mapping is not None else "n/a"
    raise ValueError(
        "Tree graph must have at least two columns (Parent, Child). "
        f"Decoded shape was {df.shape} from type {type(tree_graph).__name__}. "
        f"Available keys: {keys}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract clone proportions and tree from CONIPHER output"
    )
    parser.add_argument(
        "--CONIPHER_tree_object",
        type=str,
        required=True,
        help="Path to CONIPHER tree object (.RDS)",
    )
    parser.add_argument(
        "--CONIPHER_tree_index",
        type=int,
        default=1,
        help="Selected CONIPHER tree index (1-based)",
    )
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")
    return parser.parse_args()


def _edges_dataframe(tree_graph: Any) -> pd.DataFrame:
    df = _coerce_edges_dataframe(tree_graph)
    df["Parent"] = df["Parent"].map(normalize_cluster_id)
    df["Child"] = df["Child"].map(normalize_cluster_id)
    return df


def _normalize_index(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.index = [normalize_cluster_id(v) for v in out.index]
    return out


def _unique_cluster_values(edges: pd.DataFrame) -> List[Any]:
    values = pd.unique(edges[["Parent", "Child"]].to_numpy().ravel())
    return [normalize_cluster_id(v) for v in values if pd.notna(v)]


def infer_trunk_from_tree(tree_graph: Any) -> Any:
    edges = _edges_dataframe(tree_graph)
    parents = list(pd.unique(edges["Parent"]))
    children = set(pd.unique(edges["Child"]))
    trunk_candidates = [p for p in parents if p not in children]
    if not trunk_candidates:
        raise ValueError(
            "Could not infer trunk from selected tree graph (no root candidate found)."
        )
    return normalize_cluster_id(trunk_candidates[0])


def get_tree_level(tree_graph: pd.DataFrame, cluster: Any) -> float:
    unique_values = _unique_cluster_values(tree_graph)
    if len(set(unique_values)) == 1:
        return 1.0

    parents = list(pd.unique(tree_graph["Parent"]))
    children = set(pd.unique(tree_graph["Child"]))
    trunk_candidates = [p for p in parents if p not in children]
    if not trunk_candidates:
        return float("nan")
    trunk = trunk_candidates[0]
    cluster = normalize_cluster_id(cluster)

    if cluster == trunk:
        return 1.0

    clusters_in_tree = set(unique_values)
    if cluster not in clusters_in_tree:
        return float("nan")

    level = 1.0
    current_cluster = cluster
    while current_cluster != trunk:
        parent = tree_graph.loc[tree_graph["Child"] == current_cluster, "Parent"]
        if parent.empty:
            return float("nan")
        current_cluster = normalize_cluster_id(parent.iloc[0])
        level += 1.0
    return level


def _as_tree_list(alt_trees: Any) -> List[Any]:
    if isinstance(alt_trees, dict):
        keys = list(alt_trees.keys())
        try:
            ordered_keys = sorted(keys, key=lambda k: int(float(k)))
        except (TypeError, ValueError):
            ordered_keys = keys
        return [alt_trees[key] for key in ordered_keys]
    if isinstance(alt_trees, (list, tuple)):
        return list(alt_trees)
    raise TypeError(
        f"Unexpected type for graph_pyclone$alt_trees: {type(alt_trees).__name__}"
    )


def compute_subclone_proportions(
    tree_list: List[Any],
    ccf_cluster_table: pd.DataFrame,
    clonality_table: pd.DataFrame,
    trunk: Any,
    force_clonal_100: bool = True,
    tree_id: int = 1,
) -> pd.DataFrame:
    tree_edges = _edges_dataframe(tree_list[tree_id - 1])
    clusters_in_tree = set(_unique_cluster_values(tree_edges))

    ccf = _normalize_index(ccf_cluster_table)
    ccf = ccf.loc[ccf.index.isin(clusters_in_tree)].copy()
    ccf = ccf.apply(pd.to_numeric, errors="coerce")

    if len(clusters_in_tree) == 1:
        return (ccf > 0).astype(float) * 100.0

    region_ids = list(ccf.columns)
    clonality = _normalize_index(clonality_table)
    clonality = clonality.loc[clonality.index.isin(clusters_in_tree)].copy()

    if force_clonal_100:
        clonal_factor = (clonality == "clonal").astype(float)
        if clonal_factor.shape[1] == 1 and len(region_ids) > 1:
            one_col = clonal_factor.iloc[:, 0]
            clonal_factor = pd.concat([one_col] * len(region_ids), axis=1)
            clonal_factor.columns = region_ids
        else:
            clonal_factor = clonal_factor.reindex(columns=region_ids).fillna(0.0)
        clonal_factor = clonal_factor.reindex(index=ccf.index).fillna(0.0)
        ccf = ccf * (1.0 - clonal_factor) + 100.0 * clonal_factor

    ccf = ccf.clip(upper=100.0)
    trunk = normalize_cluster_id(trunk)
    if trunk in ccf.index:
        ccf.loc[trunk, :] = 100.0

    ccf_cluster_df = ccf.copy()
    ccf_cluster_df["cluster"] = ccf_cluster_df.index

    proportions_df = ccf.copy()
    proportions_df["cluster"] = proportions_df.index
    proportions_df.loc[:, region_ids] = 0.0

    for region_id in region_ids:
        region_ccf = ccf_cluster_df[[region_id, "cluster"]].copy()
        region_ccf.columns = ["ccf", "cluster"]
        clusters_present = set(region_ccf.loc[region_ccf["ccf"] != 0, "cluster"])
        parents_present = [
            p for p in pd.unique(tree_edges["Parent"]) if p in clusters_present
        ]
        parent_df = pd.DataFrame({"parent_node": parents_present})
        parent_df["level"] = parent_df["parent_node"].apply(
            lambda p: get_tree_level(tree_edges, p)
        )
        parent_df = parent_df.sort_values("level")

        for parent in parent_df["parent_node"]:
            children_nodes = tree_edges.loc[tree_edges["Parent"] == parent, "Child"]
            children_nodes = list(children_nodes)
            parent_row = region_ccf.loc[region_ccf["cluster"] == parent, "ccf"]
            if parent_row.empty:
                continue
            parent_ccf = float(parent_row.iloc[0])
            child_mask = region_ccf["cluster"].isin(children_nodes)
            sum_children_ccf = float(region_ccf.loc[child_mask, "ccf"].sum())

            if sum_children_ccf > parent_ccf and sum_children_ccf > 0:
                parent_proportion = 0.0
                region_ccf.loc[child_mask, "ccf"] = (
                    parent_ccf
                    * region_ccf.loc[child_mask, "ccf"]
                    / region_ccf.loc[child_mask, "ccf"].sum()
                )
            else:
                parent_proportion = parent_ccf - sum_children_ccf

            proportions_df.loc[proportions_df["cluster"] == parent, region_id] = (
                parent_proportion
            )

            for child in children_nodes:
                is_terminal = (child in set(tree_edges["Child"])) and (
                    child not in set(tree_edges["Parent"])
                )
                if is_terminal:
                    child_value = region_ccf.loc[region_ccf["cluster"] == child, "ccf"]
                    if not child_value.empty:
                        proportions_df.loc[
                            proportions_df["cluster"] == child, region_id
                        ] = float(child_value.iloc[0])

    proportions_df.index = proportions_df["cluster"]
    return proportions_df.loc[:, region_ids].copy()


def get_cp_table(
    alt_trees: List[Any],
    alt_tree_id: int,
    clonality_table: pd.DataFrame,
    ccf_cluster_table: pd.DataFrame,
    trunk: Any,
) -> pd.DataFrame:
    cp_table = compute_subclone_proportions(
        tree_list=alt_trees,
        ccf_cluster_table=ccf_cluster_table,
        clonality_table=clonality_table,
        trunk=trunk,
        force_clonal_100=True,
        tree_id=alt_tree_id,
    )
    cp_table = (cp_table / 100.0).copy()
    cp_table["clone"] = [f"clone{idx}" for idx in cp_table.index]
    cp_table.columns = [str(col).replace(".", "-") for col in cp_table.columns]
    return cp_table


def extract_tree_graph_paths(tree_graph: Any) -> List[List[Any]]:
    edges = _edges_dataframe(tree_graph)
    clones_in_tree = _unique_cluster_values(edges)
    if len(clones_in_tree) == 1:
        return [[clones_in_tree[0]]]

    parents = list(pd.unique(edges["Parent"]))
    children = list(pd.unique(edges["Child"]))
    child_set = set(children)
    parent_set = set(parents)

    trunk_candidates = [p for p in parents if p not in child_set]
    if not trunk_candidates:
        raise ValueError("Could not identify trunk node in tree graph.")
    trunk = trunk_candidates[0]

    terminal_clones = [c for c in children if c not in parent_set]
    paths: List[List[Any]] = []
    for terminal in terminal_clones:
        path = [terminal]
        current_clone = terminal
        while current_clone != trunk:
            parent = edges.loc[edges["Child"] == current_clone, "Parent"]
            if parent.empty:
                raise ValueError(
                    f"Could not find parent for clone '{current_clone}' while tracing tree."
                )
            current_clone = parent.iloc[0]
            path.append(current_clone)
        paths.append(list(reversed(path)))
    return paths


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Starting conversion of CONIPHER output...")
    tree_object = read_rds(args.CONIPHER_tree_object)
    selected_tree_index = int(args.CONIPHER_tree_index)

    graph_pyclone = get_field(tree_object, "graph_pyclone")
    alt_trees = _as_tree_list(get_field(graph_pyclone, "alt_trees"))
    number_of_trees = len(alt_trees)
    if selected_tree_index < 1 or selected_tree_index > number_of_trees:
        raise ValueError(
            f"Selected tree index {selected_tree_index} is out of bounds. "
            f"There are {number_of_trees} trees available."
        )

    selected_tree = alt_trees[selected_tree_index - 1]
    tree_paths = extract_tree_graph_paths(selected_tree)
    tree_path_clone_names = [
        [f"clone{normalize_cluster_id(node)}" for node in path] for path in tree_paths
    ]
    with (output_dir / "tree_paths.json").open("w", encoding="utf-8") as handle:
        json.dump(tree_path_clone_names, handle)

    clonality_out = get_field(tree_object, "clonality_out")
    nested_pyclone = get_field(tree_object, "nested_pyclone")
    clonality_table = to_dataframe(
        get_field(clonality_out, "clonality_table_corrected"),
        name="clonality_out$clonality_table_corrected",
    )
    ccf_cluster_table = to_dataframe(
        get_field(nested_pyclone, "ccf_cluster_table"),
        name="nested_pyclone$ccf_cluster_table",
    )
    try:
        trunk = get_field(graph_pyclone, "trunk")
    except KeyError:
        trunk = infer_trunk_from_tree(selected_tree)

    cp_table = get_cp_table(
        alt_trees, selected_tree_index, clonality_table, ccf_cluster_table, trunk
    )
    cp_table.to_csv(output_dir / "cp_table.csv", index=False)
    print("Done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
