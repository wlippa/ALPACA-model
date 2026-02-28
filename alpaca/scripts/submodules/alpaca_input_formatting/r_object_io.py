from __future__ import annotations

import inspect
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Mapping, Sequence
import warnings

import numpy as np
import pandas as pd


def _import_rdata():
    try:
        import rdata  # type: ignore
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RuntimeError(
            "Missing Python dependency 'rdata'. Install it with: pip install rdata"
        ) from exc
    return rdata


def _convert_parsed(parsed: Any, rdata_module: Any) -> Any:
    conversion = getattr(rdata_module, "conversion", None)
    if conversion is None or not hasattr(conversion, "convert"):
        return parsed

    convert_fn = conversion.convert
    kwargs: Dict[str, Any] = {}
    try:
        signature = inspect.signature(convert_fn)
    except (TypeError, ValueError):
        signature = None

    if signature is not None and "default_encoding" in signature.parameters:
        kwargs["default_encoding"] = "utf-8"

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r'Missing constructor for R class ".*"',
            category=UserWarning,
        )
        return convert_fn(parsed, **kwargs)


def _read_with_fallback(path: str | Path, read_names: list[str]) -> Any:
    rdata = _import_rdata()
    path_str = str(path)
    for fn_name in read_names:
        fn = getattr(rdata, fn_name, None)
        if callable(fn):
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r'Missing constructor for R class ".*"',
                    category=UserWarning,
                )
                return fn(path_str)

    parser = getattr(rdata, "parser", None)
    if parser is None or not hasattr(parser, "parse_file"):
        raise RuntimeError(
            "Could not find a supported rdata API for reading R files. "
            "Expected one of read_rds/read_rda or parser.parse_file."
        )

    parsed = parser.parse_file(path_str)
    return _convert_parsed(parsed, rdata)


def read_rds(path: str | Path) -> Any:
    return _read_with_fallback(path, ["read_rds", "read_rds_file"])


def read_rdata(path: str | Path) -> Any:
    obj = _read_with_fallback(
        path, ["read_rda", "read_rdata", "read_rda_file", "read_rdata_file"]
    )
    if isinstance(obj, dict):
        return obj
    # Some rdata versions may return a single object for .RData files.
    return {"object": obj}


def _to_mapping(obj: Any) -> Mapping[str, Any] | None:
    if isinstance(obj, np.void) and obj.dtype.names:
        return {str(name): obj[name] for name in obj.dtype.names}
    if isinstance(obj, Mapping):
        return obj
    if isinstance(obj, SimpleNamespace):
        return vars(obj)
    if hasattr(obj, "__dict__"):
        try:
            return vars(obj)
        except TypeError:
            return None
    return None


def _map_key(mapping: Mapping[str, Any], key: str) -> str | None:
    key_lower = key.lower()
    for existing in mapping.keys():
        if str(existing).lower() == key_lower:
            return existing
    return None


def _map_get(mapping: Mapping[str, Any], key: str, default: Any = None) -> Any:
    resolved = _map_key(mapping, key)
    if resolved is None:
        return default
    return mapping[resolved]


def _coerce_scalar(value: Any) -> Any:
    if isinstance(value, (np.generic,)):
        return value.item()
    return value


def _unwrap_singleton(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _unwrap_singleton(value.item())
        if value.size == 1:
            return _unwrap_singleton(value.reshape(-1)[0])
    if isinstance(value, (list, tuple)) and len(value) == 1:
        return _unwrap_singleton(value[0])
    return value


def _vectorize(obj: Any) -> list[Any]:
    if obj is None:
        return []
    if isinstance(obj, pd.Categorical):
        return [_coerce_scalar(v) for v in obj.tolist()]
    if isinstance(obj, pd.Series):
        return [_coerce_scalar(v) for v in obj.tolist()]
    if isinstance(obj, np.ma.MaskedArray):
        return [_coerce_scalar(v) if v is not np.ma.masked else None for v in obj.tolist()]
    if isinstance(obj, np.ndarray):
        return [_coerce_scalar(v) for v in obj.tolist()]
    if isinstance(obj, (list, tuple)):
        return [_coerce_scalar(v) for v in obj]

    mapping = _to_mapping(obj)
    if mapping is not None:
        # Factor-like representation.
        codes_key = _map_key(mapping, "codes")
        levels_key = _map_key(mapping, "levels")
        if codes_key is not None and levels_key is not None:
            codes = _vectorize(mapping[codes_key])
            levels = _vectorize(mapping[levels_key])
            out: list[Any] = []
            for code in codes:
                try:
                    idx = int(code)
                except (TypeError, ValueError):
                    out.append(None)
                    continue
                if idx <= 0 or idx > len(levels):
                    out.append(None)
                else:
                    out.append(levels[idx - 1])
            return out

        # Rle-like representation.
        values_key = _map_key(mapping, "values")
        lengths_key = _map_key(mapping, "lengths")
        if values_key is not None and lengths_key is not None:
            values = _vectorize(mapping[values_key])
            lengths = _vectorize(mapping[lengths_key])
            out: list[Any] = []
            for value, length in zip(values, lengths):
                try:
                    n = int(length)
                except (TypeError, ValueError):
                    n = 0
                if n > 0:
                    out.extend([value] * n)
            return out

        # Plain wrapper classes.
        data_key = _map_key(mapping, "data")
        if data_key is not None:
            return _vectorize(mapping[data_key])
        dot_data_key = _map_key(mapping, ".Data")
        if dot_data_key is not None:
            return _vectorize(mapping[dot_data_key])
        if len(mapping) == 1:
            only_value = next(iter(mapping.values()))
            return _vectorize(only_value)

    if isinstance(obj, (str, bytes)):
        return [obj]
    return [_coerce_scalar(obj)]


def _broadcast(values: list[Any], target_len: int) -> list[Any]:
    if target_len <= 0:
        return values
    if len(values) == target_len:
        return values
    if len(values) == 0:
        return [None] * target_len
    if len(values) == 1:
        return values * target_len
    if len(values) > target_len:
        return values[:target_len]
    # If shorter but not scalar, pad with last value to avoid shape errors.
    return values + [values[-1]] * (target_len - len(values))


def _extract_named_columns(mapping: Mapping[str, Any]) -> Dict[str, list[Any]]:
    reserved = {
        "class",
        "attributes",
        "names",
        "rownames",
        "row.names",
        "dim",
        "dimnames",
        "levels",
        "codes",
        "slotnames",
        "package",
        "version",
    }
    cols: Dict[str, list[Any]] = {}
    lengths: list[int] = []
    for key, value in mapping.items():
        key_str = str(key)
        if key_str in reserved or key_str.startswith("_"):
            continue
        vec = _vectorize(value)
        if not vec:
            continue
        cols[key_str] = vec
        lengths.append(len(vec))
    if not cols:
        return {}

    target_len = max(lengths)
    return {k: _broadcast(v, target_len) for k, v in cols.items()}


def _iter_data_entries(mapping: Mapping[str, Any]) -> Sequence[tuple[str, Any]]:
    reserved = {
        "class",
        "attributes",
        "names",
        "rownames",
        "row.names",
        "dim",
        "dimnames",
        "levels",
        "codes",
        "slotnames",
        "package",
        "version",
        "elementtype",
        "elementmetadata",
        "metadata",
    }
    entries: list[tuple[str, Any]] = []
    for key, value in mapping.items():
        key_str = str(key)
        if key_str.startswith("_"):
            continue
        if key_str.lower() in reserved:
            continue
        entries.append((key_str, value))
    return entries


def _is_granges_like(obj: Any) -> bool:
    obj = _unwrap_singleton(obj)
    mapping = _to_mapping(obj)
    if mapping is None:
        return False
    return _map_key(mapping, "ranges") is not None and _map_key(mapping, "seqnames") is not None


def _partition_to_group_names(partitioning_obj: Any, n_rows: int) -> list[Any]:
    if n_rows <= 0:
        return []

    mapping = _to_mapping(partitioning_obj)
    if mapping is None:
        return [None] * n_rows

    ends = _vectorize(_map_get(mapping, "end"))
    names = _vectorize(_map_get(mapping, "NAMES", _map_get(mapping, "names")))

    if not ends:
        return [None] * n_rows

    end_ints: list[int] = []
    for raw in ends:
        try:
            end_ints.append(int(raw))
        except (TypeError, ValueError):
            end_ints.append(0)

    if not names:
        names = [f"group_{i + 1}" for i in range(len(end_ints))]
    names = _broadcast(names, len(end_ints))

    out = [None] * n_rows
    prev = 0
    for idx, end in enumerate(end_ints):
        if end <= prev:
            continue
        # PartitioningByEnd stores 1-based inclusive end indices.
        start_ix = max(prev, 0)
        end_ix = min(end, n_rows)
        for row_ix in range(start_ix, end_ix):
            out[row_ix] = names[idx]
        prev = end_ix
    return out


def _granges_list_to_dataframe(mapping: Mapping[str, Any], *, name: str) -> pd.DataFrame:
    unlist_data = _map_get(mapping, "unlistData")
    if unlist_data is None:
        raise TypeError(f"Could not convert '{name}' as GRanges list: missing unlistData.")

    base_df = to_dataframe(unlist_data, name=f"{name}$unlistData")
    if base_df.empty:
        return base_df

    partitioning = _map_get(mapping, "partitioning")
    if partitioning is not None:
        base_df["group_name"] = _partition_to_group_names(partitioning, len(base_df))

    # Keep group_name first to mimic common refphase table layout.
    if "group_name" in base_df.columns:
        ordered = ["group_name"] + [c for c in base_df.columns if c != "group_name"]
        base_df = base_df.loc[:, ordered]
    return base_df


def _granges_sample_map_to_dataframe(mapping: Mapping[str, Any], *, name: str) -> pd.DataFrame:
    entries = _iter_data_entries(mapping)
    if not entries:
        raise TypeError(f"Could not convert '{name}' as sample->GRanges mapping.")

    if not all(_is_granges_like(value) for _, value in entries):
        raise TypeError(f"'{name}' is not a pure sample->GRanges mapping.")

    frames: list[pd.DataFrame] = []
    for sample_name, sample_obj in entries:
        sample_obj = _unwrap_singleton(sample_obj)
        frame = to_dataframe(sample_obj, name=f"{name}${sample_name}")
        if frame.empty:
            continue
        frame = frame.copy()
        frame["group_name"] = sample_name
        frames.append(frame)

    if not frames:
        return pd.DataFrame()

    out = pd.concat(frames, ignore_index=True, sort=False)
    ordered = ["group_name"] + [c for c in out.columns if c != "group_name"]
    return out.loc[:, ordered]


def _maybe_flatten_granges_sample_wide_df(df: pd.DataFrame, *, name: str) -> pd.DataFrame:
    if df.empty or len(df) != 1:
        return df

    sample_map: Dict[str, Any] = {}
    for column in df.columns:
        value = _unwrap_singleton(df.iloc[0][column])
        if not _is_granges_like(value):
            return df
        sample_map[str(column)] = value

    if not sample_map:
        return df
    return _granges_sample_map_to_dataframe(sample_map, name=name)


def _to_dataframe_from_mapping(mapping: Mapping[str, Any], *, name: str) -> pd.DataFrame:
    if _map_key(mapping, "unlistData") is not None and _map_key(mapping, "partitioning") is not None:
        return _granges_list_to_dataframe(mapping, name=name)

    data_entries = _iter_data_entries(mapping)
    if data_entries and all(_is_granges_like(value) for _, value in data_entries):
        return _granges_sample_map_to_dataframe(mapping, name=name)

    list_data_key = _map_key(mapping, "listData")
    if list_data_key is not None:
        list_data = mapping[list_data_key]
        list_mapping = _to_mapping(list_data)
        if list_mapping is not None:
            cols = _extract_named_columns(list_mapping)
            if cols:
                return pd.DataFrame(cols)
        return to_dataframe(list_data, name=f"{name}$listData")

    cols = _extract_named_columns(mapping)
    if cols:
        return pd.DataFrame(cols)

    raise TypeError(f"Could not convert '{name}' mapping to pandas DataFrame.")


def _extract_ranges(ranges_obj: Any) -> Dict[str, list[Any]]:
    mapping = _to_mapping(ranges_obj)
    if mapping is None:
        raise TypeError("Could not decode GRanges ranges slot.")

    pos = _vectorize(_map_get(mapping, "pos"))
    if pos:
        return {
            "pos": pos,
        }

    start = _vectorize(_map_get(mapping, "start"))
    width = _vectorize(_map_get(mapping, "width"))
    end = _vectorize(_map_get(mapping, "end"))

    if not end and start and width:
        n = max(len(start), len(width))
        start_b = _broadcast(start, n)
        width_b = _broadcast(width, n)
        end = []
        for s, w in zip(start_b, width_b):
            try:
                end.append(int(s) + int(w) - 1)
            except (TypeError, ValueError):
                end.append(None)
        return {"start": start_b, "width": width_b, "end": end}

    n = max(len(start), len(width), len(end))
    if n == 0:
        return {"start": [], "width": [], "end": []}
    start_b = _broadcast(start, n)
    width_b = _broadcast(width, n)
    end_b = _broadcast(end, n)
    if not width and start and end:
        width_b = []
        for s, e in zip(start_b, end_b):
            try:
                width_b.append(int(e) - int(s) + 1)
            except (TypeError, ValueError):
                width_b.append(None)
    return {"start": start_b, "width": width_b, "end": end_b}


def _granges_to_dataframe(obj: Any, *, name: str) -> pd.DataFrame:
    mapping = _to_mapping(obj)
    if mapping is None:
        raise TypeError(f"Could not convert '{name}' GRanges object.")

    seqnames = _vectorize(_map_get(mapping, "seqnames"))
    strand = _vectorize(_map_get(mapping, "strand"))
    ranges = _extract_ranges(_map_get(mapping, "ranges"))
    start = ranges.get("start", [])
    width = ranges.get("width", [])
    end = ranges.get("end", [])
    pos = ranges.get("pos", [])

    n = max(len(seqnames), len(strand), len(start), len(width), len(end), len(pos))
    if n == 0:
        return pd.DataFrame(columns=["seqnames", "start", "end", "width", "strand"])

    data = {
        "seqnames": _broadcast(seqnames, n),
        "strand": _broadcast(strand, n),
    }
    if start:
        data["start"] = _broadcast(start, n)
    if end:
        data["end"] = _broadcast(end, n)
    if width:
        data["width"] = _broadcast(width, n)
    if pos:
        data["pos"] = _broadcast(pos, n)
    frame = pd.DataFrame(data)

    meta = _map_get(mapping, "elementMetadata")
    if meta is None:
        meta = _map_get(mapping, "mcols")
    if meta is not None:
        try:
            meta_df = to_dataframe(meta, name=f"{name}$elementMetadata")
            if len(meta_df) == n:
                frame = pd.concat(
                    [frame.reset_index(drop=True), meta_df.reset_index(drop=True)],
                    axis=1,
                )
        except Exception:
            pass
    return frame


def get_field(obj: Any, field: str) -> Any:
    mapping = _to_mapping(obj)
    if mapping is not None:
        resolved = _map_key(mapping, field)
        if resolved is not None:
            return mapping[resolved]
        available = sorted(str(k) for k in mapping.keys())
        raise KeyError(f"Missing field '{field}'. Available keys: {available}")
    if hasattr(obj, field):
        return getattr(obj, field)
    raise KeyError(
        f"Missing field '{field}' in object of type {type(obj).__name__}. "
        "Expected a dict-like object or class with attributes."
    )


def to_dataframe(obj: Any, *, name: str) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return _maybe_flatten_granges_sample_wide_df(obj.copy(), name=name)

    if isinstance(obj, pd.Series):
        return obj.to_frame()

    if isinstance(obj, np.ndarray):
        if obj.dtype.names and obj.size > 0:
            first = obj.reshape(-1)[0]
            record = _unwrap_singleton(first)
            record_mapping = _to_mapping(record)
            if isinstance(record_mapping, Mapping):
                candidate_map = {
                    str(k): _unwrap_singleton(v) for k, v in record_mapping.items()
                }
                if candidate_map and all(
                    _is_granges_like(v) for v in candidate_map.values()
                ):
                    return _granges_sample_map_to_dataframe(candidate_map, name=name)
        if obj.dtype.names:
            return pd.DataFrame.from_records(obj)
        if obj.ndim == 1:
            return pd.DataFrame(obj)
        if obj.ndim == 2:
            return pd.DataFrame(obj)

    mapping = _to_mapping(obj)
    if mapping is not None:
        # GRanges-like objects returned as unconstructed S4 payloads.
        if (
            _map_key(mapping, "seqnames") is not None
            and _map_key(mapping, "ranges") is not None
            and _map_key(mapping, "strand") is not None
        ):
            return _granges_to_dataframe(obj, name=name)
        return _to_dataframe_from_mapping(mapping, name=name)

    if isinstance(obj, dict):
        return pd.DataFrame(obj)

    if hasattr(obj, "to_pandas"):
        converted = obj.to_pandas()
        if isinstance(converted, pd.Series):
            return converted.to_frame()
        if isinstance(converted, pd.DataFrame):
            return converted

    try:
        return pd.DataFrame(obj)
    except Exception as exc:  # pragma: no cover - error path
        raise TypeError(
            f"Could not convert '{name}' ({type(obj).__name__}) to pandas DataFrame."
        ) from exc


def normalize_cluster_id(value: Any) -> Any:
    if pd.isna(value):
        return value
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return value
