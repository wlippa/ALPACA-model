try:
    from importlib.metadata import version  # type: ignore
except Exception:
    try:
        from importlib_metadata import version  # type: ignore
    except Exception:
        def version(pkg_name="alpaca"):
            return "unknown"

__all__ = ["version"]
