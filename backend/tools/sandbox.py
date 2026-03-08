import logging

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


def _serialize(v):
    """Convert any value to a JSON-serializable form."""
    if v is None:
        return None
    if isinstance(v, (bool, int, float, str)):
        return v
    if isinstance(v, (list, tuple)):
        return [_serialize(i) for i in v]
    if isinstance(v, dict):
        return {str(k): _serialize(val) for k, val in v.items()}
    if isinstance(v, pd.DataFrame):
        if len(v) > 100:
            v = v.head(100)
        return v.reset_index().to_dict(orient="records")
    if isinstance(v, pd.Series):
        if len(v) > 100:
            v = v.head(100)
        return v.reset_index().to_dict(orient="records")
    if isinstance(v, np.ndarray):
        return v.tolist()
    if hasattr(v, "tolist"):
        return v.tolist()
    if hasattr(v, "item"):
        return v.item()
    return str(v)


# Keys that are pre-loaded into the exec namespace — not part of agent output
_BUILTIN_KEYS = frozenset({"datasets", "np", "pd", "stats", "__builtins__"})


def execute_sandbox(code: str, datasets: list) -> dict:
    """
    Executes Python code written by the agent.
    Available variables: datasets, np, pd, stats
    Code must set the `result` variable before finishing.

    Uses a single namespace dict for both globals and locals so that
    variable assignments inside helper functions are visible at the top level.
    """
    namespace = {
        "datasets": datasets,
        "np": np,
        "pd": pd,
        "stats": stats,
    }
    try:
        exec(compile(code, "<agent_code>", "exec"), namespace)  # noqa: S102

        result = namespace.get("result")

        # Fallback: if result not set, collect all non-builtin variables the agent created
        if result is None:
            candidates = {
                k: v for k, v in namespace.items()
                if k not in _BUILTIN_KEYS and not k.startswith("_")
            }
            if candidates:
                serialized = {k: _serialize(v) for k, v in candidates.items()}
                result = {k: v for k, v in serialized.items() if v is not None}

        if not result:
            logger.warning("Sandbox: code did not set result variable.\nCode:\n%s", code[:500])
            return {"error": "Code did not set the `result` variable. Add: result = {'key': value} at the end."}

        # Serialize the result
        if isinstance(result, dict):
            return {k: _serialize(v) for k, v in result.items()}
        return {"value": _serialize(result)}

    except Exception as e:
        logger.error("Sandbox error: %s\nCode:\n%s", e, code[:500])
        return {"error": f"Sandbox error: {type(e).__name__}: {e}", "code_attempted": code[:300]}
