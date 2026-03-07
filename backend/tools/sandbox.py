import numpy as np
import pandas as pd
from scipy import stats


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


def execute_sandbox(code: str, datasets: list) -> dict:
    """
    Executes Python code written by the agent.
    Available variables: datasets, np, pd, stats
    Code must set the `result` variable before finishing.
    """
    context = {
        "datasets": datasets,
        "np": np,
        "pd": pd,
        "stats": stats,
    }
    local_vars = {}
    try:
        exec(compile(code, "<agent_code>", "exec"), context, local_vars)  # noqa: S102
        # Try result from local_vars first, then context (in case agent assigned to global scope)
        result = local_vars.get("result", context.get("result"))
        # Fallback: collect all non-private variables if result not set
        if result is None:
            candidates = {k: v for k, v in local_vars.items() if not k.startswith("_")}
            if candidates:
                result = {k: _serialize(v) for k, v in candidates.items() if _serialize(v) is not None}
        if result is None:
            return {"error": "Code did not set the `result` variable"}
        # Serialize the result
        if isinstance(result, dict):
            return {k: _serialize(v) for k, v in result.items()}
        return {"value": _serialize(result)}
    except Exception as e:
        return {"error": f"Sandbox error: {type(e).__name__}: {e}", "code_attempted": code[:300]}
