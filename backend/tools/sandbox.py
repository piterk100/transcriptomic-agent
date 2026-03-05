import numpy as np
import pandas as pd
from scipy import stats


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
        # Fallback: if no result, collect all non-private dicts/scalars from local_vars
        if result is None:
            candidates = {k: v for k, v in local_vars.items() if not k.startswith("_") and k != "result"}
            if candidates:
                result = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in candidates.items()
                          if isinstance(v, (int, float, str, list, dict, bool))}
        if not result:
            return {"error": "Code did not set the `result` variable"}
        return result
    except Exception as e:
        return {"error": f"Sandbox error: {type(e).__name__}: {e}", "code_attempted": code[:300]}
