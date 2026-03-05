import numpy as np
import pandas as pd
from scipy import stats


def execute_sandbox(code: str, datasets: list) -> dict:
    """
    Wykonuje kod Pythona napisany przez agenta.
    Dostępne zmienne: datasets, np, pd, stats
    Kod powinien ustawić zmienną `result`.
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
            return {"error": "Kod nie ustawił zmiennej `result`"}
        return result
    except Exception as e:
        return {"error": f"Sandbox error: {type(e).__name__}: {e}", "code_attempted": code[:300]}
