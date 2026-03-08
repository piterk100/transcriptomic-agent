import asyncio
import json
import logging
import re
from typing import AsyncGenerator

import anthropic

logger = logging.getLogger(__name__)


def _repair_json(s: str) -> str:
    """Replace literal control characters inside JSON strings with proper escape sequences."""
    result = []
    in_string = False
    escape = False
    for char in s:
        if escape:
            result.append(char)
            escape = False
        elif char == "\\" and in_string:
            result.append(char)
            escape = True
        elif char == '"':
            result.append(char)
            in_string = not in_string
        elif in_string and char == "\n":
            result.append("\\n")
        elif in_string and char == "\r":
            result.append("\\r")
        elif in_string and char == "\t":
            result.append("\\t")
        else:
            result.append(char)
    return "".join(result)


def _extract_first_json_object(s: str) -> str | None:
    """
    Extract the first syntactically balanced JSON object from s.
    Unlike a greedy regex, this correctly handles nested braces and strings,
    and stops at the first complete top-level closing brace.
    """
    start = s.find("{")
    if start == -1:
        return None
    depth = 0
    in_string = False
    escape = False
    for i, char in enumerate(s[start:], start):
        if escape:
            escape = False
            continue
        if char == "\\" and in_string:
            escape = True
            continue
        if char == '"':
            in_string = not in_string
            continue
        if not in_string:
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    return s[start : i + 1]
    return None

from ..agent.system_prompt import build_system_prompt
from ..agent.seeder import generate_seeds, extract_evidence_stats
from ..tools.registry import TOOLS, CROSS_TOOL_NAMES, summarize_result
from ..tools.sandbox import execute_sandbox


async def run_agent_loop(
    datasets: list,
    max_steps: int,
    api_key: str,
) -> AsyncGenerator[dict, None]:
    """
    Async generator — yields log event dicts.
    Consumed by FastAPI StreamingResponse.
    """
    client = anthropic.AsyncAnthropic(api_key=api_key)

    all_genes = [set(ds["expr"].index) for ds in datasets]
    common_genes = all_genes[0].intersection(*all_genes[1:]) if len(all_genes) > 1 else all_genes[0] if all_genes else set()

    # ── Pre-analysis seeding ──────────────────────────────────────────────────
    loop = asyncio.get_event_loop()
    seeds, seed_summary = await loop.run_in_executor(None, generate_seeds, datasets)
    system_prompt = build_system_prompt(datasets, len(common_genes), seed_summary=seed_summary)
    messages = []
    discoveries = []
    hypotheses: list[dict] = list(seeds)   # pre-populated with S1, S2, ...
    hypo_counter = 0                        # agent-proposed: H1, H2, ...

    yield {"type": "seed", "text": f"Pre-analysis: {len(seeds)} seed hypotheses generated", "summary": seed_summary}
    for s in seeds:
        yield {"type": "hypothesis_propose", "hypothesis": dict(s)}

    for i in range(max_steps + 1):
        step_num = i + 1

        discovery_summary = (
            "Discoveries:\n" + "\n".join(f"- [{d['action']}] {d['summary']}" for d in discoveries[-8:])
            if discoveries else "First step."
        )
        hypo_summary = (
            "\nHYPOTHESES:\n" + "\n".join(f"  [{h['id']}][{h['status'].upper()}] {h['text']}" for h in hypotheses)
            if hypotheses else ""
        )

        messages.append({
            "role": "user",
            "content": f"Step {step_num}/{max_steps}. {discovery_summary}{hypo_summary}\n\nWhat will you investigate?",
        })

        yield {"type": "thinking", "text": f"Agent thinking... ({step_num}/{max_steps})"}

        raw = ""
        try:
            async with client.messages.stream(
                model="claude-sonnet-4-20250514",
                max_tokens=16000,
                system=system_prompt,
                messages=messages,
            ) as stream:
                chunk_buf = ""
                last_flush = asyncio.get_event_loop().time()
                async for text in stream.text_stream:
                    raw += text
                    chunk_buf += text
                    now = asyncio.get_event_loop().time()
                    if len(chunk_buf) >= 40 or (now - last_flush) >= 0.12:
                        yield {"type": "thought_stream", "delta": chunk_buf}
                        chunk_buf = ""
                        last_flush = now
                if chunk_buf:
                    yield {"type": "thought_stream", "delta": chunk_buf}
        except Exception as e:
            logger.error("API error at step %d: %s", step_num, e, exc_info=True)
            yield {"type": "error", "text": f"API error: {e}"}
            return

        m = _extract_first_json_object(_repair_json(raw))
        try:
            dec = json.loads(m if m else raw)
        except json.JSONDecodeError:
            logger.error("JSON parse error at step %d: %s", step_num, raw[:200])
            yield {"type": "error", "text": f"JSON parse error: {raw[:200]}"}
            messages.append({"role": "assistant", "content": raw})
            continue

        thought = dec.get("thought", "")
        action = dec.get("action", "")
        params = dec.get("params", {})
        hypo_action = dec.get("hypothesis_action")

        # Process propose BEFORE tool call
        if isinstance(hypo_action, dict) and hypo_action.get("type") == "propose" and hypo_action.get("text"):
            hypo_counter += 1
            h = {
                "id": f"H{hypo_counter}",
                "text": hypo_action["text"],
                "genes": hypo_action.get("genes", []),   # optional list of genes the hypothesis is about
                "status": "pending",
                "evidence": [],
                "proposed_at": step_num,
                "seeded_by": "llm",
            }
            hypotheses.append(h)
            yield {"type": "hypothesis_propose", "hypothesis": dict(h)}

        if action == "DONE":
            yield {"type": "done", "text": thought}
            return

        # Guard: model sometimes puts "hypothesis_action" in the action field by mistake
        if action == "hypothesis_action":
            logger.warning("Agent used 'hypothesis_action' as action at step %d — skipping tool call", step_num)
            messages.append({"role": "assistant", "content": raw})
            messages.append({"role": "user", "content": "ERROR: 'hypothesis_action' is not a valid tool name. The 'action' field must be a tool name like differential_expression, execute_code, etc. The 'hypothesis_action' is a separate JSON field. Please retry with a valid tool."})
            continue

        yield {"type": "thought", "text": thought}
        await asyncio.sleep(0)  # flush thought before running tool

        loop = asyncio.get_event_loop()
        result = None
        if action == "execute_code":
            code = params.get("code", "").strip()
            if not code:
                logger.warning("execute_code at step %d: empty code parameter", step_num)
                messages.append({"role": "assistant", "content": raw})
                messages.append({"role": "user", "content": "ERROR: execute_code requires a non-empty 'code' parameter. Please provide the Python code to execute."})
                continue
            yield {"type": "code", "code": code}
            result = await loop.run_in_executor(None, execute_sandbox, code, datasets)
            summary = f"ERROR: {result['error']}" if isinstance(result, dict) and result.get("error") else f"Code executed: {str(result)[:80]}"
            discoveries.append({"action": action, "params": params, "summary": summary, "result": result})
            yield {"type": "result", "action": action, "params": params, "result": result, "summary": summary, "isCross": False, "isDynamic": True}
        elif action in TOOLS:
            try:
                tool_fn = TOOLS[action]
                result = await loop.run_in_executor(None, lambda: tool_fn(datasets, **params))
                summary = summarize_result(action, result)
                discoveries.append({"action": action, "params": params, "summary": summary, "result": result})
                yield {
                    "type": "result", "action": action, "params": params, "result": result,
                    "summary": summary, "isCross": action in CROSS_TOOL_NAMES, "isDynamic": False,
                }
            except Exception as e:
                logger.error("Tool error [%s] at step %d: %s", action, step_num, e, exc_info=True)
                result = {"error": str(e)}
                yield {"type": "error", "text": f"{action}: {e}"}
        else:
            logger.error("Unknown tool [%s] at step %d", action, step_num)
            result = {"error": f"Unknown tool: {action}"}
            yield {"type": "error", "text": f"Unknown tool: {action}"}

        # Process evaluate AFTER tool call
        if isinstance(hypo_action, dict) and hypo_action.get("type") == "evaluate" and hypo_action.get("hypothesis_id"):
            h = next((x for x in hypotheses if x["id"] == hypo_action["hypothesis_id"]), None)
            if h:
                verdict = hypo_action.get("verdict", "uncertain")
                if verdict not in ("confirmed", "rejected", "uncertain"):
                    verdict = "uncertain"
                h["status"] = verdict
                h["evidence"].append({
                    "step": step_num,
                    "action": action,
                    "reasoning": hypo_action.get("reasoning", ""),
                    "key_stats": extract_evidence_stats(action, result or {}, h.get("genes", [])),
                })
                yield {
                    "type": "hypothesis_eval",
                    "hypothesis": dict(h),
                    "reasoning": hypo_action.get("reasoning", ""),
                }

        messages.append({"role": "assistant", "content": raw})
        result_str = json.dumps(result, default=str)[:2500] if result is not None else "null"
        messages.append({"role": "user", "content": f"Wynik {action}:\n{result_str}"})
