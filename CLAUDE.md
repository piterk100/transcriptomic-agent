# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Architecture

This is a full-stack application with two independent servers that must both run simultaneously:

- **Frontend**: React (Create React App) on port 3000, proxies `/api/*` to the backend
- **Backend**: FastAPI (Python) on port 8000, serves the AI agent and data processing

### Backend structure (`backend/`)

- `main.py` — FastAPI app with two endpoints: `POST /api/datasets` (file upload) and `POST /api/run` (streaming SSE agent loop). Datasets are stored in-memory (`_datasets` dict); they are lost on server restart.
- `agent/runner.py` — Core async generator `run_agent_loop()` that drives the agentic loop: calls Claude (`claude-sonnet-4-20250514`), parses JSON responses, dispatches tool calls, and yields SSE events.
- `agent/seeder.py` — Pre-analysis that runs before the loop: generates seed hypotheses (S1, S2, ...) by running quick cross-dataset DE on group pairs present in ≥2 datasets.
- `agent/system_prompt.py` — Builds the system prompt injected into every Claude call, including dataset descriptions, seed summary, and full tool documentation.
- `tools/single.py` — Single-dataset tools (t-test DE with BH correction, co-expression, PCA subgroups, pathway enrichment against hardcoded Hallmarks/KEGG gene sets, etc.).
- `tools/cross.py` — Cross-dataset tools (Fisher's method for meta-DE, IRM-style invariant axis via Cohen's d stability, correlation rewiring).
- `tools/registry.py` — `TOOLS` dict mapping tool names to functions; `summarize_result()` for log display.
- `tools/sandbox.py` — `execute_sandbox()` runs arbitrary Python code written by the agent using `exec()` with `datasets`, `np`, `pd`, `stats` in scope.

### Agent loop protocol

Each step the agent responds with strict JSON:
```json
{"thought": "...", "action": "tool_name", "params": {...}, "hypothesis_action": {...} | null}
```
Hypothesis IDs: `S1..Sn` = seeder-generated (pre-loop), `H1..Hn` = agent-proposed during loop.

The backend streams SSE events of type: `seed`, `thinking`, `thought`, `result`, `code`, `error`, `done`, `hypothesis_propose`, `hypothesis_eval`, `stream_end`.

### Frontend structure (`src/`)

- `App.jsx` — Single-component app. Manages dataset slots, loaded dataset metadata, group column selection, and the SSE stream reader. Three-panel layout: datasets panel (left), agent log (center), hypotheses panel (right).
- `components/DatasetSlot.jsx` — File upload UI for expression matrix + metadata CSV pair.
- `components/LogEntry.jsx` — Renders individual SSE events in the log.

### Data format

Datasets require two CSV files:
- **Expression matrix**: genes as rows, samples as columns, first column is gene names (index)
- **Metadata/phenotype**: samples as rows, columns are sample attributes (index = sample names); group columns are auto-detected as those with 2–10 unique values

## Running locally

**Backend** (requires Python 3.9+, `.venv` already exists):
```bash
source .venv/bin/activate
pip install -r backend/requirements.txt   # first time only
uvicorn backend.main:app --reload
```

**Frontend** (in a separate terminal):
```bash
npm install   # first time only
npm start
```

Set `ANTHROPIC_API_KEY` in `backend/.env` (copy from `backend/.env.example`).

## Adding a new tool

1. Implement the function in `backend/tools/single.py` or `backend/tools/cross.py`. Signature: `fn(datasets: list, **params) -> dict`. Return a dict with an `interpretation` string field.
2. Register it in `backend/tools/registry.py`: add to `TOOLS` dict and add a summary lambda in `summarize_result()`.
3. Document it in `backend/agent/system_prompt.py` under the appropriate section.
