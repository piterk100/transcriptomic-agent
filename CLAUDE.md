# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Architecture

This is a full-stack application with two independent servers that must both run simultaneously:

- **Frontend**: React (Create React App) on port 3000, proxies `/api/*` to the backend
- **Backend**: FastAPI (Python) on port 8000, serves the AI agent and data processing

### Backend structure (`backend/`)

- `main.py` — FastAPI app with two endpoints: `POST /api/datasets` (file upload) and `POST /api/run` (streaming SSE agent loop). Datasets are stored in-memory (`_datasets` dict); they are lost on server restart. `RunRequest` accepts `free_steps` (int) and `mode` (`"reproduce"` | `"explore"`); mode controls temperature (0.0 / 1.0).
- `agent/runner.py` — Core async generator `run_agent_loop()` that drives the agentic loop: calls Claude (`claude-sonnet-4-20250514`) at the configured temperature, parses JSON responses, dispatches tool calls, and yields SSE events. Writes a Markdown report to `reports/` on completion (or on step budget exhaustion).
- `agent/seeder.py` — Pre-analysis that runs before the loop: for each dataset and every pairwise group combination, runs genome-wide MWU + BH correction to find DE genes. Also runs cross-dataset DE for group pairs present in ≥2 datasets. Returns seed hypotheses (S1, S2, ...) and a `seed_data` dict used in the report.
- `agent/system_prompt.py` — Builds the system prompt injected into every Claude call, including dataset descriptions, seed summary, tool documentation, and statistical rules (no circular reasoning, no double-log, verdict criteria).
- `tools/single.py` — Single-dataset tools: `differential_expression` (MWU per gene + BH correction, consistent with seeder), `pathway_enrichment` (Hallmarks/KEGG hypergeometric, k≥3 required), `subgroup_discovery` (PCA + KMeans), `gene_network_hub` (permutation-based co-expression), and others.
- `tools/cross.py` — Cross-dataset tools: `cross_dataset_de` (Fisher's method for meta-DE), `invariant_axis` (Cohen's d stability), `cross_dataset_correlation`, `cross_dataset_rewiring`.
- `tools/registry.py` — `TOOLS` dict mapping tool names to functions; `summarize_result()` for log display.
- `tools/sandbox.py` — `execute_sandbox()` runs arbitrary Python code written by the agent using `exec()` with `datasets`, `np`, `pd`, `stats` in scope.

### Agent loop protocol

Each step the agent responds with strict JSON (action first, thought last):
```json
{"action": "tool_name", "params": {...}, "hypothesis_action": {...} | null, "thought": "..."}
```
Hypothesis IDs: `S1..Sn` = seeder-generated (pre-loop), `H1..Hn` = agent-proposed during loop.

The last step in the budget shows "FINAL STEP" instruction, prompting the agent to call DONE. If the agent exhausts the budget without DONE, a fallback report is written and a warning event is emitted.

The backend streams SSE events of type: `mode`, `seed`, `thinking`, `thought`, `thought_stream`, `result`, `code`, `error`, `done`, `hypothesis_propose`, `hypothesis_eval`, `report`, `stream_end`.

### Run modes

- **REPRODUCE** (`mode="reproduce"`, temperature=0): fully deterministic, same results on every run
- **EXPLORE** (`mode="explore"`, temperature=1): creative exploration, more varied hypotheses

### Frontend structure (`src/`)

- `App.jsx` — Single-component app. Manages dataset slots, loaded dataset metadata, group column selection, mode (REPRODUCE/EXPLORE), step count, and the SSE stream reader. Three-panel layout: datasets panel (left), agent log (center), hypotheses panel (right).
- `components/DatasetSlot.jsx` — File upload UI for expression matrix + metadata CSV pair.
- `components/LogEntry.jsx` — Renders individual SSE events in the log, including mode banner, seed pre-analysis block, hypothesis cards, and result expandable rows.

### Data format

Datasets require two CSV files:
- **Expression matrix**: genes as rows, samples as columns, first column is gene names (index). **Must be log-transformed** (log2 scale, typical range 3–14). The seeder asserts `max < 25` to catch raw counts.
- **Metadata/phenotype**: samples as rows, columns are sample attributes (index = sample names); group columns are auto-detected as those with 2–10 unique values

### Reports

Completed runs write a Markdown report to `reports/run_<timestamp>_<dataset>.md` containing:
- Pre-analysis tables (per-dataset MWU DE results + cross-dataset DE)
- All agent steps with thoughts, actions, results
- All hypotheses (S1..Sn seeds + H1..Hn agent-proposed) with evidence
- Final conclusion (agent's DONE summary or budget-exhausted warning)

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

Set `ANTHROPIC_API_KEY` in `backend/.env`.

## Language & Git conventions

- **All code, comments, strings, and documentation must be in English.** Never introduce Polish text anywhere in the codebase.
- **After every significant change, create a git commit and push to GitHub** (`git add`, `git commit`, `git push`). This allows the user to track history and revert. Do this automatically without being asked.
- GitHub repo: https://github.com/piterk100/transcriptomic-agent (branch: `main`)

## Adding a new tool

1. Implement the function in `backend/tools/single.py` or `backend/tools/cross.py`. Signature: `fn(datasets: list, **params) -> dict`. Return a dict with an `interpretation` string field.
2. Register it in `backend/tools/registry.py`: add to `TOOLS` dict and add a summary lambda in `summarize_result()`.
3. Document it in `backend/agent/system_prompt.py` under the appropriate section.
