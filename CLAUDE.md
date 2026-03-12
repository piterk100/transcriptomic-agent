# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Architecture

This is a full-stack application with two independent servers that must both run simultaneously:

- **Frontend**: React (Create React App) on port 3000, proxies `/api/*` to the backend
- **Backend**: FastAPI (Python) on port 8000, serves the AI agent and data processing

### Backend structure (`backend/`)

- `main.py` — FastAPI app. Endpoints: `POST /api/datasets` (file upload), `PATCH /api/datasets/{id}/group_col`, `DELETE /api/datasets/{id}`, `POST /api/run` (streaming SSE agent loop), `POST /api/group_mappings`, `GET /api/group_mappings`, `POST /api/datasets/upload_deg`, `GET /api/datasets/deg`. Datasets stored in-memory (`_datasets` dict). Pre-computed DEG tables stored in `deg_store` dict. Group name mappings stored in `group_mappings` dict. `RunRequest` accepts `free_steps` (int) and `mode` (`"reproduce"` | `"explore"`); mode controls temperature (0.0 / 1.0). Calls `load_dotenv()` at startup. Group column autodetection uses keyword scoring (name keywords + value keywords + fewer-unique-values preference). GEO prefix cleaning strips `"key: value"` prefixes from group values.
- `agent/runner.py` — Core async generator `run_agent_loop(datasets, max_steps, api_key, temperature, mappings, deg_datasets)`. Guards: (1) blocks DONE when any hypothesis is still PENDING (except on the last step), (2) blocks duplicate tool calls with identical params in consecutive steps, (3) blocks tools in `_ONCE_ONLY` set (`cross_dataset_de`) after first call, (4) in DEG-only mode blocks all tools not in `_DEG_ONLY_ALLOWED`. Appends GROUP MAPPINGS block to system prompt when mappings are defined. Writes a Markdown report to `reports/` on completion or budget exhaustion.
- `agent/seeder.py` — Pre-analysis that runs before the loop: for each dataset and every pairwise group combination, runs genome-wide MWU + BH correction to find DE genes. Also runs cross-dataset DE for group pairs present in ≥2 datasets. Returns seed hypotheses (S1, S2, ...) and a `seed_data` dict used in the report.
- `agent/system_prompt.py` — Builds the system prompt injected into every Claude call, including dataset descriptions, DEG dataset descriptions, seed summary, tool documentation, and rules: no circular reasoning, no double-log, verdict criteria, no duplicate tool calls, cross_dataset_de once-only, subgroup homogeneity interpretation, EXECUTE_CODE RULES (forbidden uses + self-check question), immediate hypothesis evaluation after pathway_enrichment/execute_code.
- `tools/single.py` — Single-dataset tools: `differential_expression` (MWU per gene + BH correction, consistent with seeder), `pathway_enrichment` (hypergeometric + BH against local GMT file; accepts `deg_dataset_name` to auto-extract significant genes from a DEG dataset; requires `GMT_FILE` env var), `subgroup_discovery` (PCA + KMeans; min subgroup size = max(3, 15% of samples); includes multi-seed and subsample silhouette stability), `gene_network_hub` (permutation-based co-expression), and others.
- `tools/cross.py` — Cross-dataset tools: `cross_dataset_de` (Fisher's method for meta-DE; integrates uploaded DEG tables automatically), `invariant_axis` (Cohen's d stability; warns when n < 10), `cross_dataset_correlation`, `cross_dataset_rewiring`. All group-based tools accept `mappings` param and use `resolve_group()` to map aliases to canonical names.
- `tools/deg.py` — DEG table analysis tools: `deg_voting` (per-gene vote count + direction consistency across DEG tables), `deg_biomarker_ranking` (composite score = freq × consistency × mean_abs_logFC × mean(−log10 adj_p)), `deg_cooccurrence_network` (weighted gene co-occurrence graph; warns if < 3 comparisons), `deg_direction_comparison` (concordant/discordant/specific genes between two named comparisons).
- `tools/registry.py` — `TOOLS` dict mapping tool names to functions; `DEG_TOOL_NAMES` set; `summarize_result()` for log display.
- `tools/sandbox.py` — `execute_sandbox(code, datasets, deg_datasets)` runs arbitrary Python code written by the agent using `exec()` with `datasets`, `deg_datasets`, `np`, `pd`, `stats` in scope.

### DEG-only mode

When only pre-computed DEG tables are loaded (no raw expression datasets), the agent operates in DEG-only mode:
- System prompt shows a warning banner and lists only available tools
- Runner hard-blocks any tool not in `_DEG_ONLY_ALLOWED`: `{cross_dataset_de, pathway_enrichment, execute_code, DONE} | DEG_TOOL_NAMES`
- DEG datasets are named automatically: `DEG 1`, `DEG 2`, ... (frontend assigns names sequentially)

### DEG table upload

`POST /api/datasets/upload_deg` accepts a CSV with auto-detected columns (case-insensitive, handles `.` and space separators):
- Gene: `gene`, `symbol`, `gene_symbol`, `gene_name` (skipped if `symbol` present), `id`, `geneid`; unnamed first column treated as gene index (limma format)
- logFC: `logfc`, `log2fc`, `log2foldchange`, `logfoldchange`, `lfc`, `log2_fold_change`
- adj_p (detected before p to avoid `padj` matching p bucket): `adj_p_val`, `adj_p`, `padj`, `adj_pval`, `fdr`, `q_value`, `qvalue`, `p_adj`
- p: `p_value`, `pvalue`, `pval`, `p`

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

- `App.jsx` — Single-component app. Manages dataset slots, loaded dataset metadata, DEG dataset list, group column selection, group mappings, mode (REPRODUCE/EXPLORE), step count, and the SSE stream reader. Three-panel layout: datasets panel (left), agent log (center), hypotheses panel (right). Left panel order: DATASETS → DEG TABLES → GROUP COLUMNS → GROUP MAPPINGS → MODE → STEPS → START AGENT. Agent controls shown when at least one raw dataset or DEG dataset is loaded.
- `api.js` — `setGroupMappings(mappings)`, `uploadDegDataset(file, name, groupA, groupB)`.
- `components/DatasetSlot.jsx` — File upload UI for expression matrix + metadata CSV pair.
- `components/LogEntry.jsx` — Renders individual SSE events in the log, including mode banner, seed pre-analysis block, hypothesis cards, and result expandable rows.

### Data format

Datasets require two CSV files:
- **Expression matrix**: genes as rows, samples as columns, first column is gene names (index). **Must be log-transformed** (log2 scale, typical range 3–14). The seeder asserts `max < 25` to catch raw counts.
- **Metadata/phenotype**: samples as rows, columns are sample attributes (index = sample names); group columns are auto-detected from candidates with 2–15 unique values using keyword scoring; GEO-style `"key: value"` prefixes are automatically stripped from group values

### Pathway enrichment

`pathway_enrichment` uses a local GMT file (MSigDB/KEGG/Reactome). Set `GMT_FILE=/path/to/file.gmt` in `backend/.env` before starting the backend. If the path is relative, it is resolved relative to `backend/`. See README.md for download instructions.

### Group mappings

Group name equivalences across datasets (e.g. `"normal"` = `["Endometrium-Normal", "normal endometrium"]`) can be defined via `POST /api/group_mappings`. The frontend exposes a Group Mappings section in the left panel. Mappings are passed to all cross-dataset tool calls and appended to the agent's system prompt.

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

Set `ANTHROPIC_API_KEY` and optionally `GMT_FILE` in `backend/.env`.

## Language & Git conventions

- **All code, comments, strings, and documentation must be in English.** Never introduce Polish text anywhere in the codebase.
- **After every significant change, create a git commit and push to GitHub** (`git add`, `git commit`, `git push`). This allows the user to track history and revert. Do this automatically without being asked.
- GitHub repo: https://github.com/piterk100/transcriptomic-agent (branch: `main`)

## Adding a new tool

1. Implement the function in `backend/tools/single.py`, `backend/tools/cross.py`, or `backend/tools/deg.py`. Signature: `fn(datasets: list, **params) -> dict`. Return a dict with an `interpretation` string field.
2. Register it in `backend/tools/registry.py`: add to `TOOLS` dict, add to `DEG_TOOL_NAMES` if it operates on DEG tables only, and add a summary lambda in `summarize_result()`.
3. Document it in `backend/agent/system_prompt.py` under the appropriate section.
