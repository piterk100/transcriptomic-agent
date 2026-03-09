# Transcriptomic Agent

Autonomous AI agent for discovering transcriptomic relationships across multiple datasets simultaneously.

## Getting Started

```bash
# Frontend
npm install
npm start
```

```bash
# Backend (separate terminal)
source .venv/bin/activate
pip install -r backend/requirements.txt
uvicorn backend.main:app --reload
```

Set `ANTHROPIC_API_KEY` in `backend/.env` (see `backend/.env.example`).

Open [http://localhost:3000](http://localhost:3000)

## Project Structure

```
src/
├── components/
│   ├── DatasetSlot.jsx  # file upload slot component
│   └── LogEntry.jsx     # agent log entry component
├── App.jsx              # main component, UI
└── index.js             # entry point

backend/
├── tools/
│   ├── single.py        # single-dataset tools
│   ├── cross.py         # cross-dataset tools (IRM, meta-DE, rewiring...)
│   ├── sandbox.py       # execute_code sandbox
│   └── registry.py      # tool registry + summarizeResult
├── agent/
│   ├── system_prompt.py # agent system prompt builder
│   ├── seeder.py        # pre-analysis seed hypothesis generator
│   └── runner.py        # agent loop (LLM → tool selection → result → loop)
└── main.py              # FastAPI app
```

## Adding a New Tool

1. Write a function in `backend/tools/single.py` or `backend/tools/cross.py`
   - Signature: `fn(datasets, **params) -> dict`
   - Return a plain dict with results and an `interpretation` string field

2. Register it in `backend/tools/registry.py`: add to `TOOLS` dict and add a summary lambda in `summarize_result()`

3. Document it in `backend/agent/system_prompt.py` — add the tool to the appropriate list

The agent will automatically see it and be able to call it.

## Pathway Enrichment

pathway_enrichment uses a local GMT file for fully reproducible,
publication-quality analysis.

Download GMT files from MSigDB (free, requires registration):
https://www.gsea-msigdb.org/gsea/msigdb/

Recommended files:
- h.all.v2023.2.Hs.symbols.gmt       — MSigDB Hallmarks (50 gene sets)
- c2.cp.kegg.v2023.2.Hs.symbols.gmt  — KEGG pathways (~186 gene sets)
- c2.cp.reactome.v2023.2.Hs.symbols.gmt — Reactome (~1600 gene sets)

To use multiple GMT files, concatenate them:
    cat h.all.v2023.2.Hs.symbols.gmt \
        c2.cp.kegg.v2023.2.Hs.symbols.gmt > combined.gmt

Set the path before starting the backend:
    export GMT_FILE=/path/to/combined.gmt
    uvicorn main:app --reload --port 8000

The GMT file version should be reported in your methods section, e.g.:
"Pathway enrichment was performed using a hypergeometric test with
Benjamini-Hochberg correction against MSigDB Hallmark gene sets
(v2023.2) and KEGG pathways (v2023.2)."

## Tools

### Single-dataset
| Tool | Description |
|---|---|
| `dataset_summary` | Overview of all loaded datasets |
| `top_variable_genes` | Top N genes by highest variance |
| `differential_expression` | DE between two groups (MWU + BH) |
| `gene_expression_by_group` | Mean expression of genes per group |
| `nonlinear_rule` | AND rule: GENE_A > threshold AND GENE_B <= threshold → group |
| `contextual_modules` | Context-dependent co-expression (mediated by a gene) |
| `pathway_enrichment` | Enrichment against gene sets from local GMT file (hypergeometric + BH) |
| `batch_detection` | Batch artifact detection for discovered axes |
| `subgroup_discovery` | Subgroups within a single group (PCA + KMeans) |
| `gene_network_hub` | Co-expression network hubs |

### Cross-dataset
| Tool | Description |
|---|---|
| `cross_dataset_de` | Genes DE in the same direction across ALL datasets (Fisher's method) |
| `cross_dataset_correlation` | Co-expression module replication across cohorts |
| `invariant_axis` | Biologically invariant axis (IRM-style, Cohen's d stability + bootstrap) |
| `cross_dataset_rewiring` | Correlation change between a gene pair across datasets |

### Dynamic
| Tool | Description |
|---|---|
| `execute_code` | Agent writes custom Python when no existing tool is sufficient |
