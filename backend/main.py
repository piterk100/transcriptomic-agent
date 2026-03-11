from dotenv import load_dotenv
load_dotenv()

import io
import json
import logging
import os
import re
import uuid
from io import StringIO
from typing import Optional

import pandas as pd
from dotenv import load_dotenv
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

from .agent.runner import run_agent_loop

load_dotenv()

logging.basicConfig(
    filename="backend/app.log",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory dataset store: dataset_id → dataset dict
_datasets: dict[str, dict] = {}

# Pre-computed DEG table store: dataset_name → { name, type, comparisons[] }
deg_store: dict = {}

# Group name mappings: { "canonical_name": ["alias1", "alias2", ...] }
# Example: { "normal": ["Endometrium-Normal", "normal endometrium"] }
group_mappings: dict = {}

KEYWORDS_BIO = {"source", "condition", "characteristics", "tissue", "group", "disease"}


def _detect_group_cols(meta: pd.DataFrame) -> list[str]:
    cols = []
    for col in meta.columns:
        n_unique = meta[col].nunique()
        if 2 <= n_unique <= 10:
            cols.append(col)
    return cols


def _detect_best_group_col(group_cols: list[str]) -> str:
    for col in group_cols:
        if any(kw in col.lower() for kw in KEYWORDS_BIO):
            return col
    return group_cols[0] if group_cols else ""


@app.post("/api/datasets")
async def upload_dataset(
    expr_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
    name: str = Form(...),
    dataset_id: Optional[str] = Form(None),
):
    logger.info("Upload dataset: name=%s", name)
    expr_text = (await expr_file.read()).decode("utf-8")
    meta_text = (await meta_file.read()).decode("utf-8")

    expr = pd.read_csv(StringIO(expr_text), index_col=0)
    meta = pd.read_csv(StringIO(meta_text), index_col=0)

    # Align samples
    common_samples = [s for s in expr.columns if s in meta.index]
    expr = expr[common_samples]
    meta = meta.loc[common_samples]

    # Score each candidate column by number of unique values (2-10)
    # and presence of biological keywords in column name OR values
    kw_name = ["source", "condition", "characteristics", "tissue",
               "group", "disease", "type", "status", "treatment", "diagnosis"]
    kw_values = ["endometri", "adenomyo", "control", "normal", "disease",
                 "tumor", "cancer", "healthy", "patient", "case", "eutopic",
                 "ectopic", "peritoneal", "lesion", "ovary", "uterus"]

    candidates = [c for c in meta.columns
                  if meta[c].nunique() >= 2 and meta[c].nunique() <= 15]

    if not candidates:
        gc = meta.columns[0]
    else:
        def score_col(c):
            score = 0
            # Keyword in column name
            if any(k in c.lower() for k in kw_name):
                score += 2
            # Keyword in values
            sample_vals = " ".join(meta[c].dropna().astype(str).str.lower().tolist())
            if any(k in sample_vals for k in kw_values):
                score += 3
            # Prefer fewer unique values (more group-like)
            score += max(0, 5 - meta[c].nunique())
            return score

        gc = max(candidates, key=score_col)

    group_cols = candidates
    best_col = gc

    # Clean group values — strip common GEO prefixes
    def clean_group_values(series):
        s = series.copy().astype(str).str.strip()
        # Strip "key: value" prefixes like "disease state: ", "tissue: " etc.
        s = s.str.replace(r'^[\w\s]+:\s*', '', regex=True).str.strip()
        return s

    # If after cleaning we get fewer unique values, use cleaned version
    cleaned = clean_group_values(meta[gc])
    if cleaned.nunique() < meta[gc].nunique():
        meta[gc] = cleaned

    # If every sample has a unique value, try to extract Disease/Normal/Control
    if meta[gc].nunique() == len(meta):
        extracted = meta[gc].str.extract(
            r'(disease|normal|control|endometri|adenomyo|cancer|tumor|healthy)',
            flags=re.IGNORECASE, expand=False
        ).str.lower().str.strip()
        if extracted.notna().all() and extracted.nunique() >= 2:
            meta[gc] = extracted

    groups = meta[best_col].dropna().unique().tolist() if best_col else []

    ds_id = dataset_id or str(uuid.uuid4())
    _datasets[ds_id] = {
        "id": ds_id,
        "name": name,
        "expr": expr,
        "meta": meta,
        "group_col": best_col,
        "group_cols": group_cols,
        "groups": groups,
    }

    return {
        "id": ds_id,
        "name": name,
        "gene_count": len(expr.index),
        "sample_count": len(expr.columns),
        "group_cols": group_cols,
        "group_col": gc,
        "group_col_candidates": [
            {
                "col": c,
                "unique_values": meta[c].dropna().unique().tolist()[:10]
            }
            for c in meta.columns
            if 2 <= meta[c].nunique() <= 20
            and meta[c].dtype == object
        ],
        "groups": meta[gc].unique().tolist(),
    }



@app.post("/api/datasets/upload_deg")
async def upload_deg(
    file: UploadFile = File(...),
    name: str = Form(...),
    groupA: str = Form(...),
    groupB: str = Form(...),
):
    """
    Upload a pre-computed DEG table CSV.
    Required columns (auto-detected, case-insensitive):
      - gene / gene_symbol / gene_name / id
      - logFC / log2FC / log2FoldChange / logfoldchange
      - p / pvalue / p_value / pval
      - adj_p / padj / adj.p / p_adj / fdr / q_value
    groupA and groupB are the canonical group names (matching group_mappings).
    """
    content = await file.read()

    # Read CSV — handle limma format where genes are row index
    raw = pd.read_csv(io.BytesIO(content))

    # If first column is unnamed (limma rownames exported as index), use it as gene column
    if raw.columns[0].startswith("Unnamed:"):
        raw = raw.rename(columns={raw.columns[0]: "gene"})

    # Auto-detect columns (case-insensitive, handle dots and spaces)
    def norm(s):
        return s.lower().replace(".", "_").replace(" ", "_")

    col_map: dict[str, str] = {}
    normed = {norm(c): c for c in raw.columns}

    # Gene column — check in priority order
    for candidate in ["gene", "symbol", "gene_symbol", "gene_name", "id", "geneid"]:
        if candidate in normed and "gene" not in col_map:
            # Skip "gene_name" if "symbol" is also present (gene_name has full names)
            if candidate == "gene_name" and "symbol" in normed:
                continue
            col_map["gene"] = normed[candidate]

    # logFC column
    for candidate in ["logfc", "log2fc", "log2foldchange", "logfoldchange", "lfc", "log2_fold_change"]:
        if candidate in normed and "logFC" not in col_map:
            col_map["logFC"] = normed[candidate]

    # p-value column — detect adj_p first to avoid "p" matching "padj"
    for candidate in ["adj_p_val", "adj_p", "padj", "adj_pval", "fdr", "q_value", "qvalue", "p_adj"]:
        if candidate in normed and "adj_p" not in col_map:
            col_map["adj_p"] = normed[candidate]

    for candidate in ["p_value", "pvalue", "pval", "p"]:
        if candidate in normed and "p" not in col_map:
            col_map["p"] = normed[candidate]

    missing = [k for k in ["gene", "logFC", "p", "adj_p"] if k not in col_map]
    if missing:
        raise HTTPException(
            400,
            f"Could not detect columns: {missing}. "
            f"Found columns: {raw.columns.tolist()}. "
            f"Expected: gene/symbol, logFC/log2FoldChange, p/pvalue/P.Value, "
            f"adj_p/padj/adj.P.Val/FDR",
        )

    # Standardize
    df = raw.rename(columns={v: k for k, v in col_map.items()})
    df = df[["gene", "logFC", "p", "adj_p"]].dropna()
    df["gene"] = df["gene"].astype(str).str.strip().str.upper()
    df = df.set_index("gene")

    # Remove duplicates — keep row with lowest adj_p per gene
    df = df[~df.index.duplicated(keep="first")]

    if name not in deg_store:
        deg_store[name] = {"name": name, "type": "deg", "comparisons": []}

    # Remove existing entry for same groupA/groupB if present
    deg_store[name]["comparisons"] = [
        c for c in deg_store[name]["comparisons"]
        if not (c["groupA"] == groupA and c["groupB"] == groupB)
    ]
    deg_store[name]["comparisons"].append({"groupA": groupA, "groupB": groupB, "df": df})

    return {
        "name": name,
        "type": "deg",
        "n_genes": len(df),
        "groupA": groupA,
        "groupB": groupB,
        "detected_columns": col_map,
        "comparisons": [
            {"groupA": c["groupA"], "groupB": c["groupB"], "n_genes": len(c["df"])}
            for c in deg_store[name]["comparisons"]
        ],
    }


@app.get("/api/datasets/deg")
async def list_deg_datasets():
    return {
        name: {
            "name": name,
            "type": "deg",
            "comparisons": [
                {"groupA": c["groupA"], "groupB": c["groupB"], "n_genes": len(c["df"])}
                for c in ds["comparisons"]
            ],
        }
        for name, ds in deg_store.items()
    }


class RunRequest(BaseModel):
    dataset_ids: list[str]
    group_cols: dict[str, str]  # dataset_id → chosen group_col
    free_steps: int = 6
    mode: str = "reproduce"  # "reproduce" | "explore"


@app.post("/api/run")
async def run_agent(req: RunRequest):
    api_key = os.getenv("ANTHROPIC_API_KEY", "")
    if not api_key:
        raise HTTPException(status_code=500, detail="ANTHROPIC_API_KEY not set in .env")

    datasets = []
    for ds_id in req.dataset_ids:
        if ds_id not in _datasets:
            raise HTTPException(status_code=404, detail=f"Dataset {ds_id} nie istnieje")
        ds = dict(_datasets[ds_id])
        chosen_col = req.group_cols.get(ds_id, ds["group_col"])
        if chosen_col and chosen_col in ds["meta"].columns:
            ds["group_col"] = chosen_col
            ds["groups"] = ds["meta"][chosen_col].dropna().unique().tolist()
        datasets.append(ds)

    temperature = 0.0 if req.mode == "reproduce" else 1.0
    logger.info("Run started: datasets=%s max_steps=%d mode=%s", req.dataset_ids, req.free_steps, req.mode)

    async def generate():
        async for event in run_agent_loop(datasets, req.free_steps, api_key, temperature=temperature, mappings=dict(group_mappings), deg_datasets=dict(deg_store)):
            if event.get("type") == "error":
                logger.error("Agent error: %s", event.get("text", ""))
            yield f"data: {json.dumps(event, default=str)}\n\n"
        logger.info("Run finished: datasets=%s", req.dataset_ids)
        yield "data: {\"type\":\"stream_end\"}\n\n"

    return StreamingResponse(
        generate(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
        },
    )


@app.patch("/api/datasets/{dataset_id}/group_col")
async def update_group_col(dataset_id: str, body: dict):
    if dataset_id not in _datasets:
        raise HTTPException(status_code=404, detail="Dataset not found")
    gc = body.get("group_col")
    meta = _datasets[dataset_id]["meta"]
    if gc not in meta.columns:
        raise HTTPException(status_code=400, detail=f"Column {gc} not found")
    _datasets[dataset_id]["group_col"] = gc
    _datasets[dataset_id]["groups"] = meta[gc].unique().tolist()
    return {"group_col": gc, "groups": meta[gc].unique().tolist()}


@app.post("/api/group_mappings")
async def set_group_mappings(body: dict):
    """
    Body: { "mappings": { "canonical": ["alias1", "alias2"] } }
    Example: {
        "mappings": {
            "normal": ["Endometrium-Normal", "normal endometrium"],
            "disease": ["Endometrium/Ovary-Disease", "endometrium with endometriosis"]
        }
    }
    """
    group_mappings.clear()
    group_mappings.update(body.get("mappings", {}))
    return {"status": "ok", "mappings": group_mappings}


@app.get("/api/group_mappings")
async def get_group_mappings():
    return {"mappings": group_mappings}


@app.delete("/api/datasets/{dataset_id}")
async def delete_dataset(dataset_id: str):
    _datasets.pop(dataset_id, None)
    return {"ok": True}
