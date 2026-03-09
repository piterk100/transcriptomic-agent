from dotenv import load_dotenv
load_dotenv()

import json
import logging
import os
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
        "group_col": best_col,
        "groups": groups,
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
        async for event in run_agent_loop(datasets, req.free_steps, api_key, temperature=temperature):
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


@app.delete("/api/datasets/{dataset_id}")
async def delete_dataset(dataset_id: str):
    _datasets.pop(dataset_id, None)
    return {"ok": True}
