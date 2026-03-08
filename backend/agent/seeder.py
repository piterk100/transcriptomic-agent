"""
Pre-analysis seeder: proper statistical pre-analysis run before the agent loop to
generate data-driven seed hypotheses and provide structured evidence extraction.
IDs are S1, S2, ... (distinct from agent-proposed H1, H2, ...).
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu as mwu
from statsmodels.stats.multitest import multipletests

from ..tools.cross import cross_dataset_de


def generate_seeds(datasets: list) -> tuple[list[dict], str, dict]:
    """
    Run genome-wide MWU + BH pre-analysis per dataset/group-pair and return
    (seed_hypotheses, seed_summary_text, seed_data).

    seed_hypotheses: list of hypothesis dicts with status="pending", seeded_by="auto".
    seed_summary_text: compact multi-line string injected into the system prompt.
    seed_data: full result tables for the report.
    """
    seeds: list[dict] = []
    summary_lines: list[str] = []
    seed_id = 1
    seed_data: dict = {"per_dataset_de": [], "cross_de": []}

    # Convert list → dict keyed by dataset name for iteration
    ds_dict = {ds["name"]: ds for ds in datasets}

    try:
        # ── 1. Per-dataset differential expression (MWU + BH) ────────────────
        for ds_name, ds in ds_dict.items():
            expr = ds["expr"]
            meta = ds["meta"]
            gc = ds["group_col"]
            if not gc or gc not in meta.columns:
                continue
            groups = meta[gc].dropna().unique().tolist()

            for group_a, group_b in itertools.combinations(groups, 2):
                sA = meta[meta[gc] == group_a].index.intersection(expr.columns)
                sB = meta[meta[gc] == group_b].index.intersection(expr.columns)
                if len(sA) < 3 or len(sB) < 3:
                    summary_lines.append(
                        f"  {ds_name} — {group_a} vs {group_b}: skipped (n<3)"
                    )
                    continue

                # MWU per gene
                results = []
                for gene in expr.index:
                    a_vals = expr.loc[gene, sA].values.astype(float)
                    b_vals = expr.loc[gene, sB].values.astype(float)
                    if np.std(a_vals) == 0 and np.std(b_vals) == 0:
                        continue
                    try:
                        _, p = mwu(a_vals, b_vals)
                    except Exception:
                        continue
                    logfc = float(a_vals.mean() - b_vals.mean())
                    results.append({"gene": gene, "logFC": logfc, "p": p})

                if not results:
                    continue

                df = pd.DataFrame(results)
                _, adj_p, _, _ = multipletests(df["p"].values, method="fdr_bh")
                df["adj_p"] = adj_p
                sig = df[(df["adj_p"] < 0.05) & (df["logFC"].abs() > 0.5)].sort_values("adj_p")

                n_sig = len(sig)
                top_up   = sig[sig["logFC"] > 0].head(5)["gene"].tolist()
                top_down = sig[sig["logFC"] < 0].head(5)["gene"].tolist()
                top_genes = top_up + top_down

                summary_lines.append(
                    f"  {ds_name} — {group_a} vs {group_b}: "
                    f"{n_sig} DE genes (adj_p<0.05, |logFC|>0.5)"
                    + (f", top UP: {', '.join(top_up)}" if top_up else "")
                    + (f", top DOWN: {', '.join(top_down)}" if top_down else "")
                )

                # Store full table for report
                seed_data["per_dataset_de"].append({
                    "dataset": ds_name,
                    "groupA": group_a,
                    "groupB": group_b,
                    "n_sig": n_sig,
                    "top_up": sig[sig["logFC"] > 0].head(10)[["gene", "logFC", "adj_p"]].to_dict("records"),
                    "top_down": sig[sig["logFC"] < 0].head(10)[["gene", "logFC", "adj_p"]].to_dict("records"),
                })

                if top_genes:
                    description = (
                        f"{n_sig} genes DE between {group_a} and {group_b} "
                        f"in {ds_name} (MWU + BH, adj_p<0.05, |logFC|>0.5). "
                        f"Top UP: {', '.join(top_up) or 'none'}. "
                        f"Top DOWN: {', '.join(top_down) or 'none'}."
                    )
                    seeds.append({
                        "id": f"S{seed_id}",
                        "text": description,
                        "status": "pending",
                        "evidence": [],
                        "proposed_at": 0,
                        "seeded_by": "auto",
                        "genes": top_genes,
                    })
                    seed_id += 1
                else:
                    # No DE — seed a "no signal" hypothesis worth investigating
                    top_var = expr.var(axis=1).sort_values(ascending=False).head(5).index.tolist()
                    summary_lines.append(
                        f"  {ds_name} — {group_a} vs {group_b}: "
                        f"no DE found, top variable: {', '.join(top_var)}"
                    )
                    description = (
                        f"No significant DE between {group_a} and {group_b} in {ds_name} "
                        f"despite top variable genes: {', '.join(top_var)}. "
                        f"Investigate heterogeneity, subgroups, or effect sizes."
                    )
                    seeds.append({
                        "id": f"S{seed_id}",
                        "text": description,
                        "status": "pending",
                        "evidence": [],
                        "proposed_at": 0,
                        "seeded_by": "auto",
                        "genes": top_var,
                    })
                    seed_id += 1

        # ── 2. Cross-dataset DE if >= 2 datasets ─────────────────────────────
        if len(datasets) >= 2:
            all_groups: set[str] = set()
            for ds in datasets:
                all_groups.update(ds["meta"][ds["group_col"]].dropna().unique())

            for group_a, group_b in itertools.combinations(all_groups, 2):
                n_with_pair = sum(
                    1 for ds in datasets
                    if group_a in ds["meta"][ds["group_col"]].values
                    and group_b in ds["meta"][ds["group_col"]].values
                )
                if n_with_pair < 2:
                    continue
                try:
                    result = cross_dataset_de(datasets, groupA=group_a, groupB=group_b, topN=10)
                    if "error" in result:
                        continue
                    top_up_cross   = [g["gene"] for g in result.get("top_consistent_up", [])[:3]]
                    top_down_cross = [g["gene"] for g in result.get("top_consistent_down", [])[:3]]
                    top_consistent = top_up_cross + top_down_cross
                    n_consistent   = result.get("n_genes_tested", len(top_consistent))

                    seed_data["cross_de"].append({
                        "groupA": group_a,
                        "groupB": group_b,
                        "top_up": result.get("top_consistent_up", [])[:5],
                        "top_down": result.get("top_consistent_down", [])[:5],
                        "n_tested": result.get("n_genes_tested"),
                        "interpretation": result.get("interpretation", ""),
                    })

                    summary_lines.append(
                        f"  CROSS — {group_a} vs {group_b}: "
                        f"{n_with_pair} datasets, consistent genes: "
                        + (f"UP {', '.join(top_up_cross)}" if top_up_cross else "")
                        + (f" DOWN {', '.join(top_down_cross)}" if top_down_cross else "")
                    )
                    if top_consistent:
                        description = (
                            f"{len(top_consistent)} genes consistently DE across "
                            f"{n_with_pair} datasets ({group_a} vs {group_b}). "
                            f"Top: {', '.join(top_consistent)}. "
                            f"Investigate replicability and biological significance."
                        )
                        seeds.append({
                            "id": f"S{seed_id}",
                            "text": description,
                            "status": "pending",
                            "evidence": [],
                            "proposed_at": 0,
                            "seeded_by": "auto_cross",
                            "genes": top_consistent,
                        })
                        seed_id += 1
                except Exception:
                    pass

    except Exception:
        pass  # graceful degradation: agent starts without seeds

    seed_summary = (
        "Pre-analysis statistical results:\n" + "\n".join(summary_lines)
        if summary_lines else ""
    )
    return seeds, seed_summary, seed_data


def extract_evidence_stats(action: str, result: dict, genes: list[str]) -> dict:
    """
    Extract key statistics from a tool result that are relevant to the given genes.
    Returns {} if genes is empty or if result contains an error.
    Pure function — no API calls, fully deterministic.
    """
    if not genes or not isinstance(result, dict) or "error" in result:
        return {}

    gene_set = {g.upper() for g in genes}
    stats: dict[str, dict] = {}

    if action == "differential_expression":
        for dir_key, dir_label in [("top_upregulated", "UP"), ("top_downregulated", "DOWN")]:
            for entry in result.get(dir_key, []):
                if entry.get("gene", "").upper() in gene_set:
                    stats[entry["gene"]] = {
                        "logFC": entry.get("logFC"),
                        "rbc":   entry.get("rbc"),
                        "adj_p": entry.get("adj_p"),
                        "direction": dir_label,
                    }

    elif action == "cross_dataset_de":
        for dir_key, dir_label in [("top_consistent_up", "UP"), ("top_consistent_down", "DOWN")]:
            for entry in result.get(dir_key, []):
                if entry.get("gene", "").upper() in gene_set:
                    stats[entry["gene"]] = {
                        "avg_abs_logFC":  entry.get("avg_abs_logFC"),
                        "fisher_adj_p":   entry.get("fisher_adj_p"),
                        "n_sig_datasets": entry.get("n_sig_datasets"),
                        "direction":      dir_label,
                    }

    elif action == "invariant_axis":
        for entry in result.get("top_invariant_genes", []):
            if entry.get("gene", "").upper() in gene_set:
                stats[entry["gene"]] = {
                    "invariance_score":      entry.get("invariance_score"),
                    "mean_abs_effect":       entry.get("mean_abs_effect"),
                    "p_bootstrap":           entry.get("p_bootstrap"),
                    "bootstrap_significant": entry.get("bootstrap_significant"),
                }

    elif action == "gene_expression_by_group":
        for gene, grp_data in result.get("result", {}).items():
            if gene.upper() in gene_set and isinstance(grp_data, dict):
                stats[gene] = {
                    grp: {"mean": v.get("mean"), "std": v.get("std")}
                    for grp, v in grp_data.items()
                    if isinstance(v, dict)
                }

    elif action == "pathway_enrichment":
        for entry in result.get("top_enriched", []):
            overlap = [g for g in entry.get("overlap_genes", []) if g.upper() in gene_set]
            if overlap:
                stats[entry["pathway"]] = {
                    "enrichment_fold": entry.get("enrichment_fold"),
                    "adj_p":           entry.get("adj_p"),
                    "overlap":         overlap,
                }

    elif action == "gene_network_hub":
        for hub in result.get("top_hubs", []):
            if hub.get("gene", "").upper() in gene_set:
                stats[hub["gene"]] = {"hub_degree": hub.get("degree")}
        for edge in result.get("top_edges", []):
            g1, g2 = edge.get("g1", ""), edge.get("g2", "")
            if g1.upper() in gene_set or g2.upper() in gene_set:
                stats[f"{g1}–{g2}"] = {"spearman_r": edge.get("r")}

    elif action == "subgroup_discovery":
        for dir_key in ("top_markers_sub1", "top_markers_sub2"):
            for entry in result.get(dir_key, []):
                if entry.get("gene", "").upper() in gene_set:
                    stats[entry["gene"]] = {
                        "logFC": entry.get("logFC"),
                        "adj_p": entry.get("adj_p"),
                    }

    elif action == "cross_dataset_rewiring":
        pair = result.get("gene_pair", "")
        if any(g.upper() in gene_set for g in pair.replace("—", " ").replace("–", " ").split()):
            stats["rewiring"] = {
                "max_r":              result.get("max_r"),
                "min_r":              result.get("min_r"),
                "rewiring_magnitude": result.get("rewiring_magnitude"),
            }

    # Strip None values from leaf dicts for cleanliness
    return {
        gene: {k: v for k, v in s.items() if v is not None}
        for gene, s in stats.items()
        if isinstance(s, dict)
    }
