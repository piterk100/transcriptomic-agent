"""
Pre-analysis seeder: lightweight statistical analysis run before the agent loop to
generate data-driven seed hypotheses and provide structured evidence extraction.
IDs are S1, S2, ... (distinct from agent-proposed H1, H2, ...).
"""
from __future__ import annotations

import itertools

import numpy as np

from ..tools.single import top_variable_genes
from ..tools.cross import cross_dataset_de


def generate_seeds(datasets: list) -> tuple[list[dict], str, dict]:
    """
    Run quick pre-analysis and return (seed_hypotheses, seed_summary_text, seed_data).

    seed_hypotheses: list of hypothesis dicts with status="pending", seeded_by="auto".
    seed_summary_text: compact multi-line string injected into the system prompt.
    seed_data: full result data for the report (top variable genes + cross-dataset DE tables).
    """
    seeds: list[dict] = []
    summary_lines: list[str] = []
    seed_n = 0
    seed_data: dict = {"top_variable": [], "cross_de": []}

    try:
        # ── 1. Top variable genes per dataset ────────────────────────────────
        for ds in datasets:
            try:
                topv = top_variable_genes([ds], datasetName=ds["name"], n=10)
                gene_names = [g["gene"] for g in topv[:5]]
                if gene_names:
                    summary_lines.append(
                        f"  {ds['name']} — top variable: {', '.join(gene_names)}"
                    )
                seed_data["top_variable"].append({
                    "dataset": ds["name"],
                    "genes": topv[:10],
                })
            except Exception:
                pass

        # ── 2. Cross-dataset DE for group pairs present in >= 2 datasets ─────
        group_counts: dict[str, int] = {}
        for ds in datasets:
            for g in ds.get("groups", []):
                group_counts[g] = group_counts.get(g, 0) + 1

        multi_groups = [g for g, c in group_counts.items() if c >= 2]
        candidates: list[dict] = []

        for gA, gB in itertools.combinations(multi_groups, 2):
            try:
                result = cross_dataset_de(datasets, groupA=gA, groupB=gB, topN=5)
                if "error" in result:
                    continue
                seed_data["cross_de"].append({
                    "groupA": gA,
                    "groupB": gB,
                    "top_up": result.get("top_consistent_up", [])[:5],
                    "top_down": result.get("top_consistent_down", [])[:5],
                    "n_tested": result.get("n_genes_tested"),
                    "interpretation": result.get("interpretation", ""),
                })
                for entry in result.get("top_consistent_up", [])[:3]:
                    candidates.append({**entry, "_groupA": gA, "_groupB": gB, "_dir": "UP"})
                for entry in result.get("top_consistent_down", [])[:3]:
                    candidates.append({**entry, "_groupA": gA, "_groupB": gB, "_dir": "DOWN"})
            except Exception:
                pass

        # Sort by fisher_adj_p (ascending), take top 6
        candidates.sort(key=lambda x: x.get("fisher_adj_p", 1.0))
        for entry in candidates[:6]:
            seed_n += 1
            gene = entry["gene"]
            gA, gB = entry["_groupA"], entry["_groupB"]
            direction_str = "overexpressed" if entry["_dir"] == "UP" else "underexpressed"
            n_ds = entry.get("n_datasets", len(datasets))
            adj_p = entry.get("fisher_adj_p")
            lfc = entry.get("avg_abs_logFC")

            text = (
                f"{gene} is consistently {direction_str} in {gA} vs {gB} "
                f"across all {n_ds} datasets"
                + (f" (avg_|logFC|={lfc}, fisher_adj_p={adj_p})" if lfc and adj_p else "")
            )
            seeds.append({
                "id": f"S{seed_n}",
                "text": text,
                "status": "pending",
                "evidence": [],
                "proposed_at": 0,
                "seeded_by": "auto",
                "genes": [gene],
            })
            summary_lines.append(
                f"  S{seed_n}: {gene} {entry['_dir']} {gA}↑vs↓{gB} "
                f"| fisher_adj_p={adj_p} | n_datasets={n_ds}"
            )

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
        # Include if either gene in pair matches
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
