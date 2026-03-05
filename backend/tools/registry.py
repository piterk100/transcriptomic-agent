from .single import (
    dataset_summary, top_variable_genes, differential_expression,
    gene_expression_by_group, nonlinear_rule, contextual_modules,
    pathway_enrichment, batch_detection, subgroup_discovery, gene_network_hub,
)
from .cross import (
    cross_dataset_de, cross_dataset_correlation, invariant_axis, cross_dataset_rewiring,
)

TOOLS = {
    "dataset_summary":           dataset_summary,
    "top_variable_genes":        top_variable_genes,
    "differential_expression":   differential_expression,
    "gene_expression_by_group":  gene_expression_by_group,
    "nonlinear_rule":            nonlinear_rule,
    "contextual_modules":        contextual_modules,
    "pathway_enrichment":        pathway_enrichment,
    "batch_detection":           batch_detection,
    "subgroup_discovery":        subgroup_discovery,
    "gene_network_hub":          gene_network_hub,
    "cross_dataset_de":          cross_dataset_de,
    "cross_dataset_correlation": cross_dataset_correlation,
    "invariant_axis":            invariant_axis,
    "cross_dataset_rewiring":    cross_dataset_rewiring,
}

CROSS_TOOL_NAMES = {
    "cross_dataset_de", "cross_dataset_correlation", "invariant_axis", "cross_dataset_rewiring"
}


def summarize_result(action: str, r: dict) -> str:
    if isinstance(r, dict) and r.get("error"):
        return f"BŁĄD: {r['error']}"
    try:
        m = {
            "dataset_summary":           lambda: f"{len(r)} datasetów",
            "top_variable_genes":        lambda: f"top {len(r)}: {', '.join(x['gene'] for x in r[:4])}",
            "differential_expression":   lambda: f"{r['n_significant']} DE | {r['comparison']}",
            "gene_expression_by_group":  lambda: f"{len(r.get('result', {}))} genów w {r.get('dataset', '?')}",
            "nonlinear_rule":            lambda: f"F1={r['f1']} | {r['rule']}",
            "contextual_modules":        lambda: f"spec={r['specificity']} | {r['interpretation']}",
            "pathway_enrichment":        lambda: r["interpretation"],
            "batch_detection":           lambda: r["interpretation"],
            "subgroup_discovery":        lambda: r["interpretation"],
            "gene_network_hub":          lambda: f"{r['interpretation']} | {r['n_edges']} krawędzi",
            "cross_dataset_de":          lambda: f"{r['n_consistent_genes']} spójnych genów w {len(r.get('datasets_used', []))} datasetach",
            "cross_dataset_correlation": lambda: f"replikacja min={r['min_replication']} | {r['interpretation']}",
            "invariant_axis":            lambda: f"inwariantna w {r['axis_consistent_in']} | {r['interpretation']}",
            "cross_dataset_rewiring":    lambda: f"{r['gene_pair']}: Δr={r['rewiring_magnitude']} | {r['interpretation']}",
        }
        return m[action]() if action in m else str(r)[:80]
    except Exception:
        return str(r)[:80]
