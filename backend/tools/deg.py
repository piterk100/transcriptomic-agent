from collections import Counter, defaultdict

import numpy as np

from .cross import resolve_group


def deg_voting(datasets: list, deg_datasets: dict = None,
               groupA: str = None, groupB: str = None,
               mappings: dict = None,
               adj_p_threshold: float = 0.05,
               logfc_threshold: float = 0.5,
               topN: int = 30, **_) -> dict:
    """
    For each gene, count how many DEG tables show it as significant
    (adj_p < threshold, |logFC| > threshold) with consistent direction.

    If groupA/groupB provided, only include comparisons matching those groups.
    Otherwise use all available comparisons.
    """
    mappings = mappings or {}
    deg_datasets = deg_datasets or {}

    comparisons = []
    for ds_name, ds in deg_datasets.items():
        for comp in ds["comparisons"]:
            a = resolve_group(comp["groupA"], mappings)
            b = resolve_group(comp["groupB"], mappings)
            if groupA and groupB:
                if a == groupA and b == groupB:
                    comparisons.append((ds_name, comp, 1))
                elif a == groupB and b == groupA:
                    comparisons.append((ds_name, comp, -1))
            else:
                comparisons.append((ds_name, comp, 1))

    if not comparisons:
        return {"error": "No matching DEG comparisons found"}

    n_datasets = len(comparisons)
    gene_votes: dict = {}

    for ds_name, comp, direction in comparisons:
        df = comp["df"]
        sig = df[(df["adj_p"] < adj_p_threshold) & (df["logFC"].abs() > logfc_threshold)]
        for gene, row in sig.iterrows():
            if gene not in gene_votes:
                gene_votes[gene] = {"up": 0, "down": 0, "datasets": [], "logfcs": []}
            effective_logfc = float(row["logFC"]) * direction
            if effective_logfc > 0:
                gene_votes[gene]["up"] += 1
            else:
                gene_votes[gene]["down"] += 1
            gene_votes[gene]["datasets"].append(ds_name)
            gene_votes[gene]["logfcs"].append(round(effective_logfc, 3))

    if not gene_votes:
        return {"error": "No significant genes found across DEG tables"}

    results = []
    for gene, v in gene_votes.items():
        total = v["up"] + v["down"]
        consistent = v["up"] == total or v["down"] == total
        direction_label = "UP" if v["up"] >= v["down"] else "DOWN"
        consistency_score = max(v["up"], v["down"]) / total
        results.append({
            "gene": gene,
            "n_datasets": total,
            "freq": round(total / n_datasets, 3),
            "direction": direction_label,
            "consistent": consistent,
            "consistency_score": round(consistency_score, 3),
            "mean_logFC": round(sum(v["logfcs"]) / len(v["logfcs"]), 3),
            "datasets": v["datasets"],
        })

    results.sort(key=lambda x: (-x["n_datasets"], -x["consistency_score"], -abs(x["mean_logFC"])))
    top = results[:topN]
    fully_consistent = [r for r in results if r["consistent"]]

    return {
        "n_comparisons": n_datasets,
        "n_genes_any": len(results),
        "n_genes_fully_consistent": len(fully_consistent),
        "top_genes": top,
        "interpretation": (
            f"Top: {top[0]['gene']} ({top[0]['n_datasets']}/{n_datasets} datasets, "
            f"{top[0]['direction']}, mean logFC={top[0]['mean_logFC']})"
            if top else "No consistent genes found"
        ),
    }


def deg_cooccurrence_network(datasets: list, deg_datasets: dict = None,
                              groupA: str = None, groupB: str = None,
                              mappings: dict = None,
                              adj_p_threshold: float = 0.05,
                              logfc_threshold: float = 0.5,
                              min_cooccurrence: int = 2,
                              topN_genes: int = 50, **_) -> dict:
    """
    Build co-occurrence network from DEG tables.
    Edge between gene A and gene B if both are DE in the same comparison,
    with weight = number of comparisons where they co-occur.
    Only edges with weight >= min_cooccurrence are included.
    """
    mappings = mappings or {}
    deg_datasets = deg_datasets or {}

    comparisons_genes = []
    for ds_name, ds in deg_datasets.items():
        for comp in ds["comparisons"]:
            a = resolve_group(comp["groupA"], mappings)
            b = resolve_group(comp["groupB"], mappings)
            if groupA and groupB:
                if not ((a == groupA and b == groupB) or (a == groupB and b == groupA)):
                    continue
            df = comp["df"]
            sig_genes = set(
                df[(df["adj_p"] < adj_p_threshold) & (df["logFC"].abs() > logfc_threshold)].index
            )
            if sig_genes:
                comparisons_genes.append((ds_name, sig_genes))

    if not comparisons_genes:
        return {"error": "No significant genes found in any comparison"}

    low_comparisons_warning = None
    if len(comparisons_genes) < 3:
        low_comparisons_warning = (
            f"Co-occurrence network requires >= 3 comparisons for meaningful results. "
            f"Only {len(comparisons_genes)} found. Results will show trivially connected graph."
        )

    edge_weights: dict = defaultdict(int)
    gene_freq: dict = defaultdict(int)

    for ds_name, genes in comparisons_genes:
        gene_list = sorted(genes)
        for g in gene_list:
            gene_freq[g] += 1
        for i in range(len(gene_list)):
            for j in range(i + 1, len(gene_list)):
                edge_weights[(gene_list[i], gene_list[j])] += 1

    filtered_edges = {pair: w for pair, w in edge_weights.items() if w >= min_cooccurrence}

    degree: Counter = Counter()
    for (g1, g2), w in filtered_edges.items():
        degree[g1] += w
        degree[g2] += w

    top_hubs = degree.most_common(topN_genes)
    top_edges = sorted(filtered_edges.items(), key=lambda x: -x[1])[: topN_genes * 2]

    result = {
        "n_comparisons": len(comparisons_genes),
        "n_edges_total": len(filtered_edges),
        "n_hub_genes": len(degree),
        "top_hubs": [{"gene": g, "weighted_degree": d} for g, d in top_hubs[:20]],
        "top_edges": [{"gene1": e[0][0], "gene2": e[0][1], "weight": e[1]} for e in top_edges[:30]],
        "interpretation": (
            f"Top hub: {top_hubs[0][0]} (weighted degree={top_hubs[0][1]}) "
            f"| {len(filtered_edges)} edges with co-occurrence >= {min_cooccurrence}"
            if top_hubs else "No co-occurrence network found"
        ),
    }
    if low_comparisons_warning:
        result["warning"] = low_comparisons_warning
    return result


def deg_biomarker_ranking(datasets: list, deg_datasets: dict = None,
                           groupA: str = None, groupB: str = None,
                           mappings: dict = None,
                           adj_p_threshold: float = 0.05,
                           logfc_threshold: float = 0.5,
                           topN: int = 30, **_) -> dict:
    """
    Composite biomarker score per gene:
      score = freq * consistency * mean_abs_logFC * mean_(-log10_adj_p)

    Higher score = more consistently DE with larger effect across datasets.
    Only genes significant in >= 2 datasets are ranked.
    """
    mappings = mappings or {}
    deg_datasets = deg_datasets or {}

    gene_stats: dict = {}
    n_comparisons = 0

    for ds_name, ds in deg_datasets.items():
        for comp in ds["comparisons"]:
            a = resolve_group(comp["groupA"], mappings)
            b = resolve_group(comp["groupB"], mappings)
            if groupA and groupB:
                if a == groupA and b == groupB:
                    direction = 1
                elif a == groupB and b == groupA:
                    direction = -1
                else:
                    continue
            else:
                direction = 1
            n_comparisons += 1
            df = comp["df"]
            sig = df[(df["adj_p"] < adj_p_threshold) & (df["logFC"].abs() > logfc_threshold)]
            for gene, row in sig.iterrows():
                if gene not in gene_stats:
                    gene_stats[gene] = []
                gene_stats[gene].append({
                    "logFC": float(row["logFC"]) * direction,
                    "adj_p": max(float(row["adj_p"]), 1e-300),
                })

    if not gene_stats:
        return {"error": "No significant genes found"}

    results = []
    for gene, stats in gene_stats.items():
        if len(stats) < 1:
            continue
        logfcs = [s["logFC"] for s in stats]
        adj_ps = [s["adj_p"] for s in stats]
        n = len(stats)
        freq = n / n_comparisons
        mean_abs_logfc = float(np.mean([abs(l) for l in logfcs]))
        mean_neg_log_p = float(np.mean([-np.log10(p) for p in adj_ps]))
        up = sum(1 for l in logfcs if l > 0)
        consistency = max(up, n - up) / n
        direction = "UP" if up >= n - up else "DOWN"
        score = freq * consistency * mean_abs_logfc * mean_neg_log_p
        results.append({
            "gene": gene,
            "score": round(score, 4),
            "n_datasets": n,
            "freq": round(freq, 3),
            "direction": direction,
            "consistency": round(consistency, 3),
            "mean_abs_logFC": round(mean_abs_logfc, 3),
            "mean_neg_log10_adj_p": round(mean_neg_log_p, 3),
        })

    results.sort(key=lambda x: -x["score"])
    top = results[:topN]

    return {
        "n_comparisons": n_comparisons,
        "n_genes_ranked": len(results),
        "top_biomarkers": top,
        "interpretation": (
            f"Top biomarker: {top[0]['gene']} "
            f"(score={top[0]['score']}, {top[0]['n_datasets']} datasets, "
            f"{top[0]['direction']}, consistency={top[0]['consistency']})"
            if top else "No biomarkers found"
        ),
    }


def deg_direction_comparison(datasets: list, deg_datasets: dict = None,
                              mappings: dict = None,
                              comparisonA_groupA: str = None,
                              comparisonA_groupB: str = None,
                              comparisonB_groupA: str = None,
                              comparisonB_groupB: str = None,
                              adj_p_threshold: float = 0.05,
                              logfc_threshold: float = 0.5,
                              topN: int = 20, **_) -> dict:
    """
    Compare DE signatures between two biological comparisons.
    Finds concordant (same direction), discordant (opposite), and
    comparison-specific genes.

    Example: comparisonA = endometriosis vs normal
             comparisonB = adenomyosis vs normal
    """
    mappings = mappings or {}
    deg_datasets = deg_datasets or {}

    def collect_genes(gA: str, gB: str) -> dict:
        gene_logfcs: dict = {}
        for ds_name, ds in deg_datasets.items():
            for comp in ds["comparisons"]:
                a = resolve_group(comp["groupA"], mappings)
                b = resolve_group(comp["groupB"], mappings)
                if a == gA and b == gB:
                    direction = 1
                elif a == gB and b == gA:
                    direction = -1
                else:
                    continue
                df = comp["df"]
                sig = df[(df["adj_p"] < adj_p_threshold) & (df["logFC"].abs() > logfc_threshold)]
                for gene, row in sig.iterrows():
                    if gene not in gene_logfcs:
                        gene_logfcs[gene] = []
                    gene_logfcs[gene].append(float(row["logFC"]) * direction)
        return {gene: sum(lfc) / len(lfc) for gene, lfc in gene_logfcs.items()}

    genesA = collect_genes(comparisonA_groupA, comparisonA_groupB)
    genesB = collect_genes(comparisonB_groupA, comparisonB_groupB)

    if not genesA:
        return {"error": f"No genes found for {comparisonA_groupA} vs {comparisonA_groupB}"}
    if not genesB:
        return {"error": f"No genes found for {comparisonB_groupA} vs {comparisonB_groupB}"}

    shared = set(genesA) & set(genesB)
    only_A = set(genesA) - shared
    only_B = set(genesB) - shared

    concordant, discordant = [], []
    for gene in shared:
        lfc_a, lfc_b = genesA[gene], genesB[gene]
        entry = {
            "gene": gene,
            "logFC_A": round(lfc_a, 3),
            "logFC_B": round(lfc_b, 3),
            "direction_A": "UP" if lfc_a > 0 else "DOWN",
            "direction_B": "UP" if lfc_b > 0 else "DOWN",
        }
        (concordant if (lfc_a > 0) == (lfc_b > 0) else discordant).append(entry)

    concordant.sort(key=lambda x: -(abs(x["logFC_A"]) + abs(x["logFC_B"])) / 2)
    discordant.sort(key=lambda x: -(abs(x["logFC_A"]) + abs(x["logFC_B"])) / 2)

    only_A_top = sorted(
        [{"gene": g, "logFC_A": round(genesA[g], 3)} for g in only_A],
        key=lambda x: -abs(x["logFC_A"]),
    )[:topN]
    only_B_top = sorted(
        [{"gene": g, "logFC_B": round(genesB[g], 3)} for g in only_B],
        key=lambda x: -abs(x["logFC_B"]),
    )[:topN]

    label_A = f"{comparisonA_groupA} vs {comparisonA_groupB}"
    label_B = f"{comparisonB_groupA} vs {comparisonB_groupB}"

    return {
        "comparison_A": label_A,
        "comparison_B": label_B,
        "n_genes_A": len(genesA),
        "n_genes_B": len(genesB),
        "n_shared": len(shared),
        "n_concordant": len(concordant),
        "n_discordant": len(discordant),
        "n_specific_A": len(only_A),
        "n_specific_B": len(only_B),
        "top_concordant": concordant[:topN],
        "top_discordant": discordant[:topN],
        "specific_to_A": only_A_top,
        "specific_to_B": only_B_top,
        "interpretation": (
            f"{len(concordant)} concordant, {len(discordant)} discordant genes. "
            f"Top concordant: {concordant[0]['gene']} "
            f"(logFC_A={concordant[0]['logFC_A']}, logFC_B={concordant[0]['logFC_B']})"
            if concordant else "No shared genes between comparisons"
        ),
    }
