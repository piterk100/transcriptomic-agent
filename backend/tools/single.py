import logging
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu as _mwu
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score as _silhouette_score

logger = logging.getLogger(__name__)

_GENESET_CACHE: dict = {}


def _load_gene_sets() -> dict:
    """
    Load gene sets from local GMT file.
    Path must be set via GMT_FILE environment variable.
    Returns dict: {pathway_name: [gene1, gene2, ...]}
    Raises RuntimeError if GMT_FILE is not set or file does not exist.
    """
    global _GENESET_CACHE
    if _GENESET_CACHE:
        return _GENESET_CACHE

    gmt_path = os.environ.get("GMT_FILE", None)
    if not gmt_path:
        raise RuntimeError(
            "GMT_FILE environment variable not set. "
            "Download a GMT file from https://www.gsea-msigdb.org/gsea/msigdb/ "
            "and set: export GMT_FILE=/path/to/h.all.v2026.1.Hs.symbols.gmt"
        )
    # If path is relative, resolve it relative to this file's directory (backend/)
    if not os.path.isabs(gmt_path):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        # Go up one level if this file is in backend/tools/
        backend_dir = os.path.dirname(base_dir)
        gmt_path = os.path.join(backend_dir, gmt_path)
    if not os.path.exists(gmt_path):
        raise RuntimeError(
            f"GMT file not found: {gmt_path}\n"
            f"Make sure the file exists and GMT_FILE is set correctly in .env"
        )

    gene_sets = {}
    with open(gmt_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [g.upper() for g in parts[2:] if g]
            gene_sets[name] = genes

    _GENESET_CACHE = gene_sets
    return gene_sets


def _get_ds(datasets, name):
    return next((d for d in datasets if d["name"] == name), datasets[0])


def dataset_summary(datasets, **_):
    result = []
    for ds in datasets:
        expr: pd.DataFrame = ds["expr"]
        meta: pd.DataFrame = ds["meta"]
        group_col = ds["group_col"]
        groups = meta[group_col].value_counts().to_dict() if group_col in meta.columns else {}
        result.append({
            "name": ds["name"],
            "n_genes": len(expr.index),
            "n_samples": len(expr.columns),
            "group_col": group_col,
            "groups": {str(k): int(v) for k, v in groups.items()},
        })
    return result


def top_variable_genes(datasets, datasetName=None, n=40, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]
    variances = expr.var(axis=1).nlargest(n)
    return [{"gene": gene, "variance": round(float(var), 4)} for gene, var in variances.items()]


def differential_expression(datasets, datasetName=None, groupA=None, groupB=None, topN=25, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]
    meta: pd.DataFrame = ds["meta"]
    group_col = ds["group_col"]

    # Consistent with seeder: use index.intersection(expr.columns)
    sA = meta.index[meta[group_col] == groupA].intersection(expr.columns)
    sB = meta.index[meta[group_col] == groupB].intersection(expr.columns)

    if len(sA) < 3 or len(sB) < 3:
        return {"error": f"Too few samples: {groupA}({len(sA)}), {groupB}({len(sB)})"}

    nA, nB = len(sA), len(sB)

    # Per-gene MWU using scipy — consistent with seeder.generate_seeds()
    results = []
    for gene in expr.index:
        a_vals = expr.loc[gene, sA].values.astype(float)
        b_vals = expr.loc[gene, sB].values.astype(float)
        if np.std(a_vals) == 0 and np.std(b_vals) == 0:
            continue
        try:
            res = _mwu(a_vals, b_vals, alternative="two-sided")
        except Exception:
            continue
        logfc = float(a_vals.mean() - b_vals.mean())
        rbc = float(1.0 - 2.0 * res.statistic / (nA * nB))
        results.append({"gene": gene, "logFC": logfc, "rbc": rbc, "p": float(res.pvalue)})

    if not results:
        return {"error": f"No testable genes for {groupA} vs {groupB}"}

    df = pd.DataFrame(results)
    _, adj_p, _, _ = multipletests(df["p"].values, method="fdr_bh")
    df["adj_p"] = adj_p

    n_genes_tested = len(df)
    n_raw_sig = int((df["p"] < 0.05).sum())
    n_bh_sig = int((df["adj_p"] < 0.05).sum())
    logger.info(
        "DE %s | %s vs %s: tested=%d raw_sig=%d bh_sig=%d",
        ds["name"], groupA, groupB, n_genes_tested, n_raw_sig, n_bh_sig,
    )

    sig = df[df["adj_p"] < 0.05].sort_values("adj_p")
    top_up   = sig[sig["logFC"] > 0].head(topN)[["gene", "logFC", "rbc", "adj_p"]].round({"logFC": 3, "rbc": 3, "adj_p": 4}).to_dict("records")
    top_down = sig[sig["logFC"] < 0].head(topN)[["gene", "logFC", "rbc", "adj_p"]].round({"logFC": 3, "rbc": 3, "adj_p": 4}).to_dict("records")

    return {
        "dataset": ds["name"],
        "comparison": f"{groupA} vs {groupB}",
        "test": "Mann-Whitney U + BH (scipy, per-gene)",
        "n_genes_tested": n_genes_tested,
        "n_significant": n_bh_sig,
        "top_upregulated":   top_up,
        "top_downregulated": top_down,
    }


def gene_expression_by_group(datasets, datasetName=None, genes=None, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]
    meta: pd.DataFrame = ds["meta"]
    group_col = ds["group_col"]
    genes = [g for g in (genes or []) if g in expr.index]

    result = {}
    for gene in genes:
        result[gene] = {}
        for grp, idx in meta.groupby(group_col).groups.items():
            samples = [s for s in idx if s in expr.columns]
            vals = expr.loc[gene, samples].values if samples else np.array([])
            result[gene][str(grp)] = {
                "mean": round(float(vals.mean()), 3) if len(vals) else 0,
                "std":  round(float(vals.std()),  3) if len(vals) else 0,
                "n":    len(vals),
            }
    return {"dataset": ds["name"], "result": result}


def nonlinear_rule(datasets, datasetName=None, geneHigh=None, geneLow=None, targetGroup=None, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]
    meta: pd.DataFrame = ds["meta"]
    group_col = ds["group_col"]

    if geneHigh not in expr.index or geneLow not in expr.index:
        return {"error": f"Genes not found: {geneHigh}, {geneLow}"}

    thr_h = float(np.median(expr.loc[geneHigh]))
    thr_l = float(np.median(expr.loc[geneLow]))

    samples = expr.columns.tolist()
    pos = [s for s in samples if meta.loc[s, group_col] == targetGroup] if targetGroup in meta[group_col].values else []
    neg = [s for s in samples if s not in pos]
    match = [s for s in samples if expr.loc[geneHigh, s] > thr_h and expr.loc[geneLow, s] <= thr_l]

    tp = len([s for s in match if s in pos])
    fp = len([s for s in match if s in neg])
    fn = len([s for s in pos if s not in match])
    pr = tp / (tp + fp + 1e-9)
    re = tp / (tp + fn + 1e-9)
    f1 = 2 * pr * re / (pr + re + 1e-9)

    lin_pos = expr.loc[geneHigh, pos].values - expr.loc[geneLow, pos].values if pos else np.array([0])
    lin_neg = expr.loc[geneHigh, neg].values - expr.loc[geneLow, neg].values if neg else np.array([0])
    _, lin_p = stats.ttest_ind(lin_pos, lin_neg)

    return {
        "dataset": ds["name"],
        "rule": f"{geneHigh}>{thr_h:.2f} AND {geneLow}<={thr_l:.2f}",
        "predicts": targetGroup,
        "f1": round(f1, 3),
        "precision": round(pr, 3),
        "recall": round(re, 3),
        "linear_p": round(float(lin_p), 4),
    }


def contextual_modules(datasets, datasetName=None, contextGene=None, topN=20, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]

    if contextGene not in expr.index:
        return {"error": f"Gene not found: {contextGene}"}

    ctx_scores = expr.loc[contextGene].sort_values()
    half = len(ctx_scores) // 2
    ctx_low  = ctx_scores.index[:half].tolist()
    ctx_high = ctx_scores.index[half:].tolist()

    # Select top variable genes; exclude constant genes (std=0) to avoid NaN in corrcoef
    variances = expr.var(axis=1)
    top_genes = [g for g in variances.nlargest(topN).index if variances[g] > 0]

    def avg_corr(samples):
        sub = expr.loc[top_genes, samples].T
        if len(sub) < 3 or len(top_genes) < 2:
            return 0.0
        corr_mat = np.corrcoef(sub.values.T)
        n = len(top_genes)
        upper = corr_mat[np.triu_indices(n, k=1)]
        return float(np.nanmean(np.abs(upper)))

    c_high = avg_corr(ctx_high)
    c_low  = avg_corr(ctx_low)

    return {
        "dataset": ds["name"],
        "context_gene": contextGene,
        "corr_high": round(c_high, 3),
        "corr_low":  round(c_low, 3),
        "specificity": round(abs(c_high - c_low), 3),
        "interpretation": f"Co-expression stronger when {contextGene} is {'HIGH' if c_high > c_low else 'LOW'}",
    }


def pathway_enrichment(datasets: dict, genes: list = None, deg_datasets: dict = None,
                       deg_dataset_name: str = None, adj_p_threshold: float = 0.05,
                       logfc_threshold: float = 0.5, **_) -> dict:
    """
    Hypergeometric test + BH correction against gene sets from local GMT file.
    GMT file path must be set via GMT_FILE environment variable.
    If deg_dataset_name is provided, significant DE genes are extracted automatically
    from that DEG dataset (adj_p < adj_p_threshold and |logFC| > logfc_threshold).
    """
    from scipy.stats import hypergeom
    import pandas as pd

    # Extract genes from a DEG dataset if requested
    if deg_dataset_name and deg_datasets and deg_dataset_name in deg_datasets:
        ds = deg_datasets[deg_dataset_name]
        extracted = []
        for comp in ds["comparisons"]:
            df = comp["df"]
            sig = df[(df["adj_p"] < adj_p_threshold) & (df["logFC"].abs() > logfc_threshold)]
            extracted.extend(sig.index.tolist())
        if not genes:
            genes = list(set(extracted))

    try:
        all_gene_sets = _load_gene_sets()
    except RuntimeError as e:
        return {"error": str(e)}

    N = 25000  # background: human transcriptome
    gene_set_input = {g.upper() for g in (genes or [])}
    n = len(gene_set_input)

    if n == 0:
        return {"error": "No genes provided"}

    results = []
    for pathway, pw_genes in all_gene_sets.items():
        K = len(pw_genes)
        if K == 0:
            continue
        overlap = [g for g in pw_genes if g in gene_set_input]
        k = len(overlap)
        if k < 2:
            continue
        p = float(hypergeom.sf(k - 1, N, K, n))
        enrichment = round((k / n) / (K / N), 2) if n > 0 else 0
        results.append({
            "pathway": pathway,
            "overlap_genes": overlap,
            "k": k,
            "K": K,
            "n_input": n,
            "enrichment_fold": enrichment,
            "p": round(p, 6),
        })

    if not results:
        return {
            "genes_input": n,
            "n_pathways_tested": len(all_gene_sets),
            "top_enriched": [],
            "n_significant": 0,
            "interpretation": "No significant enrichment found (k<2 for all pathways)",
        }

    df = pd.DataFrame(results)
    _, adj_p, _, _ = multipletests(df["p"].values, method="fdr_bh")
    df["adj_p"] = adj_p.round(6)
    df = df.sort_values("adj_p")

    sig = df[df["adj_p"] < 0.05]
    top = df.head(10).to_dict("records")

    return {
        "genes_input": n,
        "n_pathways_tested": len(all_gene_sets),
        "n_significant": len(sig),
        "top_enriched": top,
        "interpretation": (
            f"Top: {top[0]['pathway']} "
            f"(k={top[0]['k']}, x{top[0]['enrichment_fold']}, "
            f"adj_p={top[0]['adj_p']:.4f})"
        ) if top else "No significant enrichment found",
    }


def batch_detection(datasets, datasetName=None, genes=None, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]
    meta: pd.DataFrame = ds["meta"]
    group_col = ds["group_col"]

    avail_genes = [g for g in (genes or []) if g in expr.index]
    if not avail_genes:
        return {"error": "None of the provided genes are available in the dataset"}

    axis_score = expr.loc[avail_genes].mean(axis=0)

    bio_kw = {"source", "condition", "characteristics", "tissue", "group", "disease"}
    candidate_cols = [c for c in meta.columns if c != group_col and not any(kw in c.lower() for kw in bio_kw)]

    raw_results = []
    raw_ps = []
    for col in candidate_cols[:10]:
        groups_vals = meta[col].dropna().unique()
        if len(groups_vals) < 2 or len(groups_vals) > 8:
            continue
        # Require each group has >= 2 samples (F-test assumption)
        groups_data = [axis_score[meta.index[meta[col] == v]].values for v in groups_vals]
        groups_data = [g for g in groups_data if len(g) >= 2]
        if len(groups_data) < 2:
            continue
        f_stat, p_val = stats.f_oneway(*groups_data)
        if np.isnan(f_stat):
            continue
        raw_results.append({
            "metadata_col": col,
            "n_groups": len(groups_data),
            "f_stat": round(float(f_stat), 3),
            "p": round(float(p_val), 4),
        })
        raw_ps.append(float(p_val))

    if raw_ps:
        _, adj_ps, _, _ = multipletests(raw_ps, method="fdr_bh")
        for i, r in enumerate(raw_results):
            r["adj_p"] = round(float(adj_ps[i]), 4)
            r["suspicious"] = bool(adj_ps[i] < 0.05)

    raw_results.sort(key=lambda r: r.get("adj_p", 1.0))
    suspicious = [r for r in raw_results if r.get("suspicious", False)]
    return {
        "dataset": ds["name"],
        "n_genes_tested": len(avail_genes),
        "n_columns_tested": len(raw_results),
        "batch_signals": raw_results[:5],
        "n_suspicious": len(suspicious),
        "interpretation": (
            f"WARNING: axis correlates with {', '.join(r['metadata_col'] for r in suspicious)} — possible batch artifact"
            if suspicious else "No suspicious batch signals detected (F-test BH-adjusted)"
        ),
    }


def subgroup_discovery(datasets, datasetName=None, group=None, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]
    meta: pd.DataFrame = ds["meta"]
    group_col = ds["group_col"]

    samples = meta.index[meta[group_col] == group].tolist()
    samples = [s for s in samples if s in expr.columns]
    if len(samples) < 4:
        return {"error": f"Too few samples in group '{group}': {len(samples)}"}

    top_genes = expr[samples].var(axis=1).nlargest(100).index.tolist()
    sub_expr = expr.loc[top_genes, samples].T  # samples x genes

    pca = PCA(n_components=min(10, len(samples) - 1))
    pca_coords = pca.fit_transform(sub_expr.values)

    km = KMeans(n_clusters=2, random_state=42, n_init=10)
    labels = km.fit_predict(pca_coords)

    # Silhouette score measures cluster cohesion/separation (-1 to 1; <0.2 = poorly separated)
    try:
        sil_score = round(float(_silhouette_score(pca_coords, labels)), 3)
    except Exception:
        sil_score = 0.0

    # Silhouette stability: 10 runs with different random seeds on full data
    sil_runs = []
    for seed in range(10):
        try:
            km_r = KMeans(n_clusters=2, random_state=seed, n_init=10)
            lbl_r = km_r.fit_predict(pca_coords)
            if len(np.unique(lbl_r)) == 2:
                sil_runs.append(float(_silhouette_score(pca_coords, lbl_r)))
        except Exception:
            pass
    silhouette_mean = round(float(np.mean(sil_runs)), 3) if sil_runs else sil_score
    silhouette_std = round(float(np.std(sil_runs)), 3) if sil_runs else 0.0
    silhouette_stable = bool(silhouette_std < 0.05)

    # Subsample stability: 10 runs on 80% random subsamples
    rng_sub = np.random.default_rng(0)
    sil_sub_runs = []
    n_sub = max(4, int(0.8 * len(samples)))
    for _ in range(10):
        idx = rng_sub.choice(len(samples), size=n_sub, replace=False)
        sub_coords = pca_coords[idx]
        try:
            km_s = KMeans(n_clusters=2, random_state=42, n_init=10)
            lbl_s = km_s.fit_predict(sub_coords)
            if len(np.unique(lbl_s)) == 2:
                sil_sub_runs.append(float(_silhouette_score(sub_coords, lbl_s)))
        except Exception:
            pass
    silhouette_subsample_mean = round(float(np.mean(sil_sub_runs)), 3) if sil_sub_runs else sil_score
    silhouette_subsample_std = round(float(np.std(sil_sub_runs)), 3) if sil_sub_runs else 0.0

    sub1 = [samples[i] for i, l in enumerate(labels) if l == 0]
    sub2 = [samples[i] for i, l in enumerate(labels) if l == 1]

    n_samples = len(samples)
    if min(len(sub1), len(sub2)) < max(3, int(0.15 * n_samples)):
        return {"error": f"Subgroup too small after clustering: sub1={len(sub1)}, sub2={len(sub2)}"}

    # MWU test per marker gene + BH correction
    raw_ps, raw_markers = [], []
    for gene in top_genes:
        v1 = expr.loc[gene, sub1].values
        v2 = expr.loc[gene, sub2].values
        log_fc = float(v1.mean() - v2.mean())
        try:
            _, p = stats.mannwhitneyu(v1, v2, alternative="two-sided")
        except ValueError:
            p = 1.0
        raw_ps.append(float(p))
        raw_markers.append({"gene": gene, "logFC": round(log_fc, 3), "p": round(float(p), 4)})

    _, adj_ps, _, _ = multipletests(raw_ps, method="fdr_bh")
    for i, m in enumerate(raw_markers):
        m["adj_p"] = round(float(adj_ps[i]), 4)

    raw_markers.sort(key=lambda m: (m["adj_p"], -abs(m["logFC"])))
    sig_markers = [m for m in raw_markers if m["adj_p"] < 0.05]

    return {
        "dataset": ds["name"],
        "group": group,
        "n_samples": len(samples),
        "subgroup_1": {"n": len(sub1), "samples": sub1},
        "subgroup_2": {"n": len(sub2), "samples": sub2},
        "silhouette_score": sil_score,
        "silhouette_mean": silhouette_mean,
        "silhouette_std": silhouette_std,
        "silhouette_stable": silhouette_stable,
        "silhouette_subsample_mean": silhouette_subsample_mean,
        "silhouette_subsample_std": silhouette_subsample_std,
        "n_significant_markers": len(sig_markers),
        "top_markers_sub1": [m for m in sig_markers if m["logFC"] > 0][:10],
        "top_markers_sub2": [m for m in sig_markers if m["logFC"] < 0][:10],
        "interpretation": (
            f"Detected 2 subgroups within '{group}' ({len(sub1)} vs {len(sub2)} samples, PCA+KMeans), "
            f"{len(sig_markers)} marker genes (MWU adj_p<0.05)"
            + ("" if sil_score >= 0.2 else "; subgroups may not be well-separated (silhouette<0.2)")
        ),
    }


def gene_network_hub(datasets, datasetName=None, topN=30, corrThreshold=0.6, use_permutation_threshold=True, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]

    # Exclude constant genes (var=0) which cause NaN in corrcoef
    variances = expr.var(axis=1)
    top_genes = [g for g in variances.nlargest(200).index if variances[g] > 0][:150]
    if len(top_genes) < 2:
        return {"error": "Too few variable genes for network analysis"}
    sub = expr.loc[top_genes].values  # genes x samples

    corr_mat = np.corrcoef(sub)  # genes x genes
    np.nan_to_num(corr_mat, nan=0.0, copy=False)  # guard against any residual NaN
    abs_corr = np.abs(corr_mat)
    np.fill_diagonal(abs_corr, 0)

    n_genes = len(top_genes)
    triu_i, triu_j = np.triu_indices(n_genes, k=1)
    obs_abs = abs_corr[triu_i, triu_j]

    if use_permutation_threshold:
        # Permutation test: shuffle sample labels 200 times, count how often the null
        # |corr| >= observed |corr| for each gene pair, then FDR-correct across all pairs.
        N_PERM = 200
        rng_perm = np.random.default_rng(42)
        n_samples = sub.shape[1]
        exceed_count = np.zeros(len(triu_i), dtype=float)

        for _ in range(N_PERM):
            perm_idx = rng_perm.permutation(n_samples)
            perm_corr = np.corrcoef(sub[:, perm_idx])
            np.nan_to_num(perm_corr, nan=0.0, copy=False)
            exceed_count += np.abs(perm_corr[triu_i, triu_j]) >= obs_abs

        perm_ps = exceed_count / N_PERM
        _, adj_perm_ps, _, _ = multipletests(perm_ps, method="fdr_bh")

        # Map (i, j) → FDR-corrected permutation p-value for edge annotation
        perm_p_lookup = {
            (int(triu_i[k]), int(triu_j[k])): float(adj_perm_ps[k])
            for k in range(len(triu_i))
        }

        # Build adjacency mask from significant edges (FDR < 0.05)
        mask = np.zeros((n_genes, n_genes), dtype=bool)
        for k in np.where(adj_perm_ps < 0.05)[0]:
            i, j = int(triu_i[k]), int(triu_j[k])
            mask[i, j] = True
            mask[j, i] = True
    else:
        perm_p_lookup = {}
        mask = abs_corr >= corrThreshold

    connectivity = mask.sum(axis=1)
    top_idx = np.argsort(-connectivity)[:topN]
    top_hubs = [{"gene": top_genes[i], "degree": int(connectivity[i])} for i in top_idx]

    edge_idx = np.argwhere(mask)
    edge_idx = edge_idx[edge_idx[:, 0] < edge_idx[:, 1]]  # upper triangle
    edges = []
    for i, j in edge_idx:
        edge = {"g1": top_genes[i], "g2": top_genes[j], "r": round(float(abs_corr[i, j]), 3)}
        if perm_p_lookup:
            edge["perm_p"] = round(perm_p_lookup.get((i, j), 1.0), 4)
        edges.append(edge)
    edges.sort(key=lambda e: -e["r"])

    return {
        "dataset": ds["name"],
        "corr_threshold": corrThreshold,
        "use_permutation_threshold": use_permutation_threshold,
        "n_edges": len(edges),
        "n_genes_tested": len(top_genes),
        "top_hubs": top_hubs,
        "top_edges": edges[:20],
        "interpretation": f"Top hub: {top_hubs[0]['gene']} ({top_hubs[0]['degree']} connections)" if top_hubs else "No hubs found at this threshold",
    }
