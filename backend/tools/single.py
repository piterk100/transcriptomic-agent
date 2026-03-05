import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata as _rankdata, norm as _norm
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score as _silhouette_score

PATHWAYS = {
    "Hallmark_Hypoxia":       ["VEGFA","SLC2A1","LDHA","ENO1","HK2","ALDOA","PGAM1","TPI1","PGK1","PFKL","GAPDH","PKM","BNIP3","BNIP3L","DDIT4","P4HA1","P4HA2","PLOD1","PLOD2","COL5A1"],
    "Hallmark_EMT":           ["VIM","FN1","CDH2","TWIST1","SNAI1","SNAI2","ZEB1","ZEB2","MMP2","MMP9","ITGA5","TGFB1","ACTA2","COL1A1","COL1A2","COL3A1","SPARC","THBS1","SERPINE1","LOXL2"],
    "Hallmark_Inflammation":  ["IL6","IL1B","TNF","CXCL8","CXCL1","CXCL2","CCL2","CCL5","ICAM1","VCAM1","PTGS2","NOS2","IL6ST","NFKBIA","RELB","IRF1","STAT1","CXCL10","MX1","ISG15"],
    "Hallmark_Proliferation": ["MKI67","PCNA","MCM2","MCM6","CCNA2","CCNB1","CCNE1","CDK1","CDK2","BUB1","BUB1B","TOP2A","AURKB","AURKA","PLK1","CDC20","CDKN2A","RRM2","TYMS","E2F1"],
    "Hallmark_Apoptosis":     ["CASP3","CASP7","CASP8","CASP9","BAX","BAD","BCL2","BCL2L1","CYCS","APAF1","TP53","PUMA","NOXA","BID","FAS","FASLG","TNFRSF10A","TNFRSF10B","DIABLO","XIAP"],
    "Hallmark_IFN_gamma":     ["STAT1","IRF1","CXCL10","CXCL9","GBP1","GBP2","MX1","OAS1","IFIT1","IFIT2","IFIT3","ISG15","HLA-A","HLA-B","HLA-C","B2M","TAP1","TAPBP","PSMB9","PSMB10"],
    "Hallmark_Angiogenesis":  ["VEGFA","VEGFB","VEGFC","KDR","FLT1","ANGPT1","ANGPT2","TEK","PDGFA","PDGFB","FGF2","FGFR1","NRP1","CXCL12","ENG","PECAM1","CDH5","NOTCH1","DLL4","HIF1A"],
    "KEGG_ECM_Receptor":      ["COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","FN1","LAMA1","LAMB1","LAMC1","ITGA1","ITGA2","ITGA3","ITGA5","ITGB1","ITGB3","ITGB5","THBS1","THBS2","VTN","SPP1"],
    "KEGG_JAK_STAT":          ["JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","IL6","IL2","IL4","IL12A","IFNA1","IFNG","EPO","GH1","CISH"],
    "KEGG_TGF_beta":          ["TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7","ACVR1","ACVR2A","BMP2","BMP4","BMP7","ID1","ID2","ID3","INHBA","LTBP1","NOG"],
    "Hallmark_MYC_targets":   ["MYC","MYCN","MAX","CDK4","CDK6","RB1","E2F1","E2F2","CCND1","CCND2","CCNE1","MCM4","MCM5","CDC25A","LDHA","FASN","ODC1","SRM","NME1","NCL"],
    "Hallmark_P53_pathway":   ["TP53","CDKN1A","MDM2","BAX","PUMA","NOXA","GADD45A","GADD45B","SESN1","SESN2","RRM2B","TIGAR","DDB2","XPC","POLK","ZMAT3","FAS","TNFRSF10B","PERP","APAF1"],
}


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


def _mwu_vectorized(vals_a: np.ndarray, vals_b: np.ndarray):
    """
    Vectorized Mann-Whitney U test across all genes simultaneously.
    vals_a: (n_genes, nA), vals_b: (n_genes, nB)
    Returns U-statistics and two-sided p-values (normal approximation with tie correction).
    """
    nA, nB = vals_a.shape[1], vals_b.shape[1]
    n = nA + nB
    combined = np.concatenate([vals_a, vals_b], axis=1)  # (n_genes, n)

    # Rank each gene's values with average tie handling
    ranks = np.vstack([_rankdata(row, method="average") for row in combined])
    rank_sum_A = ranks[:, :nA].sum(axis=1)
    U = rank_sum_A - nA * (nA + 1) / 2

    # Tie correction: sum of (t^3 - t) for each tied group per gene
    tie_sum = np.array([
        np.sum((_cnts := np.unique(row, return_counts=True)[1]) ** 3 - _cnts)
        for row in combined
    ], dtype=float)
    sigma_sq = nA * nB / (n * (n - 1)) * ((n ** 3 - n) / 12 - tie_sum / 12)
    sigma = np.sqrt(np.maximum(sigma_sq, 1e-10))

    z = (U - nA * nB / 2) / sigma
    p_vals = 2 * _norm.sf(np.abs(z))
    return U, p_vals


def differential_expression(datasets, datasetName=None, groupA=None, groupB=None, topN=25, **_):
    ds = _get_ds(datasets, datasetName)
    expr: pd.DataFrame = ds["expr"]
    meta: pd.DataFrame = ds["meta"]
    group_col = ds["group_col"]

    sA = meta.index[meta[group_col] == groupA].tolist()
    sB = meta.index[meta[group_col] == groupB].tolist()
    sA = [s for s in sA if s in expr.columns]
    sB = [s for s in sB if s in expr.columns]

    if len(sA) < 2 or len(sB) < 2:
        return {"error": f"Too few samples: {groupA}({len(sA)}), {groupB}({len(sB)})"}

    nA, nB = len(sA), len(sB)
    exprA_vals = expr[sA].values
    exprB_vals = expr[sB].values

    U, p_vals = _mwu_vectorized(exprA_vals, exprB_vals)

    # Effect sizes: logFC (assumes log-normalized input) + rank-biserial correlation
    log_fc = expr[sA].mean(axis=1).values - expr[sB].mean(axis=1).values
    rbc = 1.0 - 2.0 * U / (nA * nB)  # rank-biserial correlation ∈ [-1, 1]

    _, adj_p, _, _ = multipletests(np.nan_to_num(p_vals, nan=1.0), method="fdr_bh")

    results_df = pd.DataFrame({
        "gene": expr.index,
        "logFC": np.round(log_fc, 3),
        "rbc":   np.round(rbc, 3),
        "p":     p_vals,
        "adj_p": np.round(adj_p, 4),
    }).set_index("gene")

    sig = results_df[results_df["adj_p"] < 0.05]
    top_up   = sig.nlargest(topN, "logFC")[["logFC", "rbc", "adj_p"]].reset_index().to_dict("records")
    top_down = sig.nsmallest(topN, "logFC")[["logFC", "rbc", "adj_p"]].reset_index().to_dict("records")

    return {
        "dataset": ds["name"],
        "comparison": f"{groupA} vs {groupB}",
        "test": "Mann-Whitney U (BH-adjusted, tie-corrected)",
        "n_significant": int(len(sig)),
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


def pathway_enrichment(datasets, genes=None, **_):
    from scipy.stats import hypergeom

    gene_set = {g.upper() for g in (genes or [])}
    N_BACKGROUND = 20_000  # assumed background transcriptome size
    n = len(gene_set)
    raw_ps, results = [], []

    for pathway, pw_genes in PATHWAYS.items():
        overlap = [g for g in pw_genes if g.upper() in gene_set]
        k, K = len(overlap), len(pw_genes)
        if k == 0:
            continue
        # P(X >= k) under hypergeometric null: drawing n genes from N_BACKGROUND that contains K pathway genes
        p_val = float(hypergeom.sf(k - 1, N_BACKGROUND, K, n))
        expected = max(n * K / N_BACKGROUND, 0.01)
        results.append({
            "pathway": pathway,
            "overlap_genes": overlap,
            "k": k, "K": K,
            "enrichment_fold": round(k / expected, 2),
            "p": round(p_val, 6),
        })
        raw_ps.append(p_val)

    if results:
        _, adj_ps, _, _ = multipletests(raw_ps, method="fdr_bh")
        for i, r in enumerate(results):
            r["adj_p"] = round(float(adj_ps[i]), 6)
        results.sort(key=lambda r: (r["adj_p"], -r["enrichment_fold"]))

    sig = [r for r in results if r.get("adj_p", 1.0) < 0.05]
    return {
        "genes_input": n,
        "n_pathways_tested": len(results),
        "n_significant": len(sig),
        "top_enriched": results[:8],
        "interpretation": (
            f"Top enriched: {results[0]['pathway']} "
            f"(k={results[0]['k']}, x{results[0]['enrichment_fold']}, adj_p={results[0]['adj_p']})"
            if results else "No significant enrichment"
        ),
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

    sub1 = [samples[i] for i, l in enumerate(labels) if l == 0]
    sub2 = [samples[i] for i, l in enumerate(labels) if l == 1]

    if len(sub1) < 2 or len(sub2) < 2:
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
