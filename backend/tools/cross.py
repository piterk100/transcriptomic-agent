import numpy as np
from scipy import stats
from scipy.stats import rankdata as _rankdata, norm as _norm, chi2 as _chi2
from statsmodels.stats.multitest import multipletests


def resolve_group(group_name: str, mappings: dict) -> str:
    """
    Given a group name and a mappings dict, return the canonical name
    if this group_name is listed as an alias. Otherwise return group_name.
    """
    for canonical, aliases in mappings.items():
        if group_name in aliases or group_name == canonical:
            return canonical
    return group_name


def _cohen_d(a, b):
    pooled_sd = np.sqrt((np.std(a, ddof=1) ** 2 + np.std(b, ddof=1) ** 2) / 2) + 1e-9
    return (np.mean(a) - np.mean(b)) / pooled_sd


def _mwu_cross(expr_a_vals: np.ndarray, expr_b_vals: np.ndarray):
    """Vectorized MWU with tie correction; returns (U, p_vals)."""
    nA, nB = expr_a_vals.shape[1], expr_b_vals.shape[1]
    n = nA + nB
    combined = np.concatenate([expr_a_vals, expr_b_vals], axis=1)
    ranks = np.vstack([_rankdata(row, method="average") for row in combined])
    rank_sum_A = ranks[:, :nA].sum(axis=1)
    U = rank_sum_A - nA * (nA + 1) / 2
    tie_sum = np.array([
        np.sum((_cnts := np.unique(row, return_counts=True)[1]) ** 3 - _cnts)
        for row in combined
    ], dtype=float)
    sigma = np.sqrt(np.maximum(
        nA * nB / (n * (n - 1)) * ((n ** 3 - n) / 12 - tie_sum / 12), 1e-10
    ))
    p_vals = 2 * _norm.sf(np.abs((U - nA * nB / 2) / sigma))
    return U, p_vals


def cross_dataset_de(datasets, groupA=None, groupB=None, mappings=None, topN=30, min_effect_size=0.5, **_):
    """
    Cross-dataset DE using Mann-Whitney U with:
    - per-dataset BH correction (controls within-dataset FDR)
    - Fisher's combined test across datasets for directionally consistent genes
    - global BH correction on Fisher's combined p-values
    - effect size filter: genes with avg |logFC| < min_effect_size are excluded after
      directional consistency check (reported in n_filtered_by_effect_size)
    """
    mappings = mappings or {}
    per_ds: dict[str, dict] = {}  # gene → {ds_name → {logFC, rbc, p, adj_p}}
    used = []

    for ds in datasets:
        expr = ds["expr"]
        meta = ds["meta"]
        group_col = ds["group_col"]

        resolved = meta[group_col].apply(lambda g: resolve_group(g, mappings))
        sA = [s for s in meta.index[resolved == groupA] if s in expr.columns]
        sB = [s for s in meta.index[resolved == groupB] if s in expr.columns]
        if len(sA) < 2 or len(sB) < 2:
            continue
        used.append(ds["name"])

        nA, nB = len(sA), len(sB)
        U, p_raw = _mwu_cross(expr[sA].values, expr[sB].values)
        log_fc = expr[sA].mean(axis=1).values - expr[sB].mean(axis=1).values
        rbc = 1.0 - 2.0 * U / (nA * nB)

        _, adj_p, _, _ = multipletests(np.nan_to_num(p_raw, nan=1.0), method="fdr_bh")

        for i, gene in enumerate(expr.index):
            if gene not in per_ds:
                per_ds[gene] = {}
            per_ds[gene][ds["name"]] = {
                "logFC": float(log_fc[i]),
                "rbc":   float(rbc[i]),
                "p":     float(p_raw[i]) if not np.isnan(p_raw[i]) else 1.0,
                "adj_p": float(adj_p[i]),
            }

    if not used:
        return {"error": f"No datasets found with groups {groupA}, {groupB}"}

    common = [g for g, dm in per_ds.items() if len(dm) == len(used)]

    # Build candidates: directionally consistent + minimum effect size across all datasets
    n_directionally_inconsistent = 0
    n_below_effect_size = 0
    fisher_genes, fisher_ps, entries = [], [], []
    for gene in common:
        dm = per_ds[gene]
        lfcs = [dm[d]["logFC"] for d in used]
        if not (all(x > 0 for x in lfcs) or all(x < 0 for x in lfcs)):
            n_directionally_inconsistent += 1
            continue
        avg_abs_lfc = float(np.mean(np.abs(lfcs)))
        if avg_abs_lfc < min_effect_size:
            n_below_effect_size += 1
            continue
        # Fisher's combined p-value (meta-analysis of per-dataset raw p-values)
        ps_clipped = np.clip([dm[d]["p"] for d in used], 1e-300, 1.0)
        chi2_stat = -2.0 * np.sum(np.log(ps_clipped))
        fisher_p = float(_chi2.sf(chi2_stat, df=2 * len(used)))
        fisher_genes.append(gene)
        fisher_ps.append(fisher_p)
        entries.append({
            "gene": gene,
            "direction": "UP" if lfcs[0] > 0 else "DOWN",
            "avg_abs_logFC": round(avg_abs_lfc, 3),
            "fisher_p": fisher_p,
            "n_sig_datasets": sum(1 for d in used if dm[d]["adj_p"] < 0.05),
            "n_datasets": len(used),
            "per_dataset": [{
                "dataset": d,
                "logFC": round(dm[d]["logFC"], 3),
                "rbc":   round(dm[d]["rbc"], 3),
                "adj_p": round(dm[d]["adj_p"], 4),
            } for d in used],
        })

    # Global BH correction on Fisher's combined p-values
    if fisher_ps:
        _, adj_fisher, _, _ = multipletests(fisher_ps, method="fdr_bh")
        for i, entry in enumerate(entries):
            entry["fisher_p"] = round(fisher_ps[i], 6)
            entry["fisher_adj_p"] = round(float(adj_fisher[i]), 6)
        entries.sort(key=lambda x: (x["fisher_adj_p"], -x["avg_abs_logFC"]))

    return {
        "datasets_used": used,
        "comparison": f"{groupA} vs {groupB}",
        "test": "Mann-Whitney U (per-dataset BH) + Fisher meta-analysis (global BH)",
        "n_consistent_genes": len(entries),
        "n_common_genes_tested": len(common),
        "n_filtered_by_directionality": n_directionally_inconsistent,
        "n_filtered_by_effect_size": n_below_effect_size,
        "min_effect_size_used": min_effect_size,
        "top_consistent_up":   [x for x in entries if x["direction"] == "UP"][:topN],
        "top_consistent_down": [x for x in entries if x["direction"] == "DOWN"][:topN],
    }


def cross_dataset_correlation(datasets, genes=None, **_):
    genes = genes or []
    per = []
    for ds in datasets:
        expr = ds["expr"]
        avail = [g for g in genes if g in expr.index]
        if len(avail) < 2:
            per.append({"dataset": ds["name"], "error": "Too few genes"})
            continue
        sub = expr.loc[avail]
        # Filter constant genes — corrcoef returns NaN for zero-variance rows
        variable = sub.var(axis=1) > 0
        sub = sub[variable].values
        if sub.shape[0] < 2:
            per.append({"dataset": ds["name"], "error": "Too few variable genes"})
            continue
        # Spearman rank-based correlation — robust to non-normality
        ranks = np.apply_along_axis(lambda x: _rankdata(x, method="average"), 1, sub)
        corr_mat = np.corrcoef(ranks)
        n = sub.shape[0]
        upper = np.abs(corr_mat[np.triu_indices(n, k=1)])
        upper = upper[~np.isnan(upper)]
        if len(upper) == 0:
            per.append({"dataset": ds["name"], "error": "No gene pairs available for correlation"})
            continue
        per.append({
            "dataset": ds["name"],
            "avg_module_corr": round(float(upper.mean()), 3),
            "n_genes": int(variable.sum()),
        })

    cs = [r["avg_module_corr"] for r in per if "avg_module_corr" in r]
    min_c = float(min(cs)) if cs else 0.0
    mean_c = float(np.mean(cs)) if cs else 0.0
    return {
        "genes_tested": genes,
        "per_dataset": per,
        "mean_replication": round(mean_c, 3),
        "min_replication": round(min_c, 3),
        "interpretation": "Module replicates strongly across datasets" if min_c > 0.5
            else "Module partially replicates" if min_c > 0.3
            else "Module does not replicate — possible artifact",
    }


def invariant_axis(datasets, groupA=None, groupB=None, mappings=None, topN=20, **_):
    all_genes = [set(ds["expr"].index) for ds in datasets]
    common = list(all_genes[0].intersection(*all_genes[1:])) if len(all_genes) > 1 else list(all_genes[0]) if all_genes else []
    if not common:
        return {"error": "No common genes across datasets"}

    mappings = mappings or {}
    valid_ds = []
    for ds in datasets:
        expr = ds["expr"]
        meta = ds["meta"]
        group_col = ds["group_col"]
        resolved = meta[group_col].apply(lambda g: resolve_group(g, mappings))
        sA = [s for s in meta.index[resolved == groupA] if s in expr.columns]
        sB = [s for s in meta.index[resolved == groupB] if s in expr.columns]
        if len(sA) >= 2 and len(sB) >= 2:
            valid_ds.append((ds, sA, sB))

    if len(valid_ds) < 2:
        return {"error": f"Need ≥2 datasets with groups {groupA} and {groupB}"}

    # Warn when any dataset has fewer than 10 samples in either group; reduce bootstrap N
    # to avoid overconfident estimates from permuting tiny groups.
    low_sample_warning = [
        ds["name"] for ds, sA, sB in valid_ds if min(len(sA), len(sB)) < 10
    ]
    N_BOOT = 100 if low_sample_warning else 500
    low_sample_warning_message = (
        f"Warning: datasets {low_sample_warning} have fewer than 10 samples in at least one group. "
        "Results may be unreliable."
        if low_sample_warning else ""
    )

    gene_scores = []
    for gene in common:
        effects = []
        for ds, sA, sB in valid_ds:
            expr = ds["expr"]
            d = _cohen_d(expr.loc[gene, sA].values, expr.loc[gene, sB].values)
            effects.append({"dataset": ds["name"], "d": float(d)})
        ds_arr = [e["d"] for e in effects]
        m = float(np.mean(np.abs(ds_arr)))
        cv = float(np.std(ds_arr)) / (m + 1e-9)
        consistent = all(d > 0 for d in ds_arr) or all(d < 0 for d in ds_arr)
        if not consistent:
            continue
        gene_scores.append({
            "gene": gene,
            "invariance_score": round(m / (1 + cv), 4),
            "mean_abs_effect": round(m, 3),
            "direction_consistent": True,
            "per_dataset": effects,
        })

    gene_scores.sort(key=lambda x: -x["invariance_score"])
    top_genes = [x["gene"] for x in gene_scores[:10]]

    # Bootstrap permutation test: for each top-5 gene, permute group labels N_BOOT times
    # to obtain an empirical null distribution of invariance scores.
    # N_BOOT is set to 100 (instead of 500) when any dataset has low sample counts (<5 per group)
    # to avoid overconfident p-values from tiny permutation spaces.
    rng = np.random.default_rng(42)

    for entry in gene_scores[:min(5, len(gene_scores))]:
        gene = entry["gene"]
        obs = entry["invariance_score"]
        boot_scores = []

        for _ in range(N_BOOT):
            boot_effects = []
            for ds, sA, sB in valid_ds:
                if gene not in ds["expr"].index:
                    continue
                all_s = sA + sB
                perm = rng.permutation(len(all_s))
                pA = [all_s[j] for j in perm[:len(sA)]]
                pB = [all_s[j] for j in perm[len(sA):]]
                d_b = _cohen_d(ds["expr"].loc[gene, pA].values, ds["expr"].loc[gene, pB].values)
                boot_effects.append(d_b)
            if len(boot_effects) == len(valid_ds):
                m_b = float(np.mean(np.abs(boot_effects)))
                cv_b = float(np.std(boot_effects)) / (m_b + 1e-9)
                boot_scores.append(m_b / (1 + cv_b))

        p_boot = float(np.mean(np.array(boot_scores) >= obs)) if boot_scores else 1.0
        entry["p_bootstrap"] = round(p_boot, 4)
        # Bonferroni-corrected threshold (testing up to 5 genes simultaneously)
        entry["bootstrap_significant"] = bool(p_boot < 0.01)

    # Composite axis test: Mann-Whitney U on mean expression of top invariant genes
    axis_per_ds = []
    for ds, sA, sB in valid_ds:
        expr = ds["expr"]
        avail = [g for g in top_genes if g in expr.index]
        ax_a = expr.loc[avail, sA].mean(axis=0).values
        ax_b = expr.loc[avail, sB].mean(axis=0).values
        u, p = stats.mannwhitneyu(ax_a, ax_b, alternative="two-sided")
        axis_per_ds.append({
            "dataset": ds["name"],
            "mwu_p": round(float(p), 4),
            "significant": bool(p < 0.05),
        })

    n_sig = sum(1 for x in axis_per_ds if x["significant"])
    n_boot_sig = sum(1 for e in gene_scores[:5] if e.get("bootstrap_significant", False))
    return {
        "comparison": f"{groupA} vs {groupB}",
        "n_datasets": len(valid_ds),
        "low_sample_warning": low_sample_warning,
        "low_sample_warning_message": low_sample_warning_message,
        "axis_consistent_in": f"{n_sig}/{len(valid_ds)} datasets",
        "n_bootstrap_validated": f"{n_boot_sig}/5 top genes (p_boot<0.01, Bonferroni)",
        "top_invariant_genes": gene_scores[:topN],
        "axis_per_dataset": axis_per_ds,
        "interpretation": (
            "Axis fully invariant and bootstrap-validated — strong biological signal"
            if n_sig == len(valid_ds) and n_boot_sig >= 3
            else "Axis fully invariant, weak bootstrap validation — check for batch effects"
            if n_sig == len(valid_ds)
            else "Axis partially invariant" if n_sig > len(valid_ds) / 2
            else "Axis inconsistent — possible artifact"
        ),
    }


def cross_dataset_rewiring(datasets, gene1=None, gene2=None, **_):
    """
    Uses Spearman correlation (rank-based, robust to non-normality) with p-value.
    Requires n >= 5 samples per dataset for reliable correlation estimates.
    """
    MIN_SAMPLES = 5
    per = []
    for ds in datasets:
        expr = ds["expr"]
        n = len(expr.columns)
        if gene1 not in expr.index or gene2 not in expr.index:
            per.append({"dataset": ds["name"], "error": "Gene(s) not found"})
            continue
        if n < MIN_SAMPLES:
            per.append({"dataset": ds["name"], "error": f"Too few samples: {n} < {MIN_SAMPLES}"})
            continue
        # Spearman correlation is rank-based — appropriate for non-normal expression data
        r, p = stats.spearmanr(expr.loc[gene1].values, expr.loc[gene2].values)
        per.append({
            "dataset": ds["name"],
            "spearman_r": round(float(r), 3),
            "p": round(float(p), 4),
            "n_samples": n,
        })

    valid = [x for x in per if "spearman_r" in x]
    if not valid:
        return {"error": f"No data available for {gene1}, {gene2}"}

    cs = [x["spearman_r"] for x in valid]
    rew = float(max(cs) - min(cs))
    return {
        "gene_pair": f"{gene1} — {gene2}",
        "per_dataset": per,
        "max_r": round(max(cs), 3),
        "min_r": round(min(cs), 3),
        "rewiring_magnitude": round(rew, 3),
        "interpretation": "Strong correlation rewiring across datasets" if rew > 0.6
            else "Moderate rewiring" if rew > 0.35
            else "Stable correlation across datasets",
    }
