from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_samples
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.exceptions import ConvergenceWarning
from scipy.cluster.hierarchy import leaves_list
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

import warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)



# ---------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------


@dataclass
class NMFProgram:
    """Single NMF factorization for one sample and one k."""
    w: np.ndarray           # genes x components
    h: np.ndarray           # components x cells
    genes: np.ndarray       # gene names
    cells: np.ndarray       # cell names
    d: Optional[np.ndarray] = None
    tol: Optional[float] = None
    n_iter: Optional[int] = None


@dataclass
class MultiNMFResult:
    """Analog of GeneNMF::multiNMF output (simplified)."""
    programs: Dict[str, NMFProgram]  # keys like "<sample>.k5"
    hvg_genes: List[str]
    samples: List[str]
    ks: List[int]


@dataclass
class MetaProgramResult:
    """Analog of GeneNMF::getMetaPrograms output (简化版)."""
    metaprograms_genes: Dict[str, List[str]]
    metaprograms_gene_weights: Dict[str, pd.Series]
    metaprograms_metrics: pd.DataFrame
    metaprograms_composition: pd.DataFrame
    programs_similarity: pd.DataFrame
    programs_clusters: pd.Series       # index: program_id, value: MP label
    linkage: np.ndarray                # scipy linkage matrix


# ---------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------


def _to_dense(X):
    if sp.issparse(X):
        return X.toarray()
    return np.asarray(X)


def _norm_vector(v: np.ndarray) -> np.ndarray:
    s = v.sum()
    if s > 0:
        return v / s
    return v


def weight_cumul(weights: pd.Series, weight_explained: float = 0.5) -> pd.Series:
    """
    近似 R 版 weightCumul:
    - 按权重降序
    - 累积和做一次归一化，再把最大值归一到 1
    - 取 cumulative < weight_explained 的基因
    """
    if not isinstance(weights, pd.Series):
        weights = pd.Series(weights)
    w_sorted = weights.sort_values(ascending=False)
    cs = w_sorted.cumsum()
    total = cs.sum()
    if total > 0:
        cs_norm = cs / total
        if cs_norm.max() > 0:
            cs_norm = cs_norm / cs_norm.max()
    else:
        cs_norm = cs
    selected = cs_norm < weight_explained
    return w_sorted[selected]


# ---------------------------------------------------------------------
# 对应 getDataMatrix
# ---------------------------------------------------------------------


def get_data_matrix(
    adata: ad.AnnData,
    layer: Optional[str] = None,
    use_raw: bool = False,
    genes: Optional[Sequence[str]] = None,
    center: bool = False,
    scale: bool = False,
    non_negative: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    从 AnnData 提取表达矩阵，类似 R 的 getDataMatrix：
    - 默认使用 adata.X (cells x genes)
    - 可选 layer / raw
    - 可选按 gene 子集
    - 可选 center / scale，然后裁掉负值以适配 NMF
    返回:
      X      : cells x genes (float ndarray)
      genes  : gene 名
      cells  : cell 名
    """
    if layer is not None:
        X = adata.layers[layer]
    elif use_raw and adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X

    var_names = adata.var_names

    if genes is not None:
        genes = [g for g in genes if g in var_names]
        if len(genes) == 0:
            raise ValueError("No genes left after subsetting in get_data_matrix().")
        idx = var_names.get_indexer(genes)
        X = X[:, idx]
        var_names = pd.Index(genes)

    X = _to_dense(X).astype(float)

    if center or scale:
        mean = X.mean(axis=0) if center else np.zeros(X.shape[1], dtype=float)
        std = X.std(axis=0, ddof=0) if scale else np.ones(X.shape[1], dtype=float)
        std[std == 0] = 1.0
        X = (X - mean) / std

    if non_negative:
        X[X < 0] = 0.0

    genes_out = np.array(var_names)
    cells_out = np.array(adata.obs_names)
    return X, genes_out, cells_out


# ---------------------------------------------------------------------
# 对应 findVariableFeatures_wfilters / findHVG
# ---------------------------------------------------------------------


def find_variable_features_wfilters(
    adata: ad.AnnData,
    nfeatures: int = 2000,
    genes_blocklist: Optional[Sequence[str]] = None,
    min_exp: float = 0.01,
    max_exp: float = 3.0,
    layer: Optional[str] = None,
) -> ad.AnnData:
    """
    Python 版 findVariableFeatures_wfilters：
    - 用 scanpy.pp.highly_variable_genes 先取 10000 HVG
    - 去掉 blocklist 里的基因
    - 再按平均表达过滤 (min_exp, max_exp)
    - 最后只保留 nfeatures 个 HVG
    """
    adata = adata.copy()
    n_top = min(adata.n_vars, 10000)

    sc.pp.filter_genes(adata, min_counts=5)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top,
        flavor="seurat_v3",
        layer=layer,
        # inplace=True,
    )

    hvgs = adata.var[adata.var["highly_variable"]].copy()

    if genes_blocklist is not None:
        genes_blocklist = set(genes_blocklist)
        hvgs = hvgs[~hvgs.index.isin(genes_blocklist)]

    # 平均表达
    if "means" not in hvgs.columns:
        X = adata.layers[layer] if layer is not None else adata.X
        means_all = np.asarray(X.mean(axis=0)).ravel()
        hvgs = hvgs.assign(
            means=means_all[adata.var_names.get_indexer(hvgs.index)]
        )

    hvgs = hvgs[(hvgs["means"] >= min_exp) & (hvgs["means"] <= max_exp)]

    # 选前 nfeatures
    if "highly_variable_rank" in hvgs.columns:
        hvgs = hvgs.sort_values("highly_variable_rank")
    elif "variances_norm" in hvgs.columns:
        hvgs = hvgs.sort_values("variances_norm", ascending=False)

    selected = hvgs.index[:nfeatures]

    adata.var["highly_variable"] = False
    adata.var.loc[selected, "highly_variable"] = True

    return adata


def find_hvg(
    adatas: Sequence[ad.AnnData],
    nfeatures: int = 2000,
    min_exp: float = 0.01,
    max_exp: float = 3.0,
    hvg_blocklist: Optional[Sequence[str]] = None,
    layer: Optional[str] = None,
) -> List[str]:
    """
    近似 R 版 findHVG:
    - 每个样本用上面的函数取 ncalc = min(5*nfeatures, n_genes) 个 HVG
    - 统计每个基因在多少个样本里是 HVG
    - 按出现次数排序，取前 nfeatures 作为全局 HVG
    """
    hvgs_list: List[set] = []
    for adata in adatas:
        ncalc = min(5 * nfeatures, adata.n_vars)
        tmp = find_variable_features_wfilters(
            adata,
            nfeatures=ncalc,
            genes_blocklist=hvg_blocklist,
            min_exp=min_exp,
            max_exp=max_exp,
            layer=layer,
        )
        hvgs = tmp.var.index[tmp.var["highly_variable"]].tolist()
        hvgs_list.append(set(hvgs))

    if not hvgs_list:
        raise ValueError("Empty list of AnnData objects in find_hvg().")

    all_genes = sorted(set.union(*hvgs_list))
    counts = {g: sum(g in s for s in hvgs_list) for g in all_genes}
    hvgs_sorted = sorted(all_genes, key=lambda g: -counts[g])

    return hvgs_sorted[:nfeatures]


# ---------------------------------------------------------------------
# 对应 runNMF (单 AnnData)
# ---------------------------------------------------------------------


def run_nmf(
    adata: ad.AnnData,
    k: int = 10,
    layer: Optional[str] = None,
    use_raw: bool = False,
    hvg: Optional[Sequence[str]] = None,
    new_reduction: str = "NMF",
    center: bool = False,
    scale: bool = False,
    l1: Tuple[float, float] = (0.0, 0.0),
    random_state: int = 123,
) -> ad.AnnData:
    """
    在单个 AnnData 上跑 NMF，把结果存进 adata：
      - adata.obsm["X_<lower(new_reduction)>"] : 细胞嵌入 (cells x k)
      - adata.varm["<new_reduction>_loadings"] : 基因 loading (genes x k)
    """
    if hvg is None or len(hvg) <= 1:
        if "highly_variable" not in adata.var.columns or adata.var[
            "highly_variable"
        ].sum() == 0:
            raise ValueError(
                "No HVG found. Run `scanpy.pp.highly_variable_genes` "
                "or pass `hvg` explicitly."
            )
        genes = adata.var_names[adata.var["highly_variable"].values].tolist()
    else:
        genes = list(hvg)

    X, genes_used, cells_used = get_data_matrix(
        adata,
        layer=layer,
        use_raw=use_raw,
        genes=genes,
        center=center,
        scale=scale,
        non_negative=True,
    )

    nmf_model = NMF(
        n_components=k,
        init="nndsvd",
        random_state=random_state,
        alpha_W=l1[0],
        alpha_H=l1[1],
        l1_ratio=0.0,
        max_iter=200,
    )

    W = nmf_model.fit_transform(X)  # cells x k
    H = nmf_model.components_      # k x genes

    # 存 cell embeddings
    red_key = new_reduction.lower()
    adata.obsm[f"X_{red_key}"] = W

    # 存 gene loadings（对所有基因补 0）
    loadings = np.zeros((adata.n_vars, k), dtype=float)
    idx = adata.var_names.get_indexer(genes_used)
    loadings[idx, :] = H.T
    adata.varm[f"{new_reduction}_loadings"] = loadings

    # 一些元信息
    adata.uns.setdefault(f"{new_reduction}_nmf", {})
    adata.uns[f"{new_reduction}_nmf"].update(
        {"genes": genes_used, "k": int(k), "random_state": int(random_state)}
    )

    return adata


# ---------------------------------------------------------------------
# 对应 multiNMF (多样本、多 k)
# ---------------------------------------------------------------------


def multi_nmf(
    adatas: Sequence[ad.AnnData],
    k: Sequence[int],
    nmf_layer: Optional[str] = None,
    use_raw: bool = False,
    hvg: Optional[Sequence[str]] = None,
    hvg_layer: Optional[str] = None,
    nfeatures: int = 2000,
    l1: Tuple[float, float] = (0.0, 0.0),
    min_exp: float = 0.01,
    max_exp: float = 3.0,
    center: bool = False,
    scale: bool = False,
    min_cells_per_sample: int = 10,
    hvg_blocklist: Optional[Sequence[str]] = None,
    sample_ids: Optional[Sequence[str]] = None,
    random_state: int = 123,
) -> MultiNMFResult:
    """
    近似 R 版 multiNMF：
      - 对每个样本（AnnData）和每个 k 跑一套 NMF
      - 使用全局 HVG（通过 find_hvg 得到）
      - 返回一个 MultiNMFResult，里面 programs 的 key 为 "<sample>.k<k>"
    """
    if not k:
        raise ValueError("Parameter `k` must be a non-empty list/range of integers.")

    ks = sorted(set(int(x) for x in k if int(x) >= 2))
    if len(ks) == 0:
        raise ValueError("All k values < 2; need at least 2 components for NMF.")

    # 过滤掉太小的样本
    filtered_adatas: List[ad.AnnData] = []
    filtered_samples: List[str] = []
    for i, a in enumerate(adatas):
        if a.n_obs <= min_cells_per_sample:
            continue
        filtered_adatas.append(a)
        if sample_ids is None:
            if "sample" in a.obs.columns:
                filtered_samples.append(str(a.obs["sample"].iloc[0]))
            else:
                filtered_samples.append(f"S{i+1}")
        else:
            filtered_samples.append(str(sample_ids[i]))

    if not filtered_adatas:
        raise ValueError("No samples left after filtering by `min_cells_per_sample`.")

    # 全局 HVG
    if hvg is None or len(hvg) <= 1:
        hvg = find_hvg(
            filtered_adatas,
            nfeatures=nfeatures,
            min_exp=min_exp,
            max_exp=max_exp,
            hvg_blocklist=hvg_blocklist,
            layer=hvg_layer,
        )
    else:
        hvg = list(hvg)

    programs: Dict[str, NMFProgram] = {}

    for adata, sample_id in zip(filtered_adatas, filtered_samples):
        genes = [g for g in hvg if g in adata.var_names]
        if len(genes) < 2:
            continue

        X, genes_used, cells_used = get_data_matrix(
            adata,
            layer=nmf_layer,
            use_raw=use_raw,
            genes=genes,
            center=center,
            scale=scale,
            non_negative=True,
        )

        for kk in ks:
            if kk >= X.shape[1]:
                # 组件数不能超过基因数
                continue

            model = NMF(
                n_components=kk,
                init="nndsvd",
                random_state=random_state,
                alpha_W=l1[0],
                alpha_H=l1[1],
                l1_ratio=0.0,
                max_iter=200,
            )

            W = model.fit_transform(X)   # cells x k
            H = model.components_        # k x genes

            nmf_prog = NMFProgram(
                w=H.T.copy(),             # genes x k
                h=W.T.copy(),             # k x cells
                genes=genes_used.copy(),
                cells=cells_used.copy(),
                d=None,
                tol=None,
                n_iter=getattr(model, "n_iter_", None),
            )
            run_id = f"{sample_id}.k{kk}"
            programs[run_id] = nmf_prog

    return MultiNMFResult(
        programs=programs,
        hvg_genes=list(hvg),
        samples=filtered_samples,
        ks=ks,
    )


# ---------------------------------------------------------------------
# 对应 getMetaPrograms (简化版)
# ---------------------------------------------------------------------


def get_meta_programs(
    ms_res: MultiNMFResult,
    nMP: int = 10,
    specificity_weight: float = 5.0,
    weight_explained: float = 0.5,
    max_genes: int = 200,
    metric: str = "cosine",
    min_confidence: float = 0.5,
) -> MetaProgramResult:
    """
    Python 版 meta-program 发现：
      1) 将所有 sample/k 的 NMF component 展开成 'program'
         program_id: "<sample>.k<k>.<component_idx>"
      2) 每个 program 的 gene weight 先按 specificity_weight 做特异性加权
      3) 构建程序间相似度矩阵（默认 cosine）
      4) 对程序做层次聚类为 nMP 个 cluster
      5) 对每个 cluster 取共识 gene weights（平均 + weight_cumul + confidence 过滤）
    """
    # 1. 展开成每个 program 的基因权重
    gene_program_weights: Dict[str, pd.Series] = {}

    for run_id, nmf_prog in ms_res.programs.items():
        W = np.array(nmf_prog.w, dtype=float)   # genes x k
        genes = np.array(nmf_prog.genes)

        # gene specificity: 基因在某一 component 上的相对权重最大值
        # 如果一个基因只在第 1 个因子中高表达，在其他因子中为 0，它的特异性就是 1；
        # 如果它在所有因子中均匀分布，特异性就很低
        row_sums = W.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        rownorm = W / row_sums
        spec = rownorm.max(axis=1)
        spec_w = spec ** specificity_weight  # 这会极大地惩罚那些“到处都表达”的基因
        W_weighted = W * spec_w[:, None]

        # 按列归一化
        col_sums = W_weighted.sum(axis=0, keepdims=True)
        col_sums[col_sums == 0] = 1.0
        W_weighted = W_weighted / col_sums

        for j in range(W_weighted.shape[1]):
            prog_id = f"{run_id}.{j+1}"
            gene_program_weights[prog_id] = pd.Series(W_weighted[:, j], index=genes)

    if not gene_program_weights:
        raise ValueError("No NMF programs found in MultiNMFResult.")

    program_ids = list(gene_program_weights.keys())

    # 2. 程序间相似度矩阵 J
    metric = metric.lower()
    if metric == "cosine":
        gene_table = pd.DataFrame(gene_program_weights)  # genes x programs
        gene_table = gene_table.fillna(0.0)
        J_arr = cosine_similarity(gene_table.T)         # programs x programs
        J = pd.DataFrame(J_arr, index=program_ids, columns=program_ids)
    elif metric == "jaccard":
        n = len(program_ids)
        J_arr = np.zeros((n, n), dtype=float)
        for i, pid_i in enumerate(program_ids):
            gi = set(gene_program_weights[pid_i].index)
            for j, pid_j in enumerate(program_ids):
                gj = set(gene_program_weights[pid_j].index)
                inter = len(gi & gj)
                union = len(gi) + len(gj) - inter
                J_arr[i, j] = inter / union if union > 0 else 0.0
        J = pd.DataFrame(J_arr, index=program_ids, columns=program_ids)
    else:
        raise ValueError("metric must be 'cosine' or 'jaccard'.")

    # 3. 聚类程序
    D = 1.0 - J.values
    np.fill_diagonal(D, 0.0)
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="ward")
    cl_numeric = fcluster(Z, t=nMP, criterion="maxclust")
    clusters = pd.Series(cl_numeric, index=program_ids)

    # 用 precomputed 距离计算 silhouette
    sil = silhouette_samples(D, cl_numeric, metric="precomputed")
    sil_series = pd.Series(sil, index=program_ids)
    sil_means = sil_series.groupby(clusters).mean()

    # 解析 sample 名（假设格式 "<sample>.k<k>.<idx>"）
    pattern = re.compile(r"^(.*?)\.k\d+\.\d+$")
    sample_labels = []
    for pid in program_ids:
        m = pattern.match(pid)
        if m:
            sample_labels.append(m.group(1))
        else:
            sample_labels.append(pid.split(".")[0])
    all_samples = sorted(set(sample_labels))
    sample_series = pd.Series(sample_labels, index=program_ids)

    # 4. 为每个 MP 计算共识基因 + 指标
    metaprograms_genes: Dict[str, List[str]] = {}
    metaprograms_gene_weights: Dict[str, pd.Series] = {}
    metrics_rows = []
    composition_rows = []

    for mp in range(1, nMP + 1):
        mp_label = f"MP{mp}"
        prog_in_mp = clusters.index[clusters == mp]

        if len(prog_in_mp) == 0:
            metaprograms_genes[mp_label] = []
            metaprograms_gene_weights[mp_label] = pd.Series(dtype=float)
            metrics_rows.append(
                dict(
                    sampleCoverage=0.0,
                    silhouette=0.0,
                    meanSimilarity=0.0,
                    numberGenes=0,
                    numberPrograms=0,
                )
            )
            composition_rows.append({s: 0 for s in all_samples})
            continue

        # 平均 gene weight
        df_mp = pd.DataFrame(
            {pid: gene_program_weights[pid] for pid in prog_in_mp}
        ).fillna(0.0)
        avg_weights = df_mp.mean(axis=1)

        # 计算 gene confidence：基因在多少 program 的 top80% 权重里出现
        gene_sets = []
        for pid in prog_in_mp:
            w = gene_program_weights[pid]
            top_for_prog = weight_cumul(w, weight_explained=0.8)
            gene_sets.append(set(top_for_prog.index))

        from collections import Counter

        counter = Counter()
        for sset in gene_sets:
            counter.update(sset)
        confidence = {
            g: counter[g] / float(len(gene_sets)) for g in counter.keys()
        }

        # 按平均权重 weight_cumul 再用 confidence 过滤
        candidate = weight_cumul(avg_weights, weight_explained=weight_explained)
        if min_confidence is not None:
            keep_genes = [
                g for g in candidate.index if confidence.get(g, 0.0) >= min_confidence
            ]
            candidate = candidate.loc[keep_genes]

        candidate = candidate.sort_values(ascending=False).head(max_genes)
        if candidate.sum() > 0:
            candidate = candidate / candidate.sum()

        metaprograms_genes[mp_label] = list(candidate.index)
        metaprograms_gene_weights[mp_label] = candidate

        # 样本覆盖度
        mp_samples = sample_series.loc[prog_in_mp].unique()
        sample_coverage = (
            len(mp_samples) / len(all_samples) if all_samples else 0.0
        )

        # 内部平均相似度
        if len(prog_in_mp) > 1:
            subJ = J.loc[prog_in_mp, prog_in_mp].values
            upper = subJ[np.triu_indices_from(subJ, k=1)]
            mean_sim = float(np.mean(upper))
        else:
            mean_sim = 0.0

        metrics_rows.append(
            dict(
                sampleCoverage=sample_coverage,
                silhouette=float(sil_means.get(mp, 0.0)),
                meanSimilarity=mean_sim,
                numberGenes=len(candidate),
                numberPrograms=len(prog_in_mp),
            )
        )

        # composition: 每个 sample 有多少 program 归入这个 MP
        comp = {}
        mp_sample_series = sample_series.loc[prog_in_mp]
        for s in all_samples:
            comp[s] = int((mp_sample_series == s).sum())
        composition_rows.append(comp)

    metrics_df = pd.DataFrame(metrics_rows, index=[f"MP{i}" for i in range(1, nMP + 1)])
    composition_df = pd.DataFrame(
        composition_rows,
        index=metrics_df.index,
        columns=all_samples,
    )

    # 程序 -> MP 名称映射
    program_clusters = clusters.map(lambda x: f"MP{x}")

    return MetaProgramResult(
        metaprograms_genes=metaprograms_genes,
        metaprograms_gene_weights=metaprograms_gene_weights,
        metaprograms_metrics=metrics_df,
        metaprograms_composition=composition_df,
        programs_similarity=J,
        programs_clusters=program_clusters,
        linkage=Z,
    )


# ---------------------------------------------------------------------
# 对应 dropMetaPrograms (简化版)
# ---------------------------------------------------------------------


def drop_meta_programs(
    mp_res: MetaProgramResult,
    drop_mp: Sequence[str],
) -> MetaProgramResult:
    """
    简化版 dropMetaPrograms:
    - 从 meta-program 列表和指标里删除指定 MP
    - 被删除 MP 所包含的 programs 的 cluster 标记为 "MPNA"
    """
    drop_set = set(drop_mp)
    keep_mps = [mp for mp in mp_res.metaprograms_genes.keys() if mp not in drop_set]

    metaprograms_genes = {
        mp: mp_res.metaprograms_genes[mp] for mp in keep_mps
    }
    metaprograms_gene_weights = {
        mp: mp_res.metaprograms_gene_weights[mp] for mp in keep_mps
    }
    metaprograms_metrics = mp_res.metaprograms_metrics.loc[keep_mps].copy()
    metaprograms_composition = mp_res.metaprograms_composition.loc[
        keep_mps
    ].copy()

    program_clusters = mp_res.programs_clusters.copy()
    # 把被删的 programs 标记为 MPNA
    for pid, mp in program_clusters.items():
        if mp in drop_set:
            program_clusters[pid] = "MPNA"

    # 将剩下的 MP 重新编号为 MP1..MPn（保持原有顺序）
    # keep_mps 在上面已由 mp_res.metaprograms_genes.keys() 过滤得到并保留顺序
    new_mapping = {old_mp: f"MP{i+1}" for i, old_mp in enumerate(keep_mps)}

    # 更新 metaprograms 字典的键为新名字
    metaprograms_genes = {
        new_mapping.get(mp, mp): genes for mp, genes in metaprograms_genes.items()
    }
    metaprograms_gene_weights = {
        new_mapping.get(mp, mp): weights for mp, weights in metaprograms_gene_weights.items()
    }

    # 更新 metrics / composition 的 index 名称
    metaprograms_metrics.index = [new_mapping.get(idx, idx) for idx in metaprograms_metrics.index]
    metaprograms_composition.index = metaprograms_metrics.index.copy()

    # 更新 program_clusters 中剩余 MP 的名字（MPNA 保持不变）
    program_clusters = program_clusters.map(lambda v: new_mapping.get(v, v))

    mp_res = MetaProgramResult(
        metaprograms_genes=metaprograms_genes,
        metaprograms_gene_weights=metaprograms_gene_weights,
        metaprograms_metrics=metaprograms_metrics,
        metaprograms_composition=metaprograms_composition,
        programs_similarity=mp_res.programs_similarity,
        programs_clusters=program_clusters,
        linkage=mp_res.linkage,
    )

    return mp_res

def plot_meta_programs(
    mp_res,
    similarity_cutoff: Tuple[float, float] = (0.0, 1.0),
    scale: str = "none",
    downsample: Optional[int] = None,
    showtree: bool = True,
    cmap: str = "viridis",
    annotation_colors: Optional[Dict[str, str]] = None,
    show_rownames: bool = False,
    show_colnames: bool = False,
    figsize: Tuple[float, float] = (8, 8),
    random_state: Optional[int] = 123,
    **heatmap_kws,
):
    """
    Python 版 plotMetaPrograms：对 get_meta_programs 结果画程序相似度热图。

    参数
    ----
    mp_res : MetaProgramResult
        get_meta_programs() 的返回值，要求至少有：
        - programs_similarity : pd.DataFrame，程序 × 程序 相似度
        - programs_clusters   : pd.Series，index 为程序 ID，值为 'MP1' 之类的标签
        - linkage             : scipy.cluster.hierarchy.linkage 矩阵

    similarity_cutoff : (float, float)
        相似度显示的最小/最大值，对应 R 版的 similarity.cutoff。

    scale : {"none", "row", "column"}
        是否对矩阵按行/列做 z-score，类似 pheatmap 的 scale 选项。

    downsample : int
        程序数太多时，最多显示多少个程序（按 meta-program 分层抽样）。

    showtree : bool
        是否显示树（使用 get_meta_programs 里预先算好的 linkage）。

    cmap : str
        Matplotlib colormap 名称（默认 viridis）。

    annotation_colors : dict or None
        meta-program 到颜色的映射，例如 {"MP1": "red", "MP2": "blue"}。
        如果为 None，则自动给每个 MP 分配一个颜色。

    show_rownames, show_colnames : bool
        是否显示每个程序的名字（行/列标签）。

    figsize : (float, float)
        图像大小。

    random_state : int or None
        下采样的随机种子。

    **heatmap_kws :
        透传给 seaborn.clustermap 的其他参数。
    """
    # 取程序相似度矩阵
    S = mp_res.programs_similarity.copy()
    order = leaves_list(mp_res.linkage)
    S = S.iloc[order, order]

    # drop MPNA programs
    drop_progs = mp_res.programs_clusters.index[
        mp_res.programs_clusters == "MPNA"
    ].tolist()
    S = S.drop(index=drop_progs, columns=drop_progs)

    if not isinstance(S, pd.DataFrame):
        S = pd.DataFrame(S)

    P = S.shape[0]
    prog_ids = S.index.to_list()

    # ---------- 1) 下采样（按 MP 分层抽样） ----------
    if downsample is not None and P > downsample:
        rng = np.random.default_rng(random_state)
        labels = mp_res.programs_clusters.reindex(prog_ids)
        unique_mps = labels.dropna().unique().tolist()
        per_mp = max(downsample // max(len(unique_mps), 1), 1)

        keep = []
        for mp in unique_mps:
            idx = labels[labels == mp].index.to_numpy()
            if len(idx) <= per_mp:
                keep.extend(idx.tolist())
            else:
                keep.extend(rng.choice(idx, size=per_mp, replace=False).tolist())

        # 保持原有顺序
        prog_ids = [pid for pid in prog_ids if pid in keep]
        S = S.loc[prog_ids, prog_ids]

    # ---------- 2) 相似度截断 & scaling ----------
    vmin, vmax = similarity_cutoff
    S = S.clip(lower=vmin, upper=vmax)

    mat = S.values.astype(float)
    if scale == "row":
        mean = mat.mean(axis=1, keepdims=True)
        std = mat.std(axis=1, keepdims=True) + 1e-9
        mat = (mat - mean) / std
    elif scale == "column":
        mean = mat.mean(axis=0, keepdims=True)
        std = mat.std(axis=0, keepdims=True) + 1e-9
        mat = (mat - mean) / std

    S_scaled = pd.DataFrame(mat, index=prog_ids, columns=prog_ids)

    # ---------- 3) MP 颜色注释 ----------
    mp_labels = mp_res.programs_clusters.reindex(prog_ids)
    # 按 MP 后的数字排序（若非 MP<number> 则放在后面按字母序）
    mp_vals = mp_labels.dropna().unique().tolist()
    pattern_mp = re.compile(r"^MP(\d+)$")

    def _mp_sort_key(x):
        m = pattern_mp.match(x)
        if m:
            return (0, int(m.group(1)))
        return (1, x)

    MPs = sorted(mp_vals, key=_mp_sort_key)

    if annotation_colors is None:
        pal = sns.color_palette("tab20", len(MPs))
        annotation_colors = {mp: pal[i] for i, mp in enumerate(MPs)}

    row_colors = mp_labels.map(annotation_colors)

    # ---------- 4) 树状图 / 聚类 ----------
    row_linkage = col_linkage = None
    if showtree and len(prog_ids) > 1:
        if downsample is None:
            # 没有下采样，直接复用 get_meta_programs 算好的 linkage
            row_linkage = col_linkage = None
            row_cluster = col_cluster = False
        else:
            # 下采样后重新做一次聚类
            D = 1.0 - S_scaled.values
            np.fill_diagonal(D, 0.0)
            condensed = squareform(D, checks=False)
            row_linkage = col_linkage = linkage(condensed, method="average")
            row_cluster = col_cluster = True

    # ---------- 5) 行/列名字 ----------
    ytick = prog_ids if show_rownames else False
    xtick = prog_ids if show_colnames else False

    # ---------- 6) clustermap ----------
    g = sns.clustermap(
        S_scaled,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        row_colors=row_colors,
        col_colors=row_colors,
        cmap=cmap,
        xticklabels=xtick,
        yticklabels=ytick,
        figsize=figsize,
        cbar_pos=None,
        **heatmap_kws,
    )


    heatmap = g.ax_heatmap.collections[0]

    # 创建图例的 "handles" (色块)：基于已有的 MP 标签和 annotation_colors 构造
    # MPs 在上面已由 mp_labels 推导得到
    legend_handles = [
        mpatches.Patch(color=annotation_colors.get(mp, "#808080"), label=mp)
        for mp in MPs
    ]

    # 将图例添加到 Figure 上
    # bbox_to_anchor 控制位置，loc 控制对齐方式
    g.fig.legend(
        handles=legend_handles,
        title="Metaprogram",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.98),
        frameon=False
    )

    # 在右侧添加一个新的 axes 放 colorbar
    try:
        # 尝试基于刚添加的 legend 的位置在其正下方添加 colorbar axes
        if not g.fig.legends:
            raise RuntimeError("No legend found")
        legend = g.fig.legends[-1]
        renderer = g.fig.canvas.get_renderer()
        bbox = legend.get_window_extent(renderer).transformed(g.fig.transFigure.inverted())
        pad = 0.03
        height = bbox.height / len(MPs) * 4
        left = (bbox.x1 - bbox.x0) / 5 + bbox.x0
        width = (bbox.x1 - bbox.x0) / 5
        bottom = bbox.y0 - height - pad
        if bottom < 0:
            bottom = 0.01
        cbar_ax = g.fig.add_axes([left, bottom, width, height])
    except Exception:
        # 回退到固定位置（如果无法计算 legend bbox）
        cbar_ax = g.fig.add_axes([1.02, 0.8, 0.02, 0.15])
    cbar = plt.colorbar(heatmap, cax=cbar_ax)
    # cbar.outline.set_visible(False)
    cbar.set_label("Similarity", labelpad=5)

    return g