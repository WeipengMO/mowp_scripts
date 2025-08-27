import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from typing import Union
from typing import Literal, Optional


def grouped_obs_mean(adata: ad.AnnData, group_key: str, layer: str = None) -> pd.DataFrame:
    """
    Calculate the mean expression of each gene in each group of cells.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    group_key : str
        The key of the observation grouping to consider.
    layer : str
        The layer to use. If None, use `adata.X`.
    
    Returns
    -------
    pd.DataFrame
        A DataFrame with the mean expression of each gene in each group of cells.
    """
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
        
    return out


def var_means(adata, keys: Union[list, dict], inplace=False, uns: str = 'mean_expression'):
    if isinstance(keys, list):
        keys = {'mean_expression': keys}
    
    results = []
    for k, genes in keys.items():
        common_genes = list(adata.var_names.intersection(genes))
        if len(common_genes) == 0:
            raise ValueError(f'No common genes found in {k}.')

        df = sc.get.obs_df(adata, common_genes)
        df[k] = df.mean(axis=1)
        results.append(df[[k]])
    
    df = pd.concat(results, axis=1)

    if inplace:
        adata.uns[uns] = df
    else:
        return df



def jaccard_index(set1: set, set2: set) -> float:
    set1 = set(set1)
    set2 = set(set2)

    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    return intersection / union


def _otsu_threshold(x: np.ndarray, bins: int = 256) -> float:
    """
    简单 Otsu 阈值（在非负连续值上通过直方图近似）。
    仅用于单峰/双峰明显时的粗分割。
    """
    # x = x[np.isfinite(x)]
    x = x[x >= 0]
    if x.size == 0:
        return 0.0
    hist, bin_edges = np.histogram(x, bins=bins, range=(x.min(), x.max() if x.max() > x.min() else x.min() + 1e-6))
    hist = hist.astype(float)
    p = hist / hist.sum()

    omega = np.cumsum(p)
    mu = np.cumsum(p * (bin_edges[:-1] + bin_edges[1:]) / 2.0)
    mu_t = mu[-1]

    # 类间方差
    sigma_b2 = (mu_t * omega - mu)**2 / (omega * (1 - omega) + 1e-12)
    idx = np.nanargmax(sigma_b2)
    # 取该 bin 的中心作为阈值
    thr = (bin_edges[idx] + bin_edges[idx + 1]) / 2.0
    return float(thr)


try:
    from sklearn.mixture import GaussianMixture
    _HAS_SKLEARN = True
except Exception:
    _HAS_SKLEARN = False
    

def split_cells_by_gene(
    adata: sc.AnnData,
    gene: str,
    method: Literal["auto", "median", "quantile", "otsu", "gmm"] = "auto",
    quantile: float = 0.5,
    min_cells_per_group: int = 25,
    obs_key: Optional[str] = None,
    layer: Optional[str] = None,
) -> str:
    """
    根据某个基因表达把细胞分成 high/low 两组，并把标签写入 adata.obs。

    Parameters
    ----------
    adata : AnnData
    gene : 基因名（需在 adata.var_names 中）
    method : 分组方法
        - "auto": 先尝试 Otsu，若失败则退化到 median；若安装了 sklearn 且双峰明显则用 GMM
        - "median": 中位数二分
        - "quantile": 按给定分位数二分（如 0.3/0.7 也可，两侧需自行改造）
        - "otsu": Otsu 阈值
        - "gmm": 高斯混合（2 组件）
    quantile : 用于 "quantile" 的分位数（默认 0.5）
    min_cells_per_group : 两组最少细胞数，过小会报错
    obs_key : 写入 adata.obs 的列名，默认 f"{gene}_group"
    layer : 若指定，则使用该 layer 的表达矩阵；否则用 adata.X

    Returns
    -------
    obs_key : 实际写入的列名
    """
    if (gene not in adata.var_names) and (gene not in adata.obs.columns):
        raise ValueError(f"Gene '{gene}' not in adata.")

    obs_key = obs_key or f"{gene}_group"

    x = sc.get.obs_df(adata, keys=gene, layer=layer).values

    # 选择阈值
    thr = None
    chosen = method

    def _safe_gmm_threshold(vals: np.ndarray) -> Optional[float]:
        vals = vals[np.isfinite(vals)]
        if vals.size < 100 or not _HAS_SKLEARN:
            return None
        # 避免负值对数域等问题，这里不做 log；直接在当前尺度拟合
        gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=0)
        try:
            gmm.fit(vals.reshape(-1, 1))
            means = np.sort(gmm.means_.ravel())
            # 用两个均值的中点作为切分
            return float((means[0] + means[1]) / 2.0)
        except Exception:
            return None

    if method == "median":
        thr = float(np.nanmedian(x))
    elif method == "quantile":
        thr = float(np.nanquantile(x, quantile))
    elif method == "otsu":
        thr = _otsu_threshold(x)
    elif method == "gmm":
        thr = _safe_gmm_threshold(x)
        if thr is None:
            raise RuntimeError("GMM 阈值计算失败，可能样本过少或未安装 scikit-learn。")
    elif method == "auto":
        # 先试 GMM（若可用），不行再 Otsu，最后用 median
        thr = _safe_gmm_threshold(x)
        if thr is not None:
            chosen = "gmm"
        else:
            try:
                thr = _otsu_threshold(x)
                chosen = "otsu"
            except Exception:
                thr = float(np.nanmedian(x))
                chosen = "median"
    else:
        raise ValueError("method 必须是 'auto'/'median'/'quantile'/'otsu'/'gmm' 之一。")

    # 打标签
    labels = np.where(x >= thr, "high", "low")
    adata.obs[obs_key] = pd.Categorical(labels, categories=["low", "high"], ordered=True)

    # 质量检查
    counts = adata.obs[obs_key].value_counts()
    if counts.min() < min_cells_per_group:
        raise RuntimeError(
            f"分组过于不平衡：{dict(counts)}；请尝试更换 method（如 'median' 或调整 'quantile'），"
            f"或降低 min_cells_per_group。"
        )

    print(f"[INFO] 使用 '{chosen}' 方法，阈值={thr:.4g}，分组：{dict(counts)}。写入列：{obs_key}")
    return obs_key, thr