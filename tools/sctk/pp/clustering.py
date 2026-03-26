import scanpy as sc
import numpy as np
import pandas as pd
from ..settings import configure_logger
from loguru import logger
import seaborn as sns
from matplotlib import pyplot as plt 
import multiprocessing
from joblib import Parallel, delayed
import os


def _run_leiden(adata, resolution):
    ad = sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{resolution}", copy=True)
    return ad


def leiden_iter(
        adata: sc.AnnData, 
        res_start: float = .1, 
        res_end: float = 1, 
        res_step: float = .1):
    """Iterate over Leiden resolutions.

    Parameters
    ----------
    adata
        Annotated data matrix.
    res_start
        Starting resolution.
    res_end
        Ending resolution.
    res_step
        Resolution step size.
    """

    resolutions = np.arange(res_start, res_end+res_step, res_step)
    resolutions = [round(i, 2) for i in resolutions]

    with multiprocessing.Pool(processes=len(resolutions)) as pool:
        results = pool.starmap(_run_leiden, [(adata, res) for res in resolutions])
    
    for ad, res in zip(results, resolutions):
        adata.obs = adata.obs.merge(ad.obs[[f'leiden_{res}']], left_index=True, right_index=True)


def _one_res(g, r, weights):
    import leidenalg
    
    part = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=float(r),
        weights=weights,
        n_iterations=-1,
        seed=0,
    )
    groups = np.array(part.membership).astype(str)
    return r, groups


def leiden_parallel(
        adata: sc.AnnData, 
        res_start: float = .1, 
        res_end: float = 1, 
        res_step: float = .1):
    """Leiden clustering in parallel over multiple resolutions.

    Parameters
    ----------
    adata
        Annotated data matrix.
    res_start
        Starting resolution.
    res_end
        Ending resolution.
    res_step
        Resolution step size.
    """
    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)


    resolutions = np.round(np.arange(res_start, res_end + res_step, res_step), 2)
    resolutions = resolutions[resolutions <= res_end]
    print(f"Leiden for resolutions: {', '.join(map(str, resolutions))}")

    if "connectivities" not in adata.obsp:
        raise ValueError("Please run `sc.pp.neighbors` first to compute the neighbor graph.")
    
    adjacency = sc._utils._choose_graph(adata, obsp='connectivities', neighbors_key=None)
    g = sc._utils.get_igraph_from_adjacency(adjacency, directed=True)
    weights = np.array(g.es["weight"]).astype(np.float64)
    n_jobs = min(len(resolutions), os.cpu_count() or 1)

    pairs = Parallel(n_jobs=n_jobs, prefer="processes")(delayed(_one_res)(g, r, weights) for r in resolutions)

    for r, labels in pairs:
        adata.obs[f"leiden_{r}"] = pd.Categorical(labels)


def cluster_silhouette_score(
    adata,
    cluster_key="leiden",
    obsm_key="X_scVI",
    n_dims=None,
    metric="euclidean",
    sample_size=10000,
    random_state=42,
    show=True,
):
    from sklearn.metrics import silhouette_score

    X = adata.obsm[obsm_key]

    if n_dims is not None and X.shape[1] > n_dims:
        X = X[:, :n_dims]

    prefix = f"{cluster_key}_"
    obs_keys = [c for c in adata.obs.columns if c.startswith(prefix)]

    dt_score = {}

    sample_size = sample_size if (sample_size and adata.n_obs > sample_size) else None

    for key in obs_keys:
        labels = adata.obs[key].astype(str)

        try:
            res = float(key.replace(prefix, ""))
        except ValueError:
            continue

        n_clusters = labels.nunique()
        if n_clusters <= 1 or n_clusters >= len(labels):
            dt_score[res] = np.nan
            continue

        score = silhouette_score(
            X,
            labels,
            metric=metric,
            sample_size=sample_size,
            random_state=random_state,
        )
        dt_score[res] = score
    
    if show:
        plt.figure(figsize=(4, 4))
        sns.pointplot(x=dt_score.keys(), y=dt_score.values())
        plt.xlabel('resolution')
        plt.ylabel('silhouette score')
        plt.show()

    return dt_score