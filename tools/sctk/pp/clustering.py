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


def get_shilouette_score(
        adata: sc.AnnData, 
        obsm: str = 'X_pca', 
        cluster_key: str = 'leiden',
        metric:str = 'euclidean',
        n_pcs: int = 50,
        debug: bool = False,
        show: bool = True,
    ):
    '''The function performs clustering using the Leiden algorithm on an AnnData object and calculates the silhouette score for each clustering result.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix containing the input data.

    obsm : str, optional
        The name of the observation matrix in `adata` to be used for clustering. Default is 'X_pca'.

    cluster_key : str, optional
        The key in `adata.obs` that specifies the clustering labels. Default is 'leiden'.

    metric : str, optional
        The distance metric to be used for calculating pairwise distances between observations. Default is 'euclidean'.
        Other possible values include 'manhattan' for Manhattan distance, 'cosine' for cosine similarity, and many more.

    n_pcs : int, optional
        The number of principal components to use for clustering. Default is 50.

    debug : bool, optional
        If True, enables debug logging. Default is False.

    show : bool, optional
        If True, displays a point plot of silhouette scores. Default is True.

    Returns
    -------
    dict
        A dictionary where the keys are the resolution values from the input list `ls_res` and the values are the corresponding silhouette scores calculated using the Leiden clustering algorithm.
    '''

    from sklearn.metrics import pairwise_distances, silhouette_score

    if debug:
        configure_logger(log_level="debug")
    else:
        configure_logger(log_level="info")

    _ad = adata
    cluster_key = cluster_key + '_'

    logger.debug(f"{obsm=}")
    if isinstance(obsm, str):
        obsm = _ad.obsm[obsm]
    logger.debug(f"{n_pcs=}")
    if obsm.shape[1] > n_pcs:
        obsm = obsm[:, :n_pcs]

    ar_dist = pairwise_distances(obsm, metric=metric)
    dt_score = {}
    obs_key = list(_ad.obs.filter(like=cluster_key).columns)

    dt_score = {}
    obs_key = list(adata.obs.filter(like=cluster_key).columns)
    logger.debug(f"{obs_key=}")

    for _obs_key in obs_key:
        ls_label =  adata.obs[_obs_key]
        res = float(_obs_key.replace(cluster_key, ''))
        if len(set(ls_label)) == 1:
            dt_score[res] = 0
        else:
            dt_score[res] = silhouette_score(ar_dist, adata.obs[_obs_key], metric='precomputed')

    if show:
        sns.pointplot(x=dt_score.keys(), y=dt_score.values())
        plt.xlabel('resolution')
        plt.ylabel('shilouette score')
        plt.show()

    return dt_score