from tqdm.auto import tqdm
import scanpy as sc
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from ..utils import rtools, configure_logger
from loguru import logger
import seaborn as sns
from matplotlib import pyplot as plt 
import multiprocessing


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



def clustree(
        adata: sc.AnnData,
        prefix: str = 'leiden',
        figsize: tuple = (1000, 2000),
        ):
    """Plot clustree of Leiden clusters.

    Parameters
    ----------
    adata
        Annotated data matrix.
    prefix
        Prefix of clusters. eg. 'leiden', 'louvain'.
    figsize
        Figure size.
    """
    
    clustree_r = importr("clustree")
    pattern = f"{prefix}_\d+\.\d+"
    leien_cluster_r = rtools.py2r(adata.obs.filter(regex=pattern))

    with rtools.r_inline_plot(*figsize, dpi=100):
        print(clustree_r.clustree(leien_cluster_r, prefix=f'{prefix}_'))


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