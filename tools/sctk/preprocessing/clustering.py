from tqdm import tqdm
import scanpy as sc
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from ..utils import rtools


def leiden_iter(
        adata: sc.AnnData, 
        res_start: float = .1, 
        res_end: float = 1.5, 
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

    resolution = np.arange(res_start, res_end+res_step, res_step)
    resolution = [round(i, 2) for i in resolution]

    for res in tqdm(resolution, desc="Leiden clustering"):
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")



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
    leien_cluster_r = rtools.py2r(adata.obs.filter(like=prefix))

    with rtools.r_inline_plot(*figsize, dpi=100):
        print(clustree_r.clustree(leien_cluster_r, prefix=f'{prefix}_'))