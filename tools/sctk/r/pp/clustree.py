import scanpy as sc
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from ...utils import rtools


def clustree(
        adata: sc.AnnData,
        prefix: str = 'leiden',
        figsize: tuple = (600, 600),
        dpi=75,
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

    with rtools.r_inline_plot(*figsize, dpi=dpi):
        print(clustree_r.clustree(leien_cluster_r, prefix=f'{prefix}_'))