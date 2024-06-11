import numpy as np
import pandas as pd
import anndata as ad


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



def jaccard_index(set1: set, set2: set) -> float:
    set1 = set(set1)
    set2 = set(set2)

    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    return intersection / union