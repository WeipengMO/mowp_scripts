import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from typing import Union


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