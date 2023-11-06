import scanpy as sc
import anndata
from loguru import logger
from typing import Optional, Union, Sequence, List
import numpy as np


def scanpy_pp(adata, n_top_genes=2000, n_pcs=50, n_neighbors=15, resolution=0.5, seed=1, inplace=True):
    """
    Preprocesses the AnnData object by performing the following steps:
    1. Select highly variable genes
    2. Normalize total count
    3. Log transform
    4. PCA
    5. Nearest neighbor graph
    6. UMAP
    7. Leiden clustering

    Parameters
    ----------
    adata : AnnData
        The AnnData object to be preprocessed.
    n_top_genes : int, optional
        Number of highly variable genes to select. The default is 2000.
    n_pcs : int, optional
        Number of principal components to use. The default is 50.
    n_neighbors : int, optional
        Number of neighbors to use for nearest neighbor graph. The default is 15.
    resolution : float, optional
        Resolution parameter for Leiden clustering. The default is 0.5.
    seed : int, optional
        Random seed. The default is 1.
    inplace : bool, optional
        Whether to perform the preprocessing inplace. The default is True.
    
    Returns
    -------
    adata : AnnData
        The preprocessed AnnData object.

    """

    if not inplace:
        adata = adata.copy()

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack', random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, random_state=seed)
    
    return adata if not inplace else None


def layer_pp(adata, layer=None, total=1e4, logbase=10, scale=False):
    """
    Preprocesses the AnnData object by performing the following steps:
    1. Normalize total count
    2. Log transform
    3. Scale, if specified

    Parameters
    ----------
    adata : AnnData
        The AnnData object to be preprocessed.
    layer : str, optional
        The layer to be preprocessed. The default is None.
    total : int, optional
        Total count to normalize to. The default is 1e4.
    logbase : int, optional
        Base of the log transform. The default is 10.
    scale : bool, optional
        Whether to scale the data. The default is False.
    """

    adata.layers['raw'] = adata.X.copy() if layer == None else adata.layers[layer].copy()
    adata.layers['normalize_log'] = adata.layers['raw'].copy()
    sc.pp.normalize_total(adata, total, layer='normalize_log')
    sc.pp.log1p(adata, layer='normalize_log', base=logbase)
    if scale:
        adata.layers['normalize_log_scale'] = adata.layers['normalize_log'].copy()
        sc.pp.scale(adata, layer='normalize_log_scale')

    adata.X = adata.layers['normalize_log'].copy()


def simplify_adata(
    adata,
    layers: Optional[Union[str, List[str]]] = None,
    obs_keys: Optional[Sequence[str]] = None,
    var_keys: Optional[Sequence[str]] = None,
) -> sc.AnnData:
    """
    Get a subset of adata. Only contains one layer expression matrix, and several obs information.

    Parameters
    ----------
    adata
        anndata.AnnData
    layers
        layer name or list of layer names
    obs_keys
        obs information keys
    var_keys
        var information keys

    Returns
    -------
    anndata.AnnData
        if data type of `layers` is list, all element in 'X' of returned adata will be set as 0
    """
    import scipy.sparse as ss

    if not obs_keys:
        obs_keys = []
    else:
        assert sum([x in adata.obs.columns for x in obs_keys]) == len(
            obs_keys
        ), "Requested feature not located in adata.obs"

    if not var_keys:
        var_keys = []
    else:
        assert sum([x in adata.var.columns for x in var_keys]) == len(
            var_keys
        ), "Requested feature not located in adata.var"

    if not layers:
        layers = "X"

    if isinstance(layers, list):
        dt_layerMtx = {}
        for layer in layers:
            ar_mtx = adata.X if layer == "X" else adata.layers[layer]
            dt_layerMtx[layer] = ar_mtx
        subAd = anndata.AnnData(
            ss.csr_matrix(np.zeros(adata.shape)),
            adata.obs[obs_keys],
            adata.var[var_keys],
            layers=dt_layerMtx,
        )
        logger.info('Adata.X is empty')

    elif isinstance(layers, str):
        layer = layers
        mtxAr = adata.X if layer == "X" else adata.layers[layer]
        subAd = anndata.AnnData(mtxAr, adata.obs[obs_keys], adata.var[var_keys])

    else:
        assert False, f"unsupported layers data type: {type(layers)}"

    return subAd.copy()