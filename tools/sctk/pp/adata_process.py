import anndata
import colorsys
from typing import Any, Iterator, List, Literal, Optional, Sequence, Union
from collections.abc import Mapping, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like, to_hex

import numpy as np
import scanpy as sc
import seaborn as sns
from loguru import logger
from tqdm.auto import tqdm


# ============================================================
# Preprocessing function for adata
# ============================================================

def scanpy_pp(
    adata: anndata.AnnData,
    *,
    n_top_genes: int = 3000,
    seed: int = 1,
    target_sum: float = 1e4,
    hvg_flavor: str = "seurat_v3",
    batch_key: str = None,
    counts_layer: str = "counts",
    log1p_layer: str = "log1p_norm",
    reset_from_counts: bool = True,
    inplace: bool = True,
) -> anndata.AnnData:
    """
    Standard Scanpy preprocessing workflow.

    The function performs:
    1. Save raw counts into `.layers[counts_layer]`
    2. Select highly variable genes
    3. Normalize total counts
    4. Log1p transform
    5. PCA
    6. Neighbor graph construction
    7. UMAP

    Parameters
    ----------
    adata
        AnnData object.
    n_top_genes
        Number of highly variable genes.
    n_pcs
        Number of principal components used for neighbor graph construction.
    n_neighbors
        Number of neighbors used for neighbor graph construction.
    seed
        Random seed.
    target_sum
        Target sum for total-count normalization.
    hvg_flavor
        HVG selection method.
        For "seurat_v3" and "seurat_v3_paper", raw counts are used.
    batch_key
        Optional batch key for batch-aware HVG selection.
    counts_layer
        Layer name used to store raw counts.
    log1p_layer
        Layer name used to store normalized and log-transformed expression.
        If None, no extra layer is saved.
    reset_from_counts
        If True and counts_layer already exists, reset `.X` from this layer
        before normalization.
    inplace
        Whether to modify input AnnData in place.

    Returns
    -------
    AnnData or None
        If inplace=False, returns the processed AnnData.
        If inplace=True, modifies `adata` and returns None.
    """

    if not inplace:
        adata = adata.copy()

    # -----------------------------
    # 1. Save or restore raw counts
    # -----------------------------
    if counts_layer in adata.layers:
        if reset_from_counts:
            adata.X = adata.layers[counts_layer].copy()
            print(f"Reset X from adata.layers['{counts_layer}']")
    else:
        adata.layers[counts_layer] = adata.X.copy()
        print(f"Saved X to adata.layers['{counts_layer}']")

    # -----------------------------
    # 2. Highly variable genes
    # -----------------------------
    if hvg_flavor in {"seurat_v3", "seurat_v3_paper"}:
        sc.pp.highly_variable_genes(
            adata,
            flavor=hvg_flavor,
            n_top_genes=n_top_genes,
            layer=counts_layer,
            batch_key=batch_key,
            subset=False,
        )

    # -----------------------------
    # 3. Normalize and log-transform
    # -----------------------------
    sc.pp.normalize_total(
        adata,
        target_sum=target_sum,
    )

    sc.pp.log1p(adata)

    if log1p_layer is not None:
        adata.layers[log1p_layer] = adata.X.copy()
        print(f"Saved log-normalized X to adata.layers['{log1p_layer}']")

    # For non-seurat_v3 flavors, run HVG after log-normalization
    if hvg_flavor not in {"seurat_v3", "seurat_v3_paper"}:
        sc.pp.highly_variable_genes(
            adata,
            flavor=hvg_flavor,
            n_top_genes=n_top_genes,
            batch_key=batch_key,
            subset=False,
        )

    # -----------------------------
    # 4. PCA
    # -----------------------------
    sc.tl.pca(
        adata,
        use_highly_variable=True,
        svd_solver="arpack",
        random_state=seed,
    )

    # -----------------------------
    # 5. Neighbors and UMAP
    # -----------------------------
    sc.pp.neighbors(
        adata,
        random_state=seed,
    )

    sc.tl.umap(
        adata,
        random_state=seed,
    )

    if not inplace:
        return adata

    return None


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
    if layer is None:
        layer = 'counts'

    if layer not in adata.layers:
        adata.layers['counts'] = adata.X.copy()
        logger.info(f'Save X to adat.layers["{layer}"]')
    else:
        logger.info(f'adata.layers["{layer}"] already exists')
    
    logger.info(f'Preprocessing layer: {layer}, normalize total: {total}, logbase: {logbase}')
    adata.layers['log1p_norm'] = adata.layers[layer].copy()
    sc.pp.normalize_total(adata, total, layer='log1p_norm')
    sc.pp.log1p(adata, layer='log1p_norm', base=logbase)
    if scale:
        adata.layers['scale'] = adata.layers['log1p_norm'].copy()
        sc.pp.scale(adata, layer='scale')

    adata.X = adata.layers['log1p_norm'].copy()


# ============================================================
# Color utilities for adata.obs categorical colors
# ============================================================

def _ensure_categorical_obs(
    adata: anndata.AnnData,
    label: str,
) -> None:
    """Ensure adata.obs[label] exists and is categorical."""
    if label not in adata.obs.columns:
        raise KeyError(f"`{label}` was not found in `adata.obs`.")

    if not hasattr(adata.obs[label], "cat"):
        adata.obs[label] = adata.obs[label].astype("category")


def _make_default_colors(
    n_colors: int,
    *,
    as_hex: bool = True,
    seed: int = 1,
    lightness: float = 0.65,
    saturation: float = 0.75,
) -> list:
    """
    Generate default categorical colors.

    Rules
    -----
    1. If n_colors <= 20, use matplotlib tab20.
    2. If n_colors > 20, generate HLS colors with shuffled hue order.
    3. The random order is controlled by seed, so the same n_colors and seed
       always produce the same color order.

    Parameters
    ----------
    n_colors
        Number of colors to generate.
    as_hex
        Whether to return HEX colors.
    seed
        Random seed for shuffling HLS colors.
    lightness
        Lightness value used by HLS colors.
    saturation
        Saturation value used by HLS colors.

    Returns
    -------
    list
        List of colors.
    """

    if n_colors <= 0:
        return []

    # Use tab20 for small categorical palettes
    if n_colors <= 20:
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i)[:3] for i in range(n_colors)]

    # Use shuffled HLS colors for larger categorical palettes
    else:
        rng = np.random.default_rng(seed)

        hues = np.linspace(0, 1, n_colors, endpoint=False)
        rng.shuffle(hues)

        colors = [
            colorsys.hls_to_rgb(
                hue,
                lightness,
                saturation,
            )
            for hue in hues
        ]

    if as_hex:
        return [to_hex(color) for color in colors]

    return colors


def _normalize_color(
    color: Any,
    *,
    as_hex: bool = True,
) -> str:
    """Validate and optionally convert color to HEX."""
    if not is_color_like(color):
        raise ValueError(f"Invalid matplotlib color: {color!r}")

    if as_hex:
        return to_hex(color)

    return color


def get_adata_color(
    adata: anndata.AnnData,
    label: str,
    *,
    as_hex: bool = True,
    reset_if_invalid: bool = True,
) -> dict:
    """
    Get category-color mapping from `adata.uns[f"{label}_colors"]`.

    Parameters
    ----------
    adata
        AnnData object.
    label
        Column name in `adata.obs`.
    as_hex
        Whether to return colors as HEX strings.
    reset_if_invalid
        If True, regenerate default colors when existing colors are missing
        or the color length does not match categories.

    Returns
    -------
    dict
        Mapping from category to color.
    """

    _ensure_categorical_obs(adata, label)

    categories = list(adata.obs[label].cat.categories)
    color_key = f"{label}_colors"

    has_valid_colors = (
        color_key in adata.uns
        and len(adata.uns[color_key]) == len(categories)
    )

    if not has_valid_colors:
        if not reset_if_invalid:
            raise ValueError(
                f"`adata.uns['{color_key}']` is missing or does not match "
                f"the number of categories in `adata.obs['{label}']`."
            )

        adata.uns[color_key] = _make_default_colors(
            len(categories),
            as_hex=True,
        )

    colors = list(adata.uns[color_key])

    if as_hex:
        colors = [to_hex(color) for color in colors]

    return dict(zip(categories, colors))


def set_adata_color(
    adata: anndata.AnnData,
    label: str,
    color_dict: Mapping[Any, Any] = None,
    *,
    as_hex: bool = True,
    strict: bool = True,
) -> anndata.AnnData:
    """
    Set colors for a categorical column in `adata.obs`.

    Parameters
    ----------
    adata
        AnnData object.
    label
        Column name in `adata.obs`.
    color_dict
        Mapping from category to color. It can contain all or part of the
        categories. Existing or default colors will be used for categories
        not included in `color_dict`.
    as_hex
        Whether to store colors as HEX strings.
    strict
        If True, raise an error when `color_dict` contains categories not
        present in `adata.obs[label]`.

    Returns
    -------
    AnnData
        The modified AnnData object.
    """

    _ensure_categorical_obs(adata, label)

    categories = list(adata.obs[label].cat.categories)
    color_key = f"{label}_colors"

    color_map = get_adata_color(
        adata,
        label,
        as_hex=as_hex,
        reset_if_invalid=True,
    )

    if color_dict is not None:
        unknown_categories = set(color_dict) - set(categories)

        if unknown_categories and strict:
            raise KeyError(
                f"The following categories were not found in "
                f"`adata.obs['{label}']`: {sorted(unknown_categories)}"
            )

        for category, color in color_dict.items():
            if category not in categories:
                continue

            color_map[category] = _normalize_color(
                color,
                as_hex=as_hex,
            )

    adata.uns[color_key] = [color_map[category] for category in categories]

    return adata


# ============================================================
# AnnData splitting utilities
# ============================================================

def split_adata(
    adata: anndata.AnnData,
    batch_key: str,
    *,
    copy: bool = True,
    axis: Literal[0, "cell", "obs", 1, "feature", "var"] = 0,
    disable_bar: bool = False,
    dropna: bool = True,
) -> Iterator[anndata.AnnData]:
    """
    Split an AnnData object by groups in `.obs` or `.var`.

    Parameters
    ----------
    adata
        AnnData object.
    batch_key
        Column name used for splitting.
    copy
        Whether to return copied AnnData objects.
        If False, return AnnData views.
    axis
        Split by cells or features.
        Use 0 / "cell" / "obs" for `.obs`;
        use 1 / "feature" / "var" for `.var`.
    disable_bar
        Whether to disable tqdm progress bar.
    dropna
        Whether to drop NA groups.
        This follows pandas groupby behavior. Default is True.

    Yields
    ------
    AnnData
        Subset AnnData for each group.
    """
    if axis in (0, "cell", "obs"):
        metadata = adata.obs
        slice_func = lambda positions: adata[positions, :]
    elif axis in (1, "feature", "var"):
        metadata = adata.var
        slice_func = lambda positions: adata[:, positions]
    else:
        raise ValueError(
            "`axis` must be one of 0, 'cell', 'obs', 1, 'feature', or 'var'."
        )

    if batch_key not in metadata.columns:
        raise KeyError(f"`{batch_key}` was not found in adata.{metadata_name(axis)}.")

    group_indices = metadata[batch_key].groupby(
        metadata[batch_key],
        sort=False,
        observed=True,
        dropna=dropna,
    ).indices

    for _, positions in tqdm(
        group_indices.items(),
        total=len(group_indices),
        disable=disable_bar,
    ):
        subset = slice_func(positions)
        if copy:
            subset = subset.copy()
        yield subset


def metadata_name(axis: Literal[0, "cell", "obs", 1, "feature", "var"]) -> str:
    """Return metadata slot name for error messages."""
    if axis in (0, "cell", "obs"):
        return "obs"
    if axis in (1, "feature", "var"):
        return "var"
    return "obs/var"


# ============================================================
# Legacy code, not recommended to use, will be removed in future version.
# ============================================================

def split_adata_legacy(
    adata: anndata.AnnData,
    batchKey: str,
    copy=True,
    axis: Literal[0, "cell", 1, "feature"] = 0,
    needName=False,
    disableBar=False
) -> Iterator[anndata.AnnData]:
    if axis in [0, "cell"]:
        assert batchKey in adata.obs.columns, f"{batchKey} not detected in adata"
        indexName = "index" if (not adata.obs.index.name) else adata.obs.index.name
        adata.obs["__group"] = adata.obs[batchKey]
        batchObsLs = (
            adata.obs.filter(["__group"])
            .reset_index()
            .groupby("__group")[indexName]
            .agg(list)
        )
        for batchObs in tqdm(batchObsLs, disable=disableBar):
            if needName:
                if copy:
                    yield adata[batchObs].obs.iloc[0].loc["__group"], adata[
                        batchObs
                    ].copy()
                else:
                    yield adata[batchObs].obs.iloc[0].loc["__group"], adata[batchObs]
            else:
                if copy:
                    yield adata[batchObs].copy()
                else:
                    yield adata[batchObs]

    elif axis in [1, "feature"]:
        assert batchKey in adata.var.columns, f"{batchKey} not detected in adata"
        indexName = "index" if (not adata.var.index.name) else adata.var.index.name
        adata.var["__group"] = adata.var[batchKey]
        batchVarLs = (
            adata.var.filter(["__group"])
            .reset_index()
            .groupby("__group")[indexName]
            .agg(list)
        )
        del adata.var["__group"]
        for batchVar in tqdm(batchVarLs, disable=disableBar):
            if needName:
                if copy:
                    yield adata[batchVar].var.iloc[0].loc["__group"], adata[
                        batchVar
                    ].copy()
                else:
                    yield adata[batchVar].var.iloc[0].loc["__group"], adata[batchVar]
            else:
                if copy:
                    yield adata[:, batchVar].copy()
                else:
                    yield adata[:, batchVar]

    else:
        assert False, "Unknown `axis` parameter"


def simplify_adata_legacy(
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