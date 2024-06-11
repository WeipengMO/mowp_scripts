from matplotlib import pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from ._color import ColorPalette
from typing import Union
import scanpy as sc


def adjust_heatmap_labels(ax, gene_to_show: set):
    '''
    Adjust the labels of the heatmap to show only the genes in gene_to_show

    Note
    ----
    This function is not perfect, it may not work well in some cases.
    Only test in the following case:
    ```
    ax = sc.pl.rank_genes_groups_heatmap(
        adata, n_genes=20, show_gene_labels=True, 
        key='rank_genes_groups_filtered', cmap='Reds', swap_axes=True, dendrogram=False, standard_scale='var', 
        figsize=(6, 8), show=False)
    ```
    '''
    _ax = ax['heatmap_ax']
    plt.sca(_ax)
    yticklabels = []
    xlim = _ax.get_xlim()
    for text in _ax.get_yticklabels():
        if text.get_text() in gene_to_show:
            yticklabels.append(
                plt.text(-xlim[1]/200, text.get_position()[1], text.get_text(), fontsize=8, ha='right'))
    _ax.set_yticks([])

    adjust_text(
        yticklabels, 
        ensure_inside_axes = False,
        arrowprops=dict(arrowstyle='->'),
        only_move = {'explode': 'y', 'pull': 'y', 'static': 'y', 'text': 'y'})


def adjust_groupby_ax(ax):
    '''
    Adjust the labels of the groupby_ax to show the labels in the heatmap

    Note
    ----
    This function is not perfect, it may not work well in some cases.
    Only test in the following case:
    ```
    ax = sc.pl.rank_genes_groups_heatmap(
        adata, n_genes=20, show_gene_labels=True, 
        key='rank_genes_groups_filtered', cmap='Reds', swap_axes=True, dendrogram=False, standard_scale='var', 
        figsize=(6, 8), show=False)
    ```
    '''
    _ax = ax['groupby_ax']
    plt.sca(_ax)
    xticklabels = []
    ylim = _ax.get_ylim()
    for text in _ax.get_xticklabels():
        xticklabels.append(
            plt.text(text.get_position()[0], ylim[0]*3, text.get_text(), fontsize=8))

    _ax.set_xticks([])
    _ax.set_xlabel('')

    adjust_text(
        xticklabels, 
        ensure_inside_axes=False,
        arrowprops=dict(arrowstyle='-'),
        only_move = {'explode': 'x', 'pull': 'x', 'static': 'x', 'text': 'x'},
        )


def set_adata_color(adata, key: str, color_dict: dict = None, palette: Union[list, str] = 'tab20', inplace=True):
    '''
    Set the color of adata.obs[key] based on color_dict or palette

    Parameters
    ----------
    adata
        An AnnData object.
    key
        The `adata.obs` key of the observation grouping to consider.
    color_dict
        A dictionary mapping categories to colors.
    palette
        A string or list of colors.
    inplace
        Whether to set the color in adata.uns or return the color list.
    '''
    obs_categories = adata.obs[key].astype('category').cat.categories

    if color_dict is None:
        if isinstance(palette, str):
            color_dict = dict(zip(
                obs_categories, 
                ColorPalette(len(adata.obs[key].cat.categories), palette).to_array()))
        elif isinstance(palette, list):
            color_dict = dict(zip(obs_categories, palette))

    diff_set = set(obs_categories).difference(set(color_dict))
    if len(diff_set) > 0:
        raise ValueError(f'color_dict does not contain all categories in {key}')

    if inplace:
        adata.uns[key + '_colors'] = [color_dict[x] for x in obs_categories]
    else:
        return [color_dict[x] for x in obs_categories]


def boxplot(
        adata: sc.AnnData, 
        keys: Union[str, list[str]], 
        groupby: str, 
        showfliers=False,
        layer=None,
        use_raw=False,
        palette=None,
        figsize=(4, 4),
        sns_kwargs={}):
    '''
    Boxplot of adata.obs[key] grouped by adata.obs[groupby]

    Parameters
    ----------
    adata
        An AnnData object.
    keys
        Keys for accessing variables of `.var_names` or fields of `.obs`.
    groupby
        The `adata.obs` key of the observation grouping to consider.
    showfliers
        Show the outliers.
    layer
        The layer to use.
    use_raw
        Use `adata.raw` for plotting.
    palette
        A list of colors or a color palette.
    figsize
        The size of the subfigure.
    sns_kwargs
        Other keyword arguments for seaborn.boxplot.
    '''
    from collections import OrderedDict


    if isinstance(keys, str):
        keys = [keys]
    keys = list(OrderedDict.fromkeys(keys))  # remove duplicates, preserving the order

    obs_df = sc.get.obs_df(adata, keys=[groupby] + keys, layer=layer, use_raw=use_raw)
    
    fig, axes = plt.subplots(
        nrows=1, ncols=len(keys),
        figsize=(figsize[0]*len(keys), figsize[1]),
    )

    if palette is None:
        uns_colors = adata.uns.get(groupby + '_colors')
        if uns_colors is not None:
            palette = uns_colors

    for i, key in enumerate(keys):
        if len(keys) == 1:
            ax = axes
        else:
            ax = axes[i]
        sns.boxplot(
            x=groupby, y=key,
            data=obs_df,
            ax=ax,
            linewidth=1,
            fliersize=0,
            showfliers=showfliers,
            palette=palette,
            **sns_kwargs,
        )
        ax.set_title(key)
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        sns.despine(ax=ax)

    plt.show()
    