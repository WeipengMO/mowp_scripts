from typing import Optional
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Optional
from loguru import logger
import numpy as np
from ._color import ColorPalette
from ..pp import adata_process
from ..utils import configure_logger


def percent_in_cluster(
    adata,
    groupby: str,
    label: str,
    labelColor: Optional[dict] = None,
    needCounts: bool = True,
    ax=None,
    plot_legend: bool = True,
    legend_kwargs: dict ={"bbox_to_anchor": [1, 1]},
    swapAxes: bool = False,
    sort_index: bool = False,
    index_order: list = None,
    label_order: list = None,
):
    """
    Plot the percentage of each label in each groupby category.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        The column name of the category.
    label
        The column name of the label.
    labelColor
        The color of each label.
    needCounts
        Whether to show the counts of each groupby category.
    ax
        The axes object to plot on.
    plot_legend
        Whether to plot the legend.
    legend_kwargs
        The keyword arguments for the legend.
    swapAxes
        Whether to swap the x-axis and y-axis.
    sort_index
        Whether to sort the index.
    index_order
        The order of the index.
    label_order:
        The order of the labels.

    Returns
    -------
    The axes object.

    Examples
    --------
    >>> sctk.pl.percent_in_cluster(adata, 'leiden', 'batch')
    """
    if not labelColor:
        labelColor = adata_process.get_adata_color(adata, label)

    groupbyWithLabelCountsDf = (
        adata.obs.groupby(groupby)[label].apply(lambda x: x.value_counts()).unstack()
    )

    if sort_index:
        index = sorted(groupbyWithLabelCountsDf.index)
        groupbyWithLabelCountsDf = groupbyWithLabelCountsDf.reindex(index)
    elif index_order is not None:
        groupbyWithLabelCountsDf = groupbyWithLabelCountsDf.reindex(index_order)

    if label_order is not None:
        groupbyWithLabelCountsDf = groupbyWithLabelCountsDf[label_order]
    
    groupbyWithLabelCounts_CumsumPercDf = groupbyWithLabelCountsDf.pipe(
        lambda x: x.cumsum(1).div(x.sum(1), 0) * 100
    )
    legendHandleLs = []
    legendLabelLs = []
    if not swapAxes:
        for singleLabel in groupbyWithLabelCounts_CumsumPercDf.columns[::-1]:
            ax = sns.barplot(
                x=groupbyWithLabelCounts_CumsumPercDf.index,
                y=groupbyWithLabelCounts_CumsumPercDf[singleLabel],
                color=labelColor[singleLabel],
                ax=ax,
            )
            plt.sca(ax)
            legendHandleLs.append(
                plt.Rectangle(
                    (0, 0), 1, 1, fc=labelColor[singleLabel], edgecolor="none", label=singleLabel
                )
            )
            legendLabelLs.append(singleLabel)
        legendHandleLs, legendLabelLs = legendHandleLs[::-1], legendLabelLs[::-1]
        # plt.legend(legendHandleLs, legendLabelLs, frameon=False, **legend_kwargs)
        if plot_legend:
            plt.legend(handles=legendHandleLs, frameon=False, **legend_kwargs)

        plt.xlabel(groupby.capitalize())
        plt.ylabel(f"Percentage")
        if needCounts:
            for i, label in enumerate(groupbyWithLabelCounts_CumsumPercDf.index):
                plt.text(
                    i,
                    105,
                    f"$\it{{N}}$ = {adata[adata.obs[groupby] == label].shape[0]}",
                    rotation=90,
                    ha="center",
                    va="bottom",
                )
    else:
        for singleLabel in groupbyWithLabelCounts_CumsumPercDf.columns[::-1]:
            ax = sns.barplot(
                y=groupbyWithLabelCounts_CumsumPercDf.index,
                x=groupbyWithLabelCounts_CumsumPercDf[singleLabel],
                color=labelColor[singleLabel],
                ax=ax,
            )
            plt.sca(ax)
            legendHandleLs.append(
                plt.Rectangle(
                    (0, 0), 1, 1, fc=labelColor[singleLabel], edgecolor="none", label=singleLabel
                )
            )
            legendLabelLs.append(singleLabel)

        if plot_legend:
            plt.legend(handles=legendHandleLs[::-1], frameon=False, **legend_kwargs)

        plt.ylabel(groupby.capitalize())
        plt.xlabel(f"Percentage")
        if needCounts:
            for i, label in enumerate(groupbyWithLabelCounts_CumsumPercDf.index):
                plt.text(
                    101,
                    i,
                    f"$\it{{N}}$ = {adata[adata.obs[groupby] == label].shape[0]}",
                    rotation=0,
                    ha="left",
                    va="center",
                )
                
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    return ax


def _get_group_data(
    adata, groupby: str, key: str, 
    subgroups: Dict[str, List[str]] = None, remove_suffix: bool = False):
    '''
    Get the percentage of each cell type in each group.

    Parameters
    ----------
    adata: AnnData
        Annotated data object containing the cell data.
    groupby: str
        The key of the observation grouping to consider. eg. "batch"
    key: List[str]
        Keys for accessing fields of `.obs`. eg. "cell_type"
    subgroups: Dict[str, List[str]]
        Dictionary mapping group names to a list of cell types.
    remove_suffix: bool
        Whether to remove the suffix from the index.
    '''
    
    df = adata.obs.groupby(groupby)[key].apply(lambda x: x.value_counts()).unstack()
    df = df.div(df.sum(axis=1), axis=0) * 100

    df_group = {}
    if subgroups is not None:
        assert len(subgroups) == 2, 'subgroups must have 2 groups'
        for group, sub in subgroups.items():
            df_group[group] = df.loc[sub, :].copy()
            if remove_suffix:
                # TODO: This hack to remove the suffix from the index not work for all cases
                df_group[group].index = df_group[group].index.map(lambda x: x.split('_')[0])
    else:
        df_group['all'] = df
    
    return df_group


def _reorder_data(df_group: dict, keys_to_order: list = None):
    if keys_to_order is None:
        return

    if isinstance(keys_to_order, str):
        keys_to_order = [keys_to_order]

    for key in df_group:
        df = df_group[key].loc[:, keys_to_order].copy()
        df['row_sum'] = df.sum(axis=1)
        df = df.sort_values('row_sum', ascending=False)
        index_order = df.index
        df_group[key] = df_group[key].reindex(index_order)


def _plot_grouped_bars(adata, df_group, subkeys, key, figsize, subtitle=True, debug=True):
    """
    Plot grouped bar charts with legend outside the figure.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    df_group : dict of DataFrame
        Dictionary where keys are group names and values are DataFrames for plotting.
    subkeys : dict
        Dictionary where keys are group indices and values are lists of cell types.
    key : str, optional
        Key to use for cell types in `adata.obs` .
    groupby : str, optional
        Key to group by in `adata.obs`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The resulting figure object.
    axes : numpy.ndarray of matplotlib.axes._subplots.AxesSubplot
        Array of axes objects.
    """

    if debug:
        configure_logger(log_level="debug")
    else:
        configure_logger(log_level="info")


    label_color_key = f'{key}_colors'
    if label_color_key in adata.uns:
        # Get color mapping from adata.uns
        cell_type_colors = dict(zip(adata.obs[key].cat.categories, adata.uns[label_color_key]))
    else:
        cell_type_colors = dict(zip(
            adata.obs[key].cat.categories, 
            ColorPalette(len(adata.obs[key].cat.categories), 'tab20').to_array()))

    # Determine width ratios for subplots
    width_ratios = [len(df) for df in df_group.values()]

    fig, axes = plt.subplots(
        nrows=1, ncols=len(width_ratios),
        figsize=figsize,
        gridspec_kw={'width_ratios': width_ratios},
        sharey=True
    )

    logger.debug(f"{width_ratios=}")
    
    if not isinstance(axes, np.ndarray):
        axes = [axes]

    is_first_axes = True
    for ax, (batch, df) in zip(axes, df_group.items()):
        if not is_first_axes:
            ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        is_first_axes = False
        
        for i, cell_type_list in enumerate(subkeys.values()):
            plot_df = df.loc[:, cell_type_list]
            colors = [cell_type_colors[cell_type] for cell_type in plot_df.columns]

            if i == 1:
                plot_df *= -1
            
            plot_df.plot(kind='bar', stacked=True, ax=ax, legend=False, color=colors)
            ax.axhline(0, ls='--', color='k')
            if subtitle:
                ax.set_title(batch)
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=True)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.set_xlabel('')

            sns.despine(left=True, bottom=True)
    
    axes[0].set_ylabel('Percentage')

    # Create legend
    handles, labels = [], []
    for cell_type, color in cell_type_colors.items():
        handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10))
        labels.append(cell_type)

    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.1, .9), frameon=False)

    fig.subplots_adjust(wspace=0.1)
    plt.show()


def percent_in_cluster_group(
    adata,
    groupby: str,
    key: str,
    subgroups: Dict[str, List[str]] = None,
    subkeys: Dict[str, List[str]] = None,
    label_order: str = None,
    subtitle: bool = True,
    figsize: tuple = (8, 4),
    remove_suffix: bool = False,
):
    """
    Calculate the percentage of cells in each cluster group.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    groupby : str
        Key for the column in `adata.obs` that contains the sample labels.
    key : str
        Key for the column in `adata.obs` that contains the cell labels.
    subgroups : dict
        Dictionary mapping sample group names to a list of sample labels.
    subkeys : dict
        Dictionary mapping cell group names to a list of subkeys.
    label_order : str
        The order in which the cluster groups should be plotted.
    subtitle : bool
        Whether to add a subtitle to each subplot.
    figsize : tuple
        Size of the figure.
    remove_suffix : bool
        Whether to remove the suffix from the index.
    """
    # Get the percentage of each cell type in each group
    df_group = _get_group_data(adata, groupby, key, subgroups, remove_suffix=remove_suffix)

    # Reorder the data based on the label_order
    # If subkeys is None, then label_order must be None, and ordering is by string
    if subkeys is None:
        subkeys = {key: adata.obs[key].cat.categories}
        if label_order is not None:
            raise ValueError('label_order must be None when subkeys is None')
        keys_to_order = None
    else:
        keys_to_order = subkeys[label_order] if label_order is not None else None

    _reorder_data(df_group, keys_to_order)

    # Plot the grouped bar charts
    _plot_grouped_bars(
        adata, df_group, subkeys, key=key, figsize=figsize, subtitle=subtitle)