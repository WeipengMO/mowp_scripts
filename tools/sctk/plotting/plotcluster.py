from typing import Optional
from ..preprocessing import adata_process
import matplotlib.pyplot as plt
import seaborn as sns


def percent_in_cluster(
    adata,
    groupby: str,
    label: str,
    labelColor: Optional[dict] = None,
    needCounts=True,
    ax=None,
    plot_legend=True,
    legend_kwargs={"bbox_to_anchor": [1, 1]},
    swapAxes=False,
    sort_index=False,
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


def percent_in_cluster_group(
        adata, 
        groupby: str, 
        label: str, 
        group_dict: dict, 
        figsize: tuple = (8, 3),
        sort_index: bool = True):

    """
    Plot the percentage of each label in each groupby category.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        The column name of the category. eg. 'batch'
    label
        The column name of the label. eg. 'cell_type'
    group_dict
        A dictionary with group name as key and group list as value. eg. {'group1': ['label1', 'label2'], 'group2': ['label3', 'label4']}
    figsize
        The figure size.
    sort_index
        Whether to sort the index.
    
    Examples
    --------
    >>> cd45_pos = [
            'B cells', 'NK cells', 'T cells', 'TAM', 'macrophages', 'mast cells',
            'neutrophils', 'pDCs', 'plasma cells']
        cd45_neg = [
            'endothelial cells', 'fibroblasts', 'maligant cells']
        group_dict = {
            'CD45+': cd45_pos,
            'CD45-': cd45_neg}
    >>> sctk.pl.percent_in_cluster_group(adata, 'batch', 'cell_type', group_dict)
    """

    assert len(group_dict) == 2, 'group_dict must have 2 groups'
    
    # A dataframe with groupby as index and label as columns, values are counts
    groupby_label_counts_df = (
        adata.obs.groupby(groupby)[label].apply(lambda x: x.value_counts()).unstack()
    )

    label_color = f'{label}_colors'
    if label_color in adata.uns:
        # get the same color set from adata.uns
        label_color = dict(zip(adata.obs[label].cat.categories, adata.uns[label_color]))
    else:
        label_color = None

    if sort_index:
        index = sorted(groupby_label_counts_df.index)[::-1]
        groupby_label_counts_df = groupby_label_counts_df.reindex(index)

    fig, ax = plt.subplots(1, 2, figsize=figsize, sharey=True)

    for i, group_name in enumerate(group_dict):
        _group = group_dict[group_name]
        color = [label_color[x] for x in _group]
        
        data = groupby_label_counts_df.loc[:, _group]
        data_percentage = data.div(data.sum(axis=1), axis=0) * 100
        data_percentage.plot(kind='barh', stacked=True, ax=ax[i], legend=True, color=color)
        
        ax[i].set_title(group_name)
        ax[i].set_ylabel('')
        ax[i].set_xlim(0, 100)
        ax[i].grid(False)
        ax[i].legend(loc='upper left', bbox_to_anchor=(0, -.2), ncol=2, frameon=False)

    fig.text(0.5, -.01, 'Percentage of cells', ha='center')

    plt.subplots_adjust(wspace=.1)
    
    return ax