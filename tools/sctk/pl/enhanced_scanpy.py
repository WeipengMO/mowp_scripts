from matplotlib import pyplot as plt
import seaborn as sns
from adjustText import adjust_text


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

