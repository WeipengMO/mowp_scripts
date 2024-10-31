import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
from loguru import logger


def plot_go_results(
        df: pd.DataFrame, 
        pvalue: str, 
        term: str,
        top_n: int = 10, 
        ax = None, 
        title: str = None, 
        color: str = None,
        figsize: tuple = (5, 5),
        fontsize: int = 8,
        xlabel: str = 'log(pvalue)',
        **kwargs):
    """
    Plot the top GO results.

    Parameters
    ----------
    df
        The GO results dataframe.
    pvalue
        The key of the pvalue column.
    term
        The key of the term column.
    top_n
        The number of top GO terms to plot.
    ax
        The axes object to plot on.
    title
        The title of the plot.
    color
        The color of the bars.
    
    Returns
    -------
    ax
        The axes object.
    """
    df = df.copy()
    df['-log10(pvalue)'] = -np.log10(df[pvalue])
    if term == 'index':
        df['term'] = df.index
        term = 'term'

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(
        data=df.head(top_n), x='-log10(pvalue)', y=term, 
        color=color, 
        ax=ax,
        **kwargs)
    
    for text in ax.get_yticklabels():
        plt.text(
            text.get_position()[0], 
            text.get_position()[1], 
            text.get_text(), 
            fontsize=fontsize, 
            ha='left', 
            va='center'
        )
        
    ax.set_yticks([])

    ax.set_xlabel(xlabel)
    ax.set_ylabel('')

    sns.despine(ax=ax)
    if title is not None:
        ax.set_title(title)

    return ax


def plot_volcano(
        df: pd.DataFrame,
        pvalue_key: str = 'padj',
        log2fc_key: str = 'log2FoldChange',
        pvalue_threshold: float = 0.01,
        log2fc_threshold: float = 1,
        n_top: int = 10,
        gene_list: list = None,
        ax = None,
        figsize = (4, 4),
        fontsize = 8,
        up_color: str = '#d62728',
        down_color: str = '#1f77b4',
        rest_color: str = 'lightgrey',
        set_pvalue: bool = True):
    """
    Plot a volcano plot.

    Parameters
    ----------
    df
        The dataframe containing the results.
    pvalue_key
        The key of the pvalue column.
    log2fc_key
        The key of the log2 fold change column.
    pvalue_threshold
        The pvalue threshold.
    log2fc_threshold
        The log2 fold change threshold.
    n_top
        The number of top genes to label.
    gene_list
        A list of genes to label. If None, the top genes will be labeled.
    ax
        The axes object to plot on.
    figsize
        The figure size.
    up_color
        The color of the upregulated genes.
    down_color
        The color of the downregulated genes.
    rest_color
        The color of the rest genes.
    
    Returns
    -------
    ax
        The axes object.

    Examples
    --------
    >>> sctk.pl.plot_volcano(df, log2fc_key='log2FoldChange', pvalue_key='padj', n_top=10)
    >>> sctk.pl.plot_volcano(de, log2fc_key='log2FoldChange', pvalue_key='padj', gene_list=['ISG15', 'IFIT3', 'LY6E'])
    """
    
    assert pvalue_key in df.columns, f'pvalue_key {pvalue_key} not in df.columns'
    assert log2fc_key in df.columns, f'log2fc_key {log2fc_key} not in df.columns'

    if ax is None:
        fig, ax = plt.subplots(figsize = figsize)

    df = df.copy()
    if set_pvalue:
        logger.info('Setting pvalue to min_pvalue for pvalue == 0')
        min_pvalue = df[pvalue_key][df[pvalue_key] > 0].min()
        df.loc[df[pvalue_key] == 0, pvalue_key] = min_pvalue


    df['nlog10'] = -np.log10(df[pvalue_key])
    log_pvalue_threshold = -np.log10(pvalue_threshold)
    df_up = df[(df[log2fc_key] > log2fc_threshold) & (df.nlog10 > log_pvalue_threshold)]
    df_down = df[(df[log2fc_key] < -log2fc_threshold) & (df.nlog10 > log_pvalue_threshold)]
    df_rest = df[(df[log2fc_key].abs() <= log2fc_threshold) | (df.nlog10 <= log_pvalue_threshold)]

    sns.scatterplot(
        data = df_rest, x = log2fc_key, y = 'nlog10',
        s=15, linewidth=0, ax=ax, color=rest_color)

    sns.scatterplot(
        data = df_up, x = log2fc_key, y = 'nlog10',
        s=15, linewidth=0, ax=ax, color=up_color)

    sns.scatterplot(
        data = df_down, x = log2fc_key, y = 'nlog10',
        s=15, linewidth=0, ax=ax, color=down_color)


    if gene_list is None:
        _df = df[(df[pvalue_key] < pvalue_threshold) & (abs(df[log2fc_key]) > log2fc_threshold)]
        if 'stat' in _df.columns:
            _df_sorted = _df.sort_values('stat')
            gene_list = list(_df_sorted.index[:n_top//2]) + list(_df_sorted.index[-n_top//2:])
        else:
            _df_sorted = _df.sort_values(pvalue_key)
            gene_list = list(_df_sorted.index[:n_top])
        
    if len(gene_list) > 0:
        texts = []
        df['log2FoldChange'] = df[log2fc_key]
        for item in df[df.index.isin(gene_list)].itertuples():
            texts.append(
                ax.text(
                    x = item.log2FoldChange, y = item.nlog10, s = item.Index,
                    ha='center', va='center',
                    fontsize=fontsize, weight = 'bold'))
            
        adjust_text(texts, arrowprops = dict(arrowstyle = '-', color = 'k'), ax=ax)


    ax.axhline(log_pvalue_threshold, zorder = 1, c = 'k', lw = 1, ls = '--')
    ax.axvline(log2fc_threshold, zorder = 1, c = 'k', lw = 1, ls = '--')
    ax.axvline(-log2fc_threshold, zorder = 1, c = 'k', lw = 1, ls = '--')

    plt.xlabel("$\mathrm{log_{2}}$(fold change)")
    plt.ylabel("-$\mathrm{log_{10}}$(p-value)")
    sns.despine()

    return ax