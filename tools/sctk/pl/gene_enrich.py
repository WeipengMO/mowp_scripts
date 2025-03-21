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
        gene_list_sort_key: str = None,
        ax = None,
        figsize = (4, 4),
        fontsize = 8,
        fontweight = 'bold',
        dotsize = 15,
        up_color: str = '#d62728',
        down_color: str = '#1f77b4',
        rest_color: str = 'lightgrey',
        set_pvalue: bool = True,
        control_label: str = None,
        treatment_label: str = None,
    ):
    """
    Generate a volcano plot.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the results.
    pvalue_key : str
        Column name for p-values.
    log2fc_key : str
        Column name for log2 fold changes.
    pvalue_threshold : float
        Threshold for p-values.
    log2fc_threshold : float
        Threshold for log2 fold changes.
    n_top : int
        Number of top genes to label.
    gene_list : list, optional
        List of specific genes to label. If None, top genes are labeled.
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on.
    figsize : tuple, optional
        Size of the figure.
    fontsize : int, optional
        Font size for gene labels.
    fontweight : str, optional
        Font weight for gene labels.
    up_color : str, optional
        Color for upregulated genes.
    down_color : str, optional
        Color for downregulated genes.
    rest_color : str, optional
        Color for the rest of the genes.
    set_pvalue : bool, optional
        Whether to set p-values of 0 to the minimum non-zero p-value.
    control_label : str, optional
        Label for the control group.
    treatment_label : str, optional
        Label for the treatment group.

    Examples
    --------
    >>> sctk.pl.plot_volcano(df, log2fc_key='log2FoldChange', pvalue_key='padj', n_top=10)
    >>> sctk.pl.plot_volcano(df, log2fc_key='log2FoldChange', pvalue_key='padj', gene_list=['ISG15', 'IFIT3', 'LY6E'])
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
        s=dotsize, linewidth=0, ax=ax, color=rest_color)

    sns.scatterplot(
        data = df_up, x = log2fc_key, y = 'nlog10',
        s=dotsize, linewidth=0, ax=ax, color=up_color)

    sns.scatterplot(
        data = df_down, x = log2fc_key, y = 'nlog10',
        s=dotsize, linewidth=0, ax=ax, color=down_color)


    if gene_list is None and gene_list_sort_key is not None:
        _df = df[(df[pvalue_key] < pvalue_threshold) & (abs(df[log2fc_key]) > log2fc_threshold)]
        _df_sorted = _df.sort_values(gene_list_sort_key)
        if gene_list_sort_key != pvalue_key:
            gene_list = list(_df_sorted.index[:n_top//2]) + list(_df_sorted.index[-n_top//2:])
        else:
            gene_list = list(_df_sorted[_df_sorted[log2fc_key] > 0].index[:n_top//2]) + list(_df_sorted[_df_sorted[log2fc_key] < 0].index[:n_top//2])
    else:
        logger.info('gene_list_sort_key is None, no gene_list is provided')
        
    if isinstance(gene_list, list) and len(gene_list) > 0:
        texts = []
        df['log2FoldChange'] = df[log2fc_key]
        for item in df[df.index.isin(gene_list)].itertuples():
            texts.append(
                ax.text(
                    x = item.log2FoldChange, y = item.nlog10, s = item.Index,
                    ha='center', va='center',
                    fontsize=fontsize, weight = fontweight))
            
        adjust_text(texts, arrowprops = dict(arrowstyle = '-', color = 'k'), ax=ax)


    ax.axhline(log_pvalue_threshold, zorder = 1, c = 'k', lw = 1, ls = '--')
    ax.axvline(log2fc_threshold, zorder = 1, c = 'k', lw = 1, ls = '--')
    ax.axvline(-log2fc_threshold, zorder = 1, c = 'k', lw = 1, ls = '--')

    if control_label is not None and treatment_label is not None:
        ylim = plt.ylim()
        xlim = plt.xlim()
        arrow_control_x = min(-xlim[0], xlim[1]) / 2

        plt.annotate(control_label, 
             xy=(-log2fc_threshold, ylim[1]), 
             xytext=(-arrow_control_x, ylim[1]),
             arrowprops=dict(edgecolor='k', arrowstyle='<-'),
             horizontalalignment='right', verticalalignment='center')

        plt.annotate(treatment_label, 
                    xy=(log2fc_threshold, ylim[1]), 
                    xytext=(arrow_control_x, ylim[1]),
                    arrowprops=dict(edgecolor='k', arrowstyle='<-'),
                    horizontalalignment='left', verticalalignment='center')

    plt.xlabel("$\mathrm{Log_{2}}$(fold change)")
    plt.ylabel("-$\mathrm{Log_{10}}$(p-value)")
    sns.despine()

    return ax


def plot_gsea_bar(
        enr: pd.DataFrame, 
        pvalue_key='FDR p-value', 
        pvalue_threshold=0.05,
        both_directions=True,
        n_top=10,
        fontsize=8,
        figsize = None,
        title: str = '',
        bar_color: str = 'lightblue'):
    
    enr = enr[enr[pvalue_key] < pvalue_threshold].copy()
    if both_directions:
        enr1 = enr[enr['NES'] > 0].sort_values('NES', ascending=False).head(n_top//2)
        enr2 = enr[enr['NES'] < 0].sort_values('NES', ascending=True).head(n_top//2)
        top_results = pd.concat([enr1, enr2])
    else:
        top_results = enr.sort_values('NES', ascending=False).head(n_top)
    
    if figsize is None:
        figsize = (6, 0.5 * len(top_results))
    elif isinstance(figsize, tuple) and len(figsize) == 2:
        pass
    else:
        figsize = (figsize, 0.6 * len(top_results))

    plt.figure(figsize=figsize)
    plt.barh(top_results['Term'], top_results['NES'], color=bar_color)
    plt.xlabel('Normalized Enrichment Score (NES)', fontsize=12)
    plt.title(title)
    ylim = plt.ylim()
    plt.ylim(ylim[::-1])

    ax = plt.gca()
    for text, nes in zip(ax.get_yticklabels(), top_results['NES']):
        ha = 'left' if nes > 0 else 'right'
        plt.text(
            text.get_position()[0], text.get_position()[1], text.get_text(), 
            fontsize=fontsize, ha=ha, va='center')
    ax.set_yticks([])
    ax.set_ylabel('')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    return ax


def plot_ora_bar(
        enr: pd.DataFrame, 
        pvalue_key='FDR p-value', 
        pvalue_threshold=0.05,
        both_directions=True,
        n_top=10,
        fontsize=8,
        figsize = None,
        title: str = '',
        bar_color: str = 'lightblue'):
    
    enr = enr[enr[pvalue_key] < pvalue_threshold].copy()
    # remove duplicates item
    enr.sort_values(pvalue_key, ascending=True, inplace=True)
    enr.drop_duplicates('Term', keep='first', inplace=True)
    
    if both_directions:
        enr1 = enr[enr['change'] == 'up'].sort_values('Odds ratio', ascending=False).head(n_top//2)
        enr2 = enr[enr['change'] == 'down'].sort_values('Odds ratio', ascending=False).head(n_top//2)
        enr2['Odds ratio'] = -enr2['Odds ratio']
        top_results = pd.concat([enr1, enr2])
    else:
        top_results = enr.sort_values('NES', ascending=False).head(n_top)
    
    if figsize is None:
        figsize = (6, 0.3 * len(top_results))
    elif isinstance(figsize, tuple) and len(figsize) == 2:
        pass
    else:
        figsize = (figsize, 0.3 * len(top_results))

    plt.figure(figsize=figsize)
    plt.barh(top_results['Term'], top_results['Odds ratio'], color=bar_color)
    plt.xlabel('Odds ratio', fontsize=12)
    plt.title(title)
    ylim = plt.ylim()
    plt.ylim(ylim[::-1])

    ax = plt.gca()
    for text, nes in zip(ax.get_yticklabels(), top_results['Odds ratio']):
        ha = 'left' if nes > 0 else 'right'
        plt.text(
            text.get_position()[0], text.get_position()[1], text.get_text(), 
            fontsize=fontsize, ha=ha, va='center')
    ax.set_yticks([])
    ax.set_ylabel('')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    return ax