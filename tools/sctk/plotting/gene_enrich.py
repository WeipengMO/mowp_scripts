import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_go_results(
        df: pd.DataFrame, pvalue_key: str = 'p.adjust', 
        top: int = 10, ax = None, title: str = None, color: str = None):
    """
    Plot the top GO results.

    Parameters
    ----------
    df
        The GO results dataframe.
    pvalue_key
        The key of the pvalue column.
    top
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
    df['-log10(pvalue)'] = -np.log10(df[pvalue_key])
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    if color is None:
        color = '#3274a1'
    sns.barplot(
        data=df.head(top), x='-log10(pvalue)', y='Description', 
        color=color, ax=ax)

    ax.set_xlabel('-log10(pvalue)')
    ax.set_ylabel('')
    sns.despine(ax=ax)
    if title is not None:
        ax.set_title(title)

    return ax