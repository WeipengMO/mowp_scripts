import decoupler as dc
import pandas as pd
import numpy as np
from typing import List, Union
from PyComplexHeatmap import DotClustermapPlotter
from matplotlib import pyplot as plt


def _run_ora(
        genes: list, 
        gene_sets: pd.DataFrame,
        source: str = 'geneset',
        target: str = 'genesymbol',
        threshold: float =0.05):
    
    enr_pvals = dc.get_ora_df(
        df=genes,
        net=gene_sets,
        source=source,
        target=target
    ).sort_values('FDR p-value', ascending=True)

    enr_pvals.query('`FDR p-value` < @threshold', inplace=True)

    return enr_pvals


def run_ora(
        df: Union[pd.DataFrame, List[str]],
        gene_sets: pd.DataFrame,
        source: str = 'geneset',
        target: str = 'genesymbol',
        threshold: float = 0.05):
    
    results = []
    if isinstance(df, pd.DataFrame):
        for col in df.columns:
            genes = list(df[col].values)
            enr_pvals = _run_ora(genes, gene_sets, source, target, threshold=threshold)
            enr_pvals['Name'] = col
            results.append(enr_pvals)
    elif isinstance(df, list):
        enr_pvals = _run_ora(df, gene_sets, source, target, threshold=threshold)
        enr_pvals['Name'] = 'Gene list'
        results.append(enr_pvals)

    results = pd.concat(results)
    results['-log(FDR p-value)'] = -np.log10(results['FDR p-value'])
    results['log(Odds ratio)'] = np.log10(results['Odds ratio'])
    results['log(Combined score)'] = np.log10(results['Combined score'])

    return results


def ora_dotplot(
        results: pd.DataFrame,
        x: str = 'Name',
        y: str = 'Term',
        score: str = 'log(Odds ratio)',
        pvalue: str = '-log(FDR p-value)',
        smax: float = 95,
        smin: float = 0,
        cmap: str = 'Reds',
        row_cluster: bool = True,
        col_cluster: bool = False,
        n_top: int = None,
        figsize: tuple = (3, 6),
        **kwargs):
    
    _results = results.copy()
    if n_top is not None and n_top > 0:
        _terms = _results.groupby(x).apply(
            lambda x: x.sort_values(pvalue, ascending=False).head(10)
            ).reset_index(drop=True)[y].unique()
       
        _results = _results[results[y].isin(_terms)]

    smax = np.percentile(_results[pvalue], smax)
    smin = np.percentile(_results[pvalue], smin)
    _results[pvalue] = np.clip(_results[pvalue], smin, smax)
    _results[pvalue] = np.round(_results[pvalue])

    plt.figure(figsize=figsize)
    cm = DotClustermapPlotter(
        data=_results, 
        x=x, y=y, 
        value=score, c=score, s=pvalue,
        cmap=cmap, 
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        show_rownames=True,
        show_colnames=True,
        yticklabels_kws={'labelsize': 8},
        xticklabels_kws={'labelrotation': 30},
        legend_gap=10,
        **kwargs,
    )
    plt.show()