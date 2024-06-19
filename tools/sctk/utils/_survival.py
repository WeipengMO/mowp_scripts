from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt
import seaborn as sns
from typing import List, Union
from loguru import logger


class Survival(ad.AnnData):

    def __init__(
            self, 
            exp_data: ad.AnnData, 
            meta_data: pd.DataFrame, 
            meta_index_col: str,
            transpose_exp: bool = True, 
            meta_kwargs: dict = {},
            ):
        '''
        Initialize a Survival object.

        Parameters:
        ----------
        exp_data: ad.AnnData
            A gene expression matrix.

        meta_data: pd.DataFrame
            A DataFrame containing survival data.
        
        meta_index_col: str
            The column name to be used as the index of the metadata DataFrame.
        
        transpose_exp: bool
            Whether to transpose the expression data.

        meta_kwargs: dict
            Additional keyword arguments to be passed to pd.read_csv.
        
        '''
        adata = ad.read_csv(exp_data, delimiter='\t')
        if transpose_exp:
            adata = adata.T
        
        meta_data = pd.read_csv(meta_data, sep='\t', index_col=meta_index_col, **meta_kwargs)
        sample_ids_intersect = adata.obs.index.intersection(meta_data.index)
        if len(sample_ids_intersect) == 0:
            raise ValueError('No sample IDs are intersected between expression data and metadata.')

        adata = adata[sample_ids_intersect].copy()
        adata.obs = meta_data.loc[sample_ids_intersect]

        self._init_as_actual(
                X=adata.X,
                obs=adata.obs,
                var=adata.var,
            )
            

    def group_meta(
            self, 
            groupby: Union[str, List[str]], 
            group_method: Union[str, int] = 'median',
            event: str = 'OS', 
            time: str = 'OS.time',
            time_limit: float = None):
        '''
        Group samples based on a given gene signature or a list of genes.

        Parameters:
        ----------
        groupby: str | List[str]
            A gene signature or a list of genes to be used for grouping samples.
        
        group_method: str | int
            The method to group samples. If 'median' or 'quantile', samples will be grouped based on the median or quantiles of the gene signature. If an integer is given, samples will be grouped based on the quantiles of the gene signature.
        
        event
            The event column in the metadata DataFrame.
        
        time
            The time column in the metadata DataFrame.

        time_limit
            The time limit to be used for plotting the Kaplan-Meier curve.
        '''

        if isinstance(groupby, str):
            survival_data = sc.get.obs_df(self, [groupby, event, time])
        elif isinstance(groupby, list):
            # sc.tl.score_genes(self, groupby)
            # survival_data = sc.get.obs_df(self, ['score', event, time])

            common_genes = list(self.var_names.intersection(groupby))
            if len(common_genes) == 0:
                raise ValueError('No common genes are found between the gene signature and the expression data.')
            genes_to_ignore = set(groupby) - set(common_genes)
            if len(genes_to_ignore) > 0:
                logger.info(f'Genes to ignore: {genes_to_ignore}')

            # use average expression of a gene signature as a score
            survival_data = sc.get.obs_df(self, common_genes + [event, time])
            survival_data['score'] = survival_data[common_genes].mean(axis=1)
            groupby = 'score'
        else:
            raise ValueError('groupby must be a string or a list of strings.')
        
        survival_data.dropna(inplace=True)

        if time_limit is not None:
            survival_data = survival_data[survival_data[time] <= time_limit]
        
        if group_method == 'median': 
            survival_data['group'] = pd.qcut(survival_data[groupby], 2, labels=['Low', 'High'])
            self.group_label = ['Low', 'High']
        elif group_method == 'quantile':
            survival_data['group'] = pd.qcut(survival_data[groupby], 4, labels=['Low', 'q50', 'q75', 'High'])
            self.group_label = ['Low', 'High']
        elif isinstance(group_method, int):
            self.group_label = [f'G{i}' for i in range(group_method)]
            survival_data['group'] = pd.qcut(survival_data[groupby], group_method, labels=self.group_label)

        # return survival_data

        self.uns['survival'] = survival_data
        self.event = event
        self.time = time
        self.group_method = group_method
    

    def km_plot(
            self, 
            ci_show: bool = False,
            show_censors: bool = False,
            figsize=(3, 3),
            xlabel: str = 'Time',
            ylabel: str = 'Survival probability',
            pattle = None,
            ax = None
            ):
        if pattle is None:
            pattle = sns.color_palette(['#e41a1c', '#377eb8', '#984ea3', '#ff7f00'])

        kmf = KaplanMeierFitter()
        survival_data = self.uns['survival'][self.uns['survival']['group'].isin(self.group_label)].copy()
        survival_data['group'] = survival_data['group'].astype(str)
        
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        by_group = {}
        for (group, _df), color in zip(survival_data.groupby('group'), pattle):
            by_group[group] = _df
            kmf.fit(_df[self.time], _df[self.event], label=group)
            kmf.plot_survival_function(
                ci_show=ci_show, show_censors=show_censors, color=color, ax=ax)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False)
        sns.despine()

        if self.group_method in {'median', 'quantile'}:
            results = logrank_test(
                by_group['Low'][self.time], by_group['High'][self.time],
                by_group['Low'][self.event], by_group['High'][self.event])

            print(f'Group Low: {by_group["Low"].shape[0]}')
            print(f'Group High: {by_group["High"].shape[0]}')

            # results.print_summary()
            print(f'{results.p_value=}')
            print(f'{results.test_statistic=}')

        return ax
                