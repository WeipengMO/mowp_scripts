import scanpy as sc
import memento
from tqdm.auto import tqdm
from loguru import logger
import pandas as pd
import numpy as np
import decoupler as dc
from typing import List, Union


def run_memento(
        ad: sc.AnnData,
        treatments: Union[str, List[str]],
        controls: Union[str, List[str]],
        condition_key: str,
        cell_type_key: str,
        cell_types: list = None,
        layer: str = 'counts',
        min_cells: int = 50,
        capture_rate: float = 0.25, #  0.07 for 10X v1, 0.15 for v2, and 0.25 for v3
        n_jobs: int = 32
):
    '''
    Run memento for a list of treatments and a control group.

    Parameters
    ----------
    ad : AnnData
        Annotated data object.
    treatments : list
        List of treatments to compare with control.
    controls : list
        Control group.
    condition_key : str
        Key in obs that contains the condition.
    cell_type_key : str
        Key in obs that contains the cell type.
    cell_types : list
        List of cell types to process. If None, all cell types are processed.
    layer : str
        Raw counts layer.
    min_cells : int
        Minimum number of cells to process a cell type.
    capture_rate : float
        Capture rate for the experiment.
    n_jobs : int
        Number of CPUs to use.
    '''
    result_dict = {}
    if cell_types is None:
        cell_types = list(ad.obs[cell_type_key].unique())
        logger.info(f'Cell type not specific, use cell type: {", ".join(cell_types)}')

    if isinstance(treatments, str):
        treatments = [treatments]
    
    if isinstance(controls, str):
        controls = [controls]

    for cell_type in tqdm(cell_types):
        for control in controls:
            for treatment in treatments:
                if control == treatment:
                    logger.warning(f'{control} and {treatment} are the same, skipping')
                    continue
                logger.info(f'Processing {cell_type}, {treatment} vs {control}')
                adata = ad[ad.obs[cell_type_key] == cell_type].copy()
                if len(adata) < min_cells:
                    logger.warning(f'Skipping {cell_type} due to insufficient cells')
                    continue

                value_counts = adata.obs[condition_key].value_counts()
                if value_counts[treatment] < min_cells or value_counts[control] < min_cells:
                    logger.warning(f'Skipping {cell_type} due to insufficient cells')
                    continue
                
                adata.X = adata.layers[layer]
                adata = adata[adata.obs[condition_key].isin([treatment, control])].copy()
                
                adata.obs[treatment] = adata.obs[condition_key].cat.rename_categories({control: 0, treatment: 1})

                result_1d = memento.binary_test_1d(
                    adata=adata, 
                    capture_rate=capture_rate, #  0.07 for 10X v1, 0.15 for v2, and 0.25 for v3
                    treatment_col=treatment, 
                    num_cpus=n_jobs,
                    verbose=0,
                    num_boot=5000)
                
                result_dict[f'{cell_type} | {treatment}_vs_{control}'] = result_1d
    
    return result_dict


def run_gsea(
        results_dict: dict, 
        gene_sets: pd.DataFrame, # msigdb
        stat_key: str = 'de_coef',
        source: str = 'geneset', 
        target: str = 'genesymbol', 
        enrich_dict: dict = None, 
        pvalue_threshold: str = None):

    return_results = False
    if enrich_dict is None:
        enrich_dict = {}
        return_results = True


    for res in results_dict:
        de = results_dict[res].set_index('gene')

        if pvalue_threshold is not None:
            de.query(pvalue_threshold, inplace=True)
            
        gsea_df = dc.get_gsea_df(
            df=de,
            stat=stat_key,
            net=gene_sets,
            source=source,
            target=target
        )

        gsea_df['cell_type'] = res.split(' | ')[0]
        gsea_df['compare'] = res.split(' | ')[1]
        gsea_df['change'] = gsea_df['NES'].map(lambda x: 'up' if x > 0 else 'down')

        enrich_dict[res] = gsea_df
    
    if return_results:
        return enrich_dict
    

def filter_gsea_res(enrich_dict, pvalue_key='FDR p-value', pvalue_threshold=0.05, n_top=10):
    concat_list = []
    term_list = []
    for res in enrich_dict:
        mask = enrich_dict[res][pvalue_key] < pvalue_threshold
        gsea_df = enrich_dict[res][mask].copy()

        term = list(gsea_df.query('change == "up"').sort_values(['NES'], ascending=False).head(n_top)['Term'])
        term_list.extend(term)

        term = list(gsea_df.query('change == "down"').sort_values(['NES'], ascending=True).head(n_top)['Term'])
        term_list.extend(term)

        concat_list.append(gsea_df)

    enrich_df = pd.concat(concat_list)
    term_list = list(set(term_list))
    enrich_df = enrich_df[enrich_df['Term'].isin(term_list)]
    
    return enrich_df


def run_ora(
        results_dict: dict, 
        gene_sets: pd.DataFrame, # msigdb
        foldchange_key = 'log2FoldChange',
        pvalue_key = 'pvalue',
        source: str ='geneset', 
        target: str ='genesymbol', 
        enrich_dict: dict = None, 
        fold_change_threshold: float = 1, 
        pvalue_threshold: float = 0.05):

    return_results = False
    if enrich_dict is None:
        enrich_dict = {}
        return_results = True

    for res in results_dict:
        de = results_dict[res].set_index('gene')

        enr_list = []
        for change in ('down', 'up'):
            if change == 'down':
                top_genes = de[(de[pvalue_key] < pvalue_threshold) & (de[foldchange_key] < -fold_change_threshold)]
            else:
                top_genes = de[(de[pvalue_key] < pvalue_threshold) & (de[foldchange_key] > fold_change_threshold)]

            # Run ora
            enr_pvals = dc.get_ora_df(
                df=top_genes,
                net=gene_sets,
                source=source,
                target=target
            )

            enr_pvals['cell_type'] = res.split(' | ')[0]
            enr_pvals['compare'] = res.split(' | ')[1]
            enr_pvals['change'] = change
            enr_list.append(enr_pvals)

        enrich_dict[res] = pd.concat(enr_list)
    
    if return_results:
        return enrich_dict
    

def filter_ora_res(enrich_dict, pvalue_key='FDR p-value', pvalue_threshold=0.01, n_top=5):
    concat_list = []
    term_list = []
    for res in enrich_dict:
        for fold_change in ('down', 'up'):
            mask = (enrich_dict[res][pvalue_key] < pvalue_threshold) & (enrich_dict[res]['change'] == fold_change)
            ora_df = enrich_dict[res][mask].copy()
            ora_df.sort_values(pvalue_key, ascending=True, inplace=True)

            term = ora_df.head(n_top)['Term']
            term_list.extend(term)

            concat_list.append(ora_df)

    enrich_df = pd.concat(concat_list)
    term_list = list(set(term_list))
    enrich_df = enrich_df[enrich_df['Term'].isin(term_list)]
    
    return enrich_df