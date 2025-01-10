import numpy as np
import pandas as pd
import scanpy as sc
from typing import Union, List
from loguru import logger


def create_pseudo_bulk(
        cell_subset: sc.AnnData, 
        sample_key: str, 
        condition_key: str,
        n_cells: int = 30,
        n_pseudo_replicates: int = 1,
        counts_key: str = 'counts'
        ) -> List[sc.AnnData]:
    """
    Create pseudo-bulk samples from a subset of cells.

    Parameters
    ----------
    cell_subset
        An AnnData object containing a subset of cells.
    sample_key
        The key in `cell_subset.obs` containing the sample name.
    condition_key
        The key in `cell_subset.obs` containing the condition name.
    n_cells
        The number of cells to sample per donor.
    n_pseudo_replicates
        The number of pseudo-replicates to create for each donor.
    counts_key
        The key in `cell_subset.layers` containing the counts matrix.

    Returns
    -------
    pb
        An AnnData objects containing a pseudo-bulk sample.

    Examples
    --------
    >>> cell_subset = adata[adata.obs['cell_type'] == 'B cells']
    >>> pb_samples = create_pseudo_bulk(cell_subset, 'sample', 'condition')
    """

    assert counts_key in cell_subset.layers.keys(), f"counts_key {counts_key} not found in cell_subset.layers.keys()"
    cell_subset.X = cell_subset.layers[counts_key]

    pbs = []
    for sample in cell_subset.obs[sample_key].unique():
        samp_cell_subset = cell_subset[cell_subset.obs[sample_key] == sample]
        
        indices = list(samp_cell_subset.obs_names)
        if len(indices) < n_cells * n_pseudo_replicates:
            logger.info(f"Skipping group {sample} with {len(indices)} cells, which is less than {n_cells * n_pseudo_replicates}")
            continue
        
        np.random.seed(1)
        np.random.shuffle(indices)
        indices = np.array_split(np.array(indices), n_pseudo_replicates) #change number here for number of replicates deisred
        
        for i, pseudo_rep in enumerate(indices):
            rep_adata = sc.AnnData(
                X = samp_cell_subset[indices[i]].X.sum(axis = 0),
                var = samp_cell_subset[indices[i]].var)

            rep_adata.obs_names = [sample + '_' + str(i)]
            rep_adata.obs['condition'] = samp_cell_subset.obs[condition_key].iloc[0]
            rep_adata.obs['replicate'] = i
            rep_adata.obs['psbulk_n_cells'] = len(indices[i])
            rep_adata.obs['psbulk_counts'] = rep_adata.X.sum()

            pbs.append(rep_adata)
        
    pb = sc.concat(pbs)
    # sc.pp.filter_genes(pb, min_cells = 1)

    return pb


def run_deseq2(
        pb: sc.AnnData,
        condition_key: str = 'condition',
        treatment_key: str = 'treat',
        control_key: str = 'ctrl',
        ref_level: list = None,
        plot_pca: bool = False,
        lfc_shrink: bool = True,
        n_cpus: int = 8):
    """
    Run DESeq2 on a pseudo-bulk sample.

    Parameters
    ----------
    pb
        An AnnData object containing a pseudo-bulk sample.
    condition_key
        The key in `pb.obs` containing the condition column.
    treatment_key
        The key in `pb.obs.condition_key` containing the treatment name.
    control_key
        The key in `pb.obs.condition_key` containing the control name.
    plot_pca
        Whether to plot PCA.
    lfc_shrink
        Whether to shrink log fold changes.
    n_cpus
        The number of CPUs to use.
    
    Returns
    -------
    de
        A dataframe containing the DESeq2 results.
    dds
        The DESeq2 object.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    
    
    counts = pd.DataFrame(pb.X, columns = pb.var_names, index=pb.obs_names) # need to do this to pass var names
    inference = DefaultInference(n_cpus=n_cpus)

    dds = DeseqDataSet(
        counts = counts,
        metadata=pb.obs,
        design_factors=condition_key,
        ref_level=ref_level,
        refit_cooks=True,
        inference=inference)

    dds.deseq2()

    if plot_pca:
        sc.tl.pca(dds)
        sc.pl.pca(dds, color = condition_key, size = 200)
    
    stat_res = DeseqStats(dds, contrast=(condition_key, treatment_key, control_key))
    stat_res.summary()

    if lfc_shrink:
        _coeff = f'{condition_key}_{treatment_key}_vs_{control_key}'
        stat_res.lfc_shrink(coeff=_coeff)
    
    de = stat_res.results_df
    de = de.sort_values('stat', ascending = False)

    return de, dds
