import numpy as np
import pandas as pd
import scanpy as sc
from typing import Union, List
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from loguru import logger
from ..utils import rtools
from ..utils.rtools import rcontext
import rpy2.robjects as ro
from rpy2.robjects.packages import importr


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

    pbs = []
    for sample in cell_subset.obs[sample_key].unique():
        samp_cell_subset = cell_subset[cell_subset.obs[sample_key] == sample]
        
        samp_cell_subset.X = samp_cell_subset.layers[counts_key] #make sure to use raw data
        
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
    
    counts = pd.DataFrame(pb.X, columns = pb.var_names, index=pb.obs_names) # need to do this to pass var names
    inference = DefaultInference(n_cpus=n_cpus)

    dds = DeseqDataSet(
        counts = counts,
        metadata=pb.obs,
        design_factors=condition_key,
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


@rcontext
def run_edgeR(
    pb: sc.AnnData, 
    condition_key: str = 'condition', 
    treatment_key: str = 'treat', 
    control_key: str = 'ctrl',
    rEnv=None) -> pd.DataFrame: 
    """
    Run edgeR on a pseudo-bulk sample.

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
    rEnv
        rpy2.robjects.r environment to use. If None, a new one will be created.

    Returns
    -------
    res
        A dataframe containing the edgeR results.
    
    """

    importr('edgeR')
    R = ro.r
    R('''
        fit_model <- function(adata_, label){
            # create an edgeR object with counts and grouping factor
            y <- DGEList(assay(adata_, "X"), group = adata_[[label]])
            # filter out genes with low counts
            print("Dimensions before subsetting:")
            print(dim(y))
            print("")
            keep <- filterByExpr(y)
            y <- y[keep, , keep.lib.sizes=FALSE]
            print("Dimensions after subsetting:")
            print(dim(y))
            print("")
            # normalize
            y <- calcNormFactors(y)
            # create a vector that is concatentation of condition and cell type that we will later use with contrasts
            group <- adata_[[label]]
            design <- model.matrix(~ 0 + group)
            # estimate dispersion
            y <- estimateDisp(y, design = design)
            # fit the model
            fit <- glmQLFit(y, design)
            
            return(list("fit"=fit, "design"=design, "y"=y))}
    ''')

    obs_keep = [condition_key]
    pb.obs = pb.obs[obs_keep]

    adata_sce = rtools.py2r(pb)
    rEnv['outs'] = rEnv['fit_model'](adata_sce, condition_key)

    contrasts = f'group{treatment_key}-group{control_key}'
    R(f'''
        fit <- outs$fit
        y <- outs$y
        myContrast <- makeContrasts('{contrasts}', levels = y$design)
        qlf <- glmQLFTest(fit, contrast=myContrast)
        # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
        tt <- topTags(qlf, n = Inf)
        tt <- tt$table
    ''')

    res = rtools.r2py(rEnv['tt'])

    return res
    