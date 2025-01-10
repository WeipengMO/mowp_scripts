from ...utils import rtools
from ...utils.rtools import rcontext
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import scanpy as sc
from loguru import logger
import os
import numpy as np


@rcontext
def run_inferCNV_r(
        adata_merge: sc.AnnData,
        gene_order_file: str,
        ref_group_names: list = ['neutrophils', 'macrophages'],
        layer: str = 'counts', 
        cell_type: str = 'cell_type',
        cutoff: float = 0.1,
        out_dir: str = 'inferCNV_out',
        cluster_by_groups: bool = True,
        HMM: bool = False,
        denoise: bool = True,
        num_threads: int = 8,
        rEnv=None,
        ):
    """Run inferCNV.

    Parameters
    ----------
    adata_merge
        AnnData object with merged samples.
    gene_order_file
        Path to gene order file.
    ref_group_names
        A list containing the classifications of the reference (normal) cells to use for infering cnv
    layer
        Adata Layer to use.
    cell_type
        Adata.obs Column name of cell type annotation.
    cutoff
        cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir
        path to directory to deposit outputs
    cluster_by_groups
        If observations are defined according to groups (ie. patients), each group of cells will be clustered separately. (default=False, instead will use k_obs_groups setting)
    HMM
        when set to True, runs HMM to predict CNV level (default: False)
    denoise
        when set to True, runs denoising algorithm (default: True)
    num_threads
        number of threads to use (default: 24)
    rEnv
        rpy2.robjects.r environment to use. If None, a new one will be created.
    """

    if out_dir[0] != '/':
        cwd = os.getcwd()
        out_dir = os.path.join(cwd, out_dir)

    logger.info(f'output directory: {out_dir}')

    if isinstance(adata_merge.layers[layer], np.ndarray):
        raw_counts_matrix = adata_merge.layers[layer].T
    else:
        raw_counts_matrix = adata_merge.layers[layer].A.T
    
    annotations_file = adata_merge.obs[[cell_type]]

    rEnv['raw_counts_matrix'] = rtools.py2r(raw_counts_matrix)
    del raw_counts_matrix
    rEnv['var_names'] = rtools.py2r(adata_merge.var_names)
    rEnv['obs_names'] = rtools.py2r(adata_merge.obs_names)
    ro.r('rownames(raw_counts_matrix) <- var_names')
    ro.r('colnames(raw_counts_matrix) <- obs_names')

    rEnv['annotations_file'] = rtools.py2r(annotations_file)
    del annotations_file
    rEnv['gene_order_file'] = gene_order_file
    rEnv['ref_group_names'] = ro.r.c(*ref_group_names)

    inferCNV_r = importr('infercnv')
    R = ro.r

    logger.info('Running inferCNV CreateInfercnvObject...')
    R(
        '''
        infercnv_obj <- infercnv::CreateInfercnvObject(
            raw_counts_matrix=raw_counts_matrix,
            annotations_file=annotations_file,
            gene_order_file=gene_order_file,
            delim="\t",
            ref_group_names=ref_group_names)
        '''
    )

    logger.info('Running inferCNV run...')
    rEnv['cluster_by_groups'] = cluster_by_groups
    rEnv['HMM'] = HMM
    rEnv['denoise'] = denoise
    rEnv['cutoff'] = cutoff
    rEnv['out_dir'] = out_dir
    rEnv['num_threads'] = num_threads
    R(
        f'''
        infercnv_obj <- infercnv::run(
                infercnv_obj,
                cutoff=cutoff, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                out_dir=out_dir,
                cluster_by_groups=cluster_by_groups, 
                plot_steps=F,
                denoise=denoise,
                HMM=HMM,
                no_prelim_plot=T,
                write_expr_matrix=T,
                num_threads = num_threads)
        '''
    )