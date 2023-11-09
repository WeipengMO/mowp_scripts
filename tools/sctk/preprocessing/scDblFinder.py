import scanpy as sc
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from ..utils import rtools
from loguru import logger


def scDblFinder(
    adata: sc.AnnData,
    batch_key: str = None,
    dbr: int = None,
    n_jobs: int = 4,
    drop_doublet: bool = True,
    ):
    """
    Detects doublets in single-cell RNA-seq data using the scDblFinder R package.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str, optional
        Key for batch information in adata.obs.
    dbr : int, optional
        The expected number of doublets in the dataset.
    n_jobs : int, optional
        Number of cores to use for parallel processing.
    droup_doublet : bool, optional
        Drop doublets from the AnnData object.

    Examples
    --------
    >>> sctk.pp.scDblFinder(adata)
    
    """
def scDblFinder(
    adata: sc.AnnData,
    batch_key: str = None,
    dbr: int = None,
    n_jobs: int = 4,
    drop_doublet: bool = True,
):

    scDblFinder = importr('scDblFinder')
    BiocParallel = importr('BiocParallel')

    with ro.local_context() as rEnv:
        if n_jobs:
            rEnv['n_jobs'] = n_jobs

        if batch_key is not None:
            rEnv['batch_key'] = batch_key
        else:
            rEnv['batch_key'] = ro.r('NULL')
        
        if dbr is not None:
            rEnv['dbr'] = dbr
        else:
            rEnv['dbr'] = ro.r('NULL')

        rEnv['data_mat'] = rtools.py2r(adata.X.T)
        rEnv['colData'] = rtools.py2r(adata.obs)


        ro.r(
            '''
            set.seed(1)

            sce <- SingleCellExperiment(assays = list(counts = data_mat), colData = colData)

            sce <- scDblFinder(
                sce,
                samples=batch_key,
                dbr=dbr,
                BPPARAM=MulticoreParam(n_jobs),
                )

            doublet_score = sce$scDblFinder.score
            doublet_class = sce$scDblFinder.class
            ''')
        
        doublet_score = rtools.r2py(rEnv['doublet_score'])
        doublet_class = rtools.r2py(rEnv['doublet_class'])
            
        adata.obs["scDblFinder_score"] = doublet_score
        adata.obs["scDblFinder_class"] = doublet_class

        if drop_doublet:
            logger.info(f"Before dropping: {len(adata)}")
            adata = adata[adata.obs['scDblFinder_class'] == 'singlet'].copy()
            logger.info(f"After dropping: {len(adata)}")
        else:
            doublet_info = adata.obs.scDblFinder_class.value_counts().to_dict()
            logger.info(', '.join(f'{k}: {v}' for k, v in doublet_info.items()))
            logger.info(f"Anndata still contains doublets")
        
        return adata
