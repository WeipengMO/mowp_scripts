from ...utils import rtools
from ...utils.rtools import rcontext
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import scanpy as sc
import pandas as pd


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
    