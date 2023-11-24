from math import ceil
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import scipy.sparse as ss
from ..utils.rtools import r2py, py2r, r_inline_plot
import pandas as pd


def AUCell_r(
    adata,
    genesets: dict,
    layer: str = 'counts',
    n_jobs: int = 5,
    label="AUCell",
    n_cols=3,
    show=False,
):
    """
    Run AUCell analysis on a single AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    genesets : dict
        Dictionary of gene sets.
    layer : str, optional
        Layer to use for analysis, by default 'counts'.
    n_jobs : int, optional
        Number of cores to use for parallel processing, by default 5.
    label : str, optional
        Label for the output, by default "AUCell".
    n_cols : int, optional
        Number of columns for the plot, only works when show is True, by default 3.
    show : bool, optional
        Whether to show the plot, by default False.
    """

    rBase = importr("base")
    dplyr = importr("dplyr")
    aucell = importr("AUCell")
    BiocParallel = importr("BiocParallel")

    R = ro.r
    rEnv = ro.globalenv

    def get_threshold(objR_name, geneCate):
        df = r2py(R(f"as.data.frame({objR_name}$`{geneCate}`$aucThr$thresholds)"))
        df = df.assign(geneCate=geneCate)

        return df

    
    def get_assignment(objR_name, geneCate):
        # print(f"c({objR_name}$`{geneCate}`$assignment)")
        cell_name = list(r2py(R(f"{objR_name}$`{geneCate}`$assignment")))

        return cell_name

    n_rows = ceil(len(genesets) / n_cols)
    var_keys = adata.var.index.to_list()  # gene names
    obs_keys = adata.obs.index.to_list()  # cell names
    mtx = adata.layers[layer].T

    # if ss.issparse(mtx) & forceDense:
    #     mtx = mtx.A
    mtxR = py2r(mtx)
    del mtx

    fc_list2R = lambda x: R.unlist(R.list(x))
    obs_keys = fc_list2R(obs_keys)
    var_keys = fc_list2R(var_keys)
    genesets_r= R.list(**{x: fc_list2R(y) for x, y in genesets.items()})

    rEnv["mtxR"] = mtxR
    rEnv["obs_keys"] = obs_keys
    rEnv["var_keys"] = var_keys
    rEnv["genesets"] = genesets_r
    rEnv["n_jobs"] = n_jobs
    # rEnv["aucMaxRank"] = aucMaxRank
    rEnv["n_cols"] = n_cols
    rEnv["n_rows"] = n_rows

    R(
        """
        set.seed(0)
        BPPARAM=BiocParallel::MulticoreParam(n_jobs)
        rownames(mtxR) <- var_keys
        colnames(mtxR) <- obs_keys
        # cells_rankings <- AUCell_buildRankings(mtxR, nCores=n_jobs, plotStats=TRUE)
        # cells_AUC <- AUCell_calcAUC(dtR_genes, cells_rankings, aucMaxRank=aucMaxRank)
        cells_AUC <- AUCell_run(mtxR, genesets)
        """)
    
    if show:
        with r_inline_plot(width=600):
            R(
                """
                par(mfrow=c(n_rows, n_cols))
                cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=T, assign=T)
                """)
    else:
        R('cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F)')
    df_auc = r2py(R("as.data.frame(cells_AUC@assays@data$AUC)")).T
    df_aucThreshold = pd.concat(
        [get_threshold("cells_assignment", x) for x in genesets.keys()]
    )
    adata.obsm[label] = df_auc.copy()
    adata.uns[f'{label}_threshold'] = df_aucThreshold.copy()
    adata.uns[f'{label}_assignment'] = {x: get_assignment("cells_assignment", x) for x in genesets.keys()}