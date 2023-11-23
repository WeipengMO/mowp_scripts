import scanpy as sc
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from ..utils import rtools
from loguru import logger


def remove_ambient(
    adata: sc.AnnData, 
    adata_raw: sc.AnnData,
    soupX_ambient_frac: bool = True,
    ):
    """Remove ambient RNA from the data by `SoupX`

    Parameters
    ----------
        adata: Annotated data of filtered_feature_bc_matrix.
        adata_raw: Annotated data of raw_feature_bc_matrix.
        soupX_ambient_frac: Whether to add the ambient fraction to the obs.
    """

    SoupX = importr('SoupX')
    assert (adata.var.index == adata_raw.var.index).all(), 'The order of genes in the raw and filtered matrices must be the same'

    logger.info("Preprocessing data for SoupX")
    # normalize and log1p transform it
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)

    # run PCA and clustering
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")

    # Preprocess variables for SoupX
    soupx_groups = adata_pp.obs["soupx_groups"]
    del adata_pp

    # Prepare the data for SoupX
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T

    # raw gene by cells matrix
    data_tod = adata_raw.X.T
    del adata_raw

    # Run SoupX
    logger.info("Running SoupX")
    with ro.local_context() as rEnv:
        rEnv['data'] = rtools.py2r(data)
        rEnv['data_tod'] = rtools.py2r(data_tod)
        rEnv['cells'] = rtools.py2r(cells)
        rEnv['genes'] = rtools.py2r(genes)
        rEnv['soupx_groups'] = ro.r.c(**soupx_groups)
        
        ro.r(
            '''
            set.seed(0)
            # specify row and column names of data
            rownames(data) = genes
            colnames(data) = cells
            # ensure correct sparse format for table of counts and table of droplets
            data <- as(data, "sparseMatrix")
            data_tod <- as(data_tod, "sparseMatrix")

            # Generate SoupChannel Object for SoupX 
            sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

            # Add extra meta data to the SoupChannel object
            soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
            sc = setSoupProfile(sc, soupProf)
            # Set cluster information in SoupChannel
            sc = setClusters(sc, soupx_groups)

            # Estimate contamination fraction
            sc  = autoEstCont(sc, doPlot=FALSE)
            # Infer corrected table of counts and rount to integer
            out = adjustCounts(sc, roundToInt = TRUE)
        ''')

        adata.layers["soupX_counts"] = rtools.r2py(rEnv['out']).T
        adata.X = adata.layers["soupX_counts"]

        if soupX_ambient_frac:
            adata.obs['soupX_ambient_frac'] = 1 - (adata.layers['soupX_counts'].sum(1).A.reshape(-1) / adata.layers['counts'].sum(1).A.reshape(-1))