import anndata as ad
import os
from functools import wraps
import scanpy as sc
from loguru import logger
from .common import datasetdir


def check_datasetdir_exists(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        if not os.path.isdir(datasetdir):
            logger.info(f"Creating '{os.path.realpath(datasetdir)}' directory")
            os.mkdir(datasetdir)
        return f(*args, **kwargs)

    return wrapper


@check_datasetdir_exists
def NeurIPS_filtered() -> ad.AnnData:
    """\
    10x Multiome data set generated for a single cell data integration challenge at 
    the NeurIPS conference 2021. This dataset captures single-cell multiomics data 
    from bone marrow mononuclear cells of 12 healthy human donors measured at 
    four different sites to obtain nested batch effects. 

    """
    filename = f'{datasetdir}/NeurIPS_filtered_feature_bc_matrix.h5ad'
    url = "https://figshare.com/ndownloader/files/39546196"

    adata = sc.read_10x_h5(
        filename=filename,
        backup_url=url,
    )

    return adata


@check_datasetdir_exists
def NeurIPS_raw() -> ad.AnnData:
    """\
    raw_feature_bc_matrix for NeurIPS datasets

    """
    filename = f'{datasetdir}/NeurIPS_raw_feature_bc_matrix.h5ad'
    url = "https://figshare.com/ndownloader/files/39546217"

    adata = sc.read_10x_h5(
        filename=filename,
        backup_url=url,
    )

    return adata


@check_datasetdir_exists
def s4d8_clustered() -> ad.AnnData:
    """\
    Sample site4-donor8 from the NeurIPS human bone marrow dataset.

    """
    filename = f'{datasetdir}/s4d8_clustered.h5ad'
    url = "https://figshare.com/ndownloader/files/41436666"

    adata = sc.read(
        filename=filename,
        backup_url=url,
    )

    return adata


@check_datasetdir_exists
def lung_atlas() -> ad.AnnData:
    """\
    lung atlas integration task from the scIB manuscript.
    https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/harmonization.html

    """
    filename = f'{datasetdir}/lung_atlas.h5ad'
    url = "https://figshare.com/ndownloader/files/24539942"

    adata = sc.read(
        filename=filename,
        backup_url=url,
    )

    return adata