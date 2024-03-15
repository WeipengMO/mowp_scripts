import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from loguru import logger


def read_st_pipeline(stdata: str, library_id: str, img_path: str):
    '''
    Read spatial transcriptomics data from st_pipeline. 
    https://github.com/di-0579/Spatial_epigenome-transcriptome_co-sequencing

    Parameters
    ----------
    stdata : `str`
        Path to the stdata file.
    library_id : `str`
        Library id.
    img_path : `str`
        Path to the image folder.
    '''
    from squidpy._constants._pkg_constants import Key
    from squidpy.read._utils import _load_image
    import json
    from pathlib import Path


    df = pd.read_csv(stdata, sep='\t', index_col=0)
    adata = ad.AnnData(
        X=csc_matrix(df.values.astype('float32')), 
        obs=pd.DataFrame(index=df.index),
        var=pd.DataFrame(index=df.columns))
    del df

    path = Path(img_path)

    tissue_positions_file = (
        path / "tissue_positions.csv"
        if (path / "tissue_positions.csv").exists()
        else path / "tissue_positions_list.csv"
    )
    coords = pd.read_csv(
        tissue_positions_file,
        header=1 if tissue_positions_file.name == "tissue_positions.csv" else None,
    )
    coords.columns = ['barcode', "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
    coords.index = coords.apply(lambda x: f'{x.array_row}x{x.array_col}', axis=1)

    if len(set(adata.obs_names) - set(coords.index)) > 0:
        logger.warning('Some barcodes in adata are not in coords, please check')
        return coords

    coords = coords.loc[adata.obs_names]
    adata.obs = coords

    adata.uns[Key.uns.spatial] = {library_id: {"metadata": {}}}
    adata.uns[Key.uns.spatial][library_id][Key.uns.image_key] = {
        res: _load_image(path/ f"tissue_{res}_image.png") for res in ["hires", "lowres"]
    }
    adata.uns[Key.uns.spatial][library_id]["scalefactors"] = json.loads(
        (path/ "scalefactors_json.json").read_bytes()
    )
    adata.obsm[Key.obsm.spatial] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    adata = adata[adata.obs['in_tissue'].astype(bool)].copy()

    return adata