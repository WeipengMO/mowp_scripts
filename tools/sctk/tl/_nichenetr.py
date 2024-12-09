from ..utils import rtools
from ..utils.rtools import rcontext
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import scanpy as sc

R = ro.r


@rcontext
def _load_prior_knowledge(organism, rEnv=None):
    import yaml
    with open('tl/nichenet_config.ymal') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    config = config[organism]
    lr_network, ligand_target_matrix, weighted_networks = (
        config['lr_network'], 
        config['ligand_target_matrix'], 
        config['weighted_networks'])

    R(f'''
        lr_network <- readRDS('{lr_network}')
        ligand_target_matrix <- readRDS('{ligand_target_matrix}')
        weighted_networks <- readRDS('{weighted_networks}')
    ''')



@rcontext
def run_nichenetr(
        adata: sc.AnnData,
        organism: str,
        sender_celltypes: list,
        receiver_celltypes: str,
        condition_colname: str,
        condition_oi: str,
        condition_reference: str,
        expression_pct = 0.1,
        cell_type_key: str = None,
        geneset: str = 'DE',
        lfc_cutoff: float = 0.25,
        rEnv=None,
        ):
    """Run Nichenet analysis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cell_type_key : str
        Key in adata.obs that contains cell type information.
    organism : str
        Organism of the data. Must be either 'human' or 'mouse'.
    sender_celltypes : list
        List of sender cell types.
    receiver_celltypes : str
        Receiver cell type.
    condition_colname : str
        Column in adata.obs that contains condition information.
    condition_oi : str
        Condition of interest.
    condition_reference : str
        Reference condition.
    geneset : str
        Geneset to use. Must be either 'DE', 'up' or 'down'.
    lfc_cutoff : float
        Log fold change cutoff.
    rEnv : dict, optional
        Environment to run the R code in.
    """
    importr('Seurat')
    importr('nichenetr')
    importr('tidyverse')

    if isinstance(adata, sc.AnnData):
        if cell_type_key is None:
            raise ValueError('cell_type_key must be specified.')
        
        seuratObj = rtools.ad2so(adata, layer='counts')
        rEnv['seuratObj'] = seuratObj
        R(f"Idents(object = seuratObj) <- '{cell_type_key}'")
    elif isinstance(adata, ro.methods.RS4) and rtools.r2py(R('class')(adata))[0] == 'Seurat':
        seuratObj = adata
        rEnv['seuratObj'] = seuratObj
    else:
        raise ValueError('adata must be either a Seurat object or an AnnData object.')
    

    rEnv['sender_celltypes'] = rtools.py2r(sender_celltypes)
    rEnv['receiver_celltypes'] = rtools.py2r(receiver_celltypes)
    rEnv['condition_colname'] = rtools.py2r(condition_colname)
    rEnv['condition_oi'] = rtools.py2r(condition_oi)
    rEnv['condition_reference'] = rtools.py2r(condition_reference)
    rEnv['expression_pct'] = rtools.py2r(expression_pct)
    rEnv['geneset'] = rtools.py2r(geneset)
    rEnv['lfc_cutoff'] = rtools.py2r(lfc_cutoff)

    if organism not in ['human', 'mouse']:
        raise ValueError('organism must be either "human" or "mouse"')
    _load_prior_knowledge(organism, rEnv=rEnv)

    R('''
        nichenet_output <- nichenet_seuratobj_aggregate(
            seurat_obj = seuratObj, 
            sender = sender_celltypes, 
            receiver = receiver_celltypes, 
            condition_colname = condition_colname,
            condition_oi = condition_oi,
            condition_reference = condition_reference,
            expression_pct = expression_pct,
            ligand_target_matrix = ligand_target_matrix,
            lr_network = lr_network,
            weighted_networks = weighted_networks,
            geneset = geneset,
            lfc_cutoff = lfc_cutoff,
        )
    ''')

    return rEnv['nichenet_output']
