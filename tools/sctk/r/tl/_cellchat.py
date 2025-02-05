from ...utils import rtools
from ...utils.rtools import rcontext
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
R = ro.r

@rcontext
def adata_to_cellchat(adata, cell_type, rEnv=None):
    importr('CellChat')
    importr('patchwork')
    mat = rtools.py2r(adata.layers['log1p_norm'].T)
    meta = rtools.py2r(adata.obs)
    cell_ident = rtools.py2r(adata.var_names)

    rEnv['mat'] = mat
    rEnv['meta'] = meta
    rEnv['cell_ident'] = cell_ident
    rEnv['cell_type'] = cell_type

    R('''
    colnames(mat) <- rownames(meta)
    rownames(mat) <- cell_ident
    cellchat <- createCellChat(object = mat, meta = meta, group.by = cell_type)
      ''')
    
    return rEnv['cellchat']