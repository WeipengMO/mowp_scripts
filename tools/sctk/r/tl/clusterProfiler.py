import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from ...utils import rtools


def go_analysis(
        gene_list: list, 
        species: str = 'human',
        ont: str = None,
        keytype: str = 'SYMBOL',
        ) -> dict:
    """
    Perform GO analysis using clusterProfiler.

    Parameters
    ----------
    gene_list
        A list of gene symbols.
    species
        The species of the gene list. One of 'human', 'mouse'.
    ont
        The ontology to use. One of 'MF', 'BP', 'CC'. If None, all ontologies will be used.
    keytype
        The keytype of the gene list.
    
    Returns
    -------
    go_enrich
        A dict of enriched GO terms.
    
    Examples
    --------
    >>> gene_list = ['TP53', 'BRCA1', 'BRCA2', 'CDKN2A', 'CDKN2B']
    >>> go_enrich = go_analysis(gene_list, species='human')
    """

    if species == 'human':
        orgDb = 'org.Hs.eg.db'
    elif species == 'mouse':
        orgDb = 'org.Mm.eg.db'
    else:
        raise ValueError(f'Unsupported species: {species}')
    
    with ro.local_context() as rEnv:
        clusterProfiler = importr('clusterProfiler')
        rEnv['gene_list'] = ro.r.c(*gene_list)
        rEnv['orgDb'] = rtools.py2r(orgDb)
        rEnv['keytype'] = rtools.py2r(keytype)
        if ont is None:
            ont = ['MF', 'BP', 'CC']
        elif isinstance(ont, str):
            ont = [ont]
        else:
            raise ValueError(f'Unsupported ont: {ont}')

        go_enrich = {}
        for _ont in ont:
            rEnv['ont'] = rtools.py2r(_ont)
            ro.r(
                f"""
                go_enrich <- enrichGO(gene          = gene_list,
                                    OrgDb         = orgDb,
                                    ont           = ont,
                                    keyType       = keytype)
                go_enrich <- as.data.frame(go_enrich)
                """
            )
            _go_enrich = rtools.r2py(rEnv['go_enrich'])
            go_enrich[_ont] = _go_enrich
    
    return go_enrich