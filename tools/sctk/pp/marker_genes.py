from typing import Union, List, Dict
import scanpy as sc

def find_marker_genes_in_data(adata, marker_genes: Union[List[str], Dict[str, List[str]]]):
    """
    Find marker genes in the adata.var.index.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    marker_genes : list or dict
        List of marker genes or dictionary of cell types and marker genes.
    """
    if isinstance(marker_genes, dict):
        marker_genes_in_data = dict()
        for ct, markers in marker_genes.items():
            markers_found = list()
            for marker in markers:
                if marker in adata.var.index:
                    markers_found.append(marker)
            if len(markers_found) > 0:
                marker_genes_in_data[ct] = markers_found
    
    elif isinstance(marker_genes, list):
        marker_genes_in_data = list()
        for marker in marker_genes:
            if marker in adata.var.index:
                marker_genes_in_data.append(marker)
    
    else:
        raise ValueError("marker_genes must be a list or a dictionary.")
    
    return marker_genes_in_data


def get_rank_genes_from_groups(
    adata, 
    groupby: str = 'leiden', 
    key: str = 'dea_leiden_filtered',
    n_genes=None,
    print_rank_genes=True,
    return_rank_genes=False):
    '''Get rank genes from `sc.tl.filter_rank_genes_groups`

    Parameters
    ----------
    adata: AnnData
    groupby: str
        The key of the observation grouping to consider.
    key: str
        The key of the ranking to consider.
    n_genes: int
        Number of genes to show.
    print_rank_genes: bool
        Whether print rank genes.
    return_rank_genes: bool
        Whether return rank genes.
    
    Examples
    --------
    >>> get_rank_genes_from_groups(adata, groupby='leiden', key='dea_leiden_filtered', n_genes=5)
    >>> 1: CXADR, EGFR, ESRP1, ITGA2, BAIAP2L1
        2: AIF1, CYBB, ITGB2, FAM49A, CD86
        3: SKAP1, CD3D, TRAC, CD247, TRBC2
        4: PRRX1, FBN1, C1R, FAP, BICC1
        5: MZB1, IGHGP, DERL3, FKBP11, IGHG2
        6: PECAM1, CALCRL, ADGRL4, LDB2, SPARCL1
        7: BCL2A1, AQP9, CSF3R, LUCAT1, PLEK
        8: TPSB2, CPA3, TPSAB1, HPGDS, LTC4S

    '''
    rank_genes = {}
    for i in adata.obs[groupby].cat.categories:
        _data = sc.get.rank_genes_groups_df(adata, group=i, key=key).dropna()['names'].values
        if len(_data) > 0:
            if n_genes is not None:
                rank_genes[i] = _data[: n_genes]
                if print_rank_genes:
                    print(f"{i}: {', '.join(_data[: n_genes])}")
            else:
                rank_genes[i] = _data
                if print_rank_genes:
                    print(f"{i}: {', '.join(_data)}")
    
    if return_rank_genes:
        return rank_genes