from collections import namedtuple

CC_genes = namedtuple('CC_genes', ['s_genes', 'g2m_genes'])
cc_genes_human = CC_genes(
    ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG',
    'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP',
    'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76',
    'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2',
    'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2',
    'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8'],
    
    ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A',
    'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF',
    'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB',
    'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP',
    'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1',
    'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR',
    'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2',
    'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
)

# convert by biomaRt
cc_genes_mouse = CC_genes(
    ['Mcm4', 'Rrm2', 'Chaf1b', 'Gmnn', 'Uhrf1', 'Mcm2', 'Msh2', 'Cdc45',
    'Exo1', 'Slbp', 'Rad51ap1', 'Ubr7', 'Hells', 'Cdc6', 'Fen1',
    'Dscc1', 'Nasp', 'Tyms', 'Wdr76', 'Gins2', 'Dtl', 'Mcm6', 'Tipin',
    'Ung', 'Usp1', 'Mcm5', 'Pola1', 'Blm', 'Clspn', 'Rpa2', 'Cdca7',
    'Prim1', 'Ccne2', 'Casp8ap2', 'Brip1', 'Rad51', 'Rrm1', 'E2f8',
    'Pcna', 'Rfc2'],

    ['Cbx5', 'Gtse1', 'Dlgap5', 'Cdc25c', 'Smc4', 'Ctcf', 'Cdca2',
    'Kif20b', 'Cdca3', 'Anp32e', 'Cks2', 'Ttk', 'Gas2l3', 'Top2a',
    'Ncapd2', 'Tacc3', 'Ckap2', 'Hmgb2', 'Cdk1', 'G2e3', 'Tpx2',
    'Aurkb', 'Lbr', 'Kif11', 'Cenpe', 'Kif2c', 'Kif23', 'Cenpa',
    'Rangap1', 'Hjurp', 'Cenpf', 'Cks1b', 'Ect2', 'Birc5', 'Ccnb2',
    'Cdc20', 'Cdca8', 'Hmmr', 'Ndc80', 'Anln', 'Psrc1', 'Nek2',
    'Nusap1', 'Ube2c', 'Mki67', 'Bub1', 'Aurka', 'Ckap2l', 'Nuf2',
    'Ckap5', 'Tubb4b']
)


def load_cell_cycle_genes(species: str = 'human') -> tuple:
    """Load cell cycle genes for human or mouse.
    Args:
        species: 'human' or 'mouse'.
    Returns:
        A tuple of two lists: (S phase genes, G2M phase genes).
    """
    if species == 'human':
        return cc_genes_human.s_genes, cc_genes_human.g2m_genes
    elif species == 'mouse':
        return cc_genes_mouse.s_genes, cc_genes_mouse.g2m_genes
    else:
        raise ValueError(f'Unknown species: {species}')