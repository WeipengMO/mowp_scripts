import requests
from loguru import logger
import json


def convert_gene_symbol(genes: str, species='human'):
    """
    Convert gene symbol to specified species

    Parameters
    ----------
    genes : str
        Gene symbol to convert
    species : str
        Species to convert to
    
    Returns
    -------
    convert_genes : list
        List of gene symbols
    """
    def run_mygene(genes, species='human'):
        mg = mygene.MyGeneInfo()
        if isinstance(genes, str):
            genes = [genes]

        query = mg.querymany(genes, scopes='symbol', species=species)
        
        convert_genes = []
        for entry in query:
            if 'symbol' in entry:
                convert_genes.append(entry['symbol'])

        return convert_genes

    def run_mygene_request(gene_symbol, species='human'):
        api_url = "http://mygene.info/v3/query"
        params = {
            "q": f"symbol:{gene_symbol}",
            "species": species,
            "size": 1
        }
        response = requests.get(api_url, params=params)
        
        if response.status_code == 200:
            data = response.json()
            if data["hits"]:
                gene_symbol = data['hits'][0]['symbol']
                return gene_symbol
            else:
                return None
        else:
            raise f"Error: Unable to reach mygene.info API. Status code: {response.status_code}"
    
    if isinstance(genes, str):
        genes = [genes]
    
    try:
        import mygene
        return run_mygene(genes, species)
    except ImportError:
        logger.info("mygene package not found. use request API instead")
        convert_genes = []
        no_hit_genes = []
        for gene in genes:
            _convert_gene = run_mygene_request(gene, species)
            if _convert_gene is not None:
                convert_genes.append(_convert_gene)
            else:
                no_hit_genes.append(gene)
        
        logger.info(f"Genes not found: {','.join(no_hit_genes)}")

        return convert_genes
    

def ensemble_to_symbol(ensemble_id: list):
    """
    Convert ensemble id to gene symbol

    Parameters
    ----------
    ensemble_id : list
        List of ensemble id
    
    Returns
    -------
    output : dict
        Dictionary of ensemble id and gene symbol
    """
    if not isinstance(ensemble_id, list):
        ensemble_id = list(ensemble_id)
        
    url = "https://biotools.fr/human/ensembl_symbol_converter/"
    body = {'api': 1, 'ids': json.dumps(ensemble_id)}
    r = requests.post(url, data = body)
    output = r.json()

    return output