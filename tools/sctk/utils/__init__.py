from . import rtools
import requests
import json


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
    url = "https://biotools.fr/human/ensembl_symbol_converter/"
    body = {'api': 1, 'ids': json.dumps(ensemble_id)}
    r = requests.post(url, data = body)
    output = r.json()

    return output