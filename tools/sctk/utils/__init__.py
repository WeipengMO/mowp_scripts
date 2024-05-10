from . import rtools
import requests
import json
from loguru import logger
import sys


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


def configure_logger(log_level: str):
    """
    Configures the logger with the specified log level.

    Parameters
    ----------
    log_level (str): The desired log level ('info', 'debug', 'warning', etc.).
    """
    valid_log_levels = {"trace", "debug", "info", "success", "warning", "error", "critical"}
    if log_level.lower() not in valid_log_levels:
        raise ValueError(f"Invalid log level '{log_level}'. Valid levels are: {', '.join(valid_log_levels)}")

    logger.remove()  # Remove any existing handlers
    logger.add(sys.stdout, level=log_level.upper())  # Add a new handler for console output