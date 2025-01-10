from loguru import logger
import sys


def set_seed(seed=0):
    import os
    import numpy as np
    import random
    
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
    
    try:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.enabled = False
    except ModuleNotFoundError:
        logger.warning("torch not found, skipping seed setting for torch")


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
