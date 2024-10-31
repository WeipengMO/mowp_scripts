from loguru import logger

from . import rtools
from ._genes import (
    convert_gene_symbol,
    ensemble_to_symbol
)
from .settings import (
    configure_logger,
    set_seed
)
try:
    from ._survival import Survival
except ImportError:
    logger.warning("Could not import Survival class. Please install lifelines package: pip install lifelines")
    pass
