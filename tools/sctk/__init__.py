#!/usr/bin/env python3
# coding: utf-8

from loguru import logger
from . import settings
settings.set_seed()

from .common import version as __version__
from . import pp, tl, pl, read, sample_data

try:
    import rpy2.robjects as ro
    R = ro.r
    seed = 0
    R("set.seed")(seed)

    from . import r, utils

except ModuleNotFoundError:
    logger.warning("R module could not be imported. Please install rpy2 (pip install rpy2).")
