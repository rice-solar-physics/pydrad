"""
pydrad--a Python package for configuring, parsing and visualizing
HYDRAD simulations.
"""
import logging
from astropy.logger import AstropyLogger


def _init_log():
    """
    Initializes the SunPy log.
    In most circumstances this is called automatically when importing
    SunPy. This code is based on that provided by Astropy see
    "licenses/ASTROPY.rst".
    """
    orig_logger_cls = logging.getLoggerClass()
    logging.setLoggerClass(AstropyLogger)
    try:
        log = logging.getLogger('pydrad')
        log._set_defaults()
    finally:
        logging.setLoggerClass(orig_logger_cls)
    return log


log = logging.getLogger()
log = _init_log()
