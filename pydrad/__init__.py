"""
pydrad--a Python package for configuring, parsing and visualizing
HYDRAD simulations.
"""
import logging

from astropy.logger import AstropyLogger

from pydrad.version import version as __version__


def _init_log():
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

__all__ = ['__version__', 'log']
