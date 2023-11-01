"""
Custom Jinja2 filters
"""
import numpy as np
import astropy.units as u
import plasmapy.particles
from jinja2 import Undefined

from .util import MissingParameter

__all__ = [
    'units_filter',
    'log10_filter',
    'get_atomic_symbol',
    'get_atomic_number',
    'sort_elements',
    'is_required',
    'sci_notation',
]


def units_filter(quantity, unit):
    """
    Convert quantity to given units and extract value
    """
    if not isinstance(quantity, u.Quantity):
        raise u.UnitsError(
            f'Value must be a quantity with units compatible with {unit}')
    return quantity.to(unit).value


def log10_filter(value):
    return np.log10(value)


def get_atomic_symbol(element):
    if type(element) is str:
        element = element.capitalize()
    return plasmapy.particles.atomic_symbol(element).lower()


def get_atomic_number(element):
    if type(element) is str:
        element = element.capitalize()
    return plasmapy.particles.atomic_number(element)


def sort_elements(elements):
    return sorted(elements, key=get_atomic_number)


def is_required(value):
    if isinstance(value, Undefined):
        raise MissingParameter('Parameter required for configuration')
    else:
        return value


def sci_notation(value, sig_figs=8):
    """
    Print value in scientific notation with number of signficant
    digits specified by `sig_figs`. This filter is required as
    there are known issues when supplying some quantities
    (e.g. loop length) at arbitrarily high precision
    """
    format_str = f'1.{sig_figs}e'
    return f'{value:{format_str}}'
