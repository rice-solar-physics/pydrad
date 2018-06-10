"""
Custom Jinja2 filters
"""
import numpy as np
import astropy.units as u
import plasmapy.atomic

from .util import MissingParameter

__all__ = ['units_filter', 'log10_filter', 'get_atomic_symbol', 'get_atomic_number',
           'sort_elements', 'is_required']


def units_filter(quantity, unit):
    """
    Convert quantity to given units and extract value
    """
    if not isinstance(quantity, u.Quantity):
        raise u.UnitsError(f'Value must be a quantity with units compatible with {unit}')
    return quantity.to(unit).value


def log10_filter(value):
    return np.log10(value)


def get_atomic_symbol(element):
    if type(element) is str:
        element = element.capitalize()
    return plasmapy.atomic.atomic_symbol(element).lower()


def get_atomic_number(element):
    if type(element) is str:
        element = element.capitalize()
    return plasmapy.atomic.atomic_number(element)


def sort_elements(elements):
    return sorted(elements, key=get_atomic_number)


def is_required(value):
    # Check for missing values that are not False booleans
    if not value and value is not False:
        raise MissingParameter('Parameter required for configuration')
    else:
        return value
