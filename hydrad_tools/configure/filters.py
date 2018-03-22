"""
Custom Jinja2 filters
"""
import numpy as np
import astropy.units as u
import plasmapy.atomic

__all__ = ['units_filter', 'log10_filter', 'get_atomic_symbol', 'get_atomic_number',
           'sort_elements']


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