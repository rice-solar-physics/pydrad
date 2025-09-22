"""
Data related to configuration
"""
import asdf
import copy
import importlib


def get_defaults():
    filename = importlib.resources.files(
        'pydrad' ) / 'configure/data/defaults.asdf'
    with asdf.open(filename) as af:
        return copy.deepcopy(dict(af.tree))
