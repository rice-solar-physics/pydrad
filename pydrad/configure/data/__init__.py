"""
Data related to configuration
"""
import copy
import importlib

import asdf


def get_defaults():
    filename = importlib.resources.files(
        'pydrad' ) / 'configure/data/defaults.asdf'
    with asdf.open(filename) as af:
        return copy.deepcopy(dict(af.tree))
