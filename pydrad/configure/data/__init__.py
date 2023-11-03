"""
Data related to configuration
"""
import copy

import asdf
import pkg_resources


def get_defaults():
    filename = pkg_resources.resource_filename(
        'pydrad',
        'configure/data/defaults.asdf',
    )
    with asdf.open(filename) as af:
        return copy.deepcopy(dict(af.tree))
