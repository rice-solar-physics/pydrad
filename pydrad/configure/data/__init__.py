"""
Data related to configuration
"""
import copy
import pkg_resources

import asdf


def get_defaults():
    filename = pkg_resources.resource_filename(
        'pydrad',
        'configure/data/defaults.asdf',
    )
    with asdf.open(filename) as af:
        return copy.deepcopy(dict(af.tree))
