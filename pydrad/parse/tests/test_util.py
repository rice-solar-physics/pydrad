"""
Unit tests for reader functions
"""
import plasmapy.particles

from pydrad.parse.util import read_ine_file, read_phy_file


def test_read_phy_file(strand):
    tab = read_phy_file(strand[0]._phy_filename)
    assert len(tab.colnames) == 11
    assert tab['coordinate'].shape == strand[0].grid_centers.shape


def test_read_ine_file(strand):
    tab = read_ine_file(strand[0]._ine_filename, strand[0].grid_centers.shape[0])
    n_columns = sum([plasmapy.particles.atomic_number(el)+1
                     for el in strand.config['radiation']['elements_nonequilibrium']])
    assert len(tab.colnames) == n_columns  # NEI elements modeled are H and C: Z_H+1 + Z_C+1=9
    assert tab['hydrogen_1'].shape == strand[0].grid_centers.shape
