"""
Unit tests for reader functions
"""
import pytest

from pydrad.parse.util import read_ine_file, read_phy_file


@pytest.fixture
def profile(strand):
    return strand[0]


def test_read_phy_file(profile):
    tab = read_phy_file(profile._phy_filename)
    assert len(tab.colnames) == 11
    assert tab['coordinate'].shape == profile.grid_centers.shape


def test_read_ine_file(profile):
    tab = read_ine_file(profile._ine_filename, profile.grid_centers.shape[0])
    assert len(tab.colnames) == 9  # NEI elements modeled are H and C: Z_H+1 + Z_C+1=9
    assert tab['hydrogen_1'].shape == profile.grid_centers.shape
