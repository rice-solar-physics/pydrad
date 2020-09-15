"""
Test parsing HYDRAD results
"""
import pytest
import astropy.units as u

from pydrad.parse import Strand


@pytest.fixture
def strand(hydrad):
    return Strand(hydrad)


def test_parse_initial_conditions(strand):
    assert hasattr(strand, 'initial_conditions')

@pytest.mark.parametrize(
    'quantity',
    ['coordinate',
     'grid_edges',
     'grid_widths',
     'grid_centers',
     'electron_temperature',
     'ion_temperature',
     'electron_density',
     'ion_density',
     'electron_pressure',
     'ion_pressure',
     'velocity'
     ]
)
def test_has_quantity(strand, quantity):
    for p in strand:
        assert hasattr(p, quantity)
        assert isinstance(getattr(p, quantity), u.Quantity)


def test_emission_measure(strand):
    em, bins = strand[0].column_emission_measure()
    assert isinstance(em, u.Quantity)
    assert isinstance(bins, u.Quantity)
    assert len(bins) == len(em) + 1
