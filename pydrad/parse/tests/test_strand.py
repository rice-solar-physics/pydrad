"""
Test parsing HYDRAD results
"""
import h5py
import pytest
import astropy.units as u

from pydrad.parse import Strand


VAR_NAMES = [
    'coordinate',
    'grid_edges',
    'grid_widths',
    'grid_centers',
    'electron_temperature',
    'ion_temperature',
    'electron_density',
    'ion_density',
    'electron_pressure',
    'ion_pressure',
    'velocity',
    'electron_heating_term',
    'hydrogen_heating_term',
    'radiative_loss_term',
    'sound_speed',
    'electron_conduction',
    'ion_conduction',
]


@pytest.fixture
def strand(hydrad):
    return Strand(hydrad)


def test_parse_initial_conditions(strand):
    assert hasattr(strand, 'initial_conditions')


@pytest.mark.parametrize('quantity', VAR_NAMES)
def test_has_quantity(strand, quantity):
    for p in strand:
        assert hasattr(p, quantity)
        assert isinstance(getattr(p, quantity), u.Quantity)


def test_time_arrays_same(hydrad, strand):
    """Check that reading time arrays different ways yields same answer"""
    strand2 = Strand(hydrad, read_from_cfg=True)
    assert u.allclose(strand.time, strand2.time, rtol=0.0, atol=1e-2*u.s)


def test_to_hdf5(strand, tmp_path):
    filename = tmp_path / 'hydrad_results.h5'
    strand.to_hdf5(filename, *VAR_NAMES)
    with h5py.File(filename, 'r') as hf:
        t = u.Quantity(hf['time'], hf['time'].attrs['unit'])
        assert (strand.time == t).all()
        for i, _ in enumerate(strand.time):
            p = strand[i]
            for v in VAR_NAMES:
                ds = hf[f'index{i}/{v}']
                assert (getattr(p, v) == u.Quantity(ds, ds.attrs['unit'])).all()


def test_emission_measure(strand):
    em, bins = strand[0].column_emission_measure()
    assert isinstance(em, u.Quantity)
    assert isinstance(bins, u.Quantity)
    assert len(bins) == len(em) + 1
