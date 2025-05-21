"""
Test parsing HYDRAD results
"""
import astropy.units as u
import h5py
import numpy as np
import pytest

from pydrad.parse import Strand

VAR_NAMES = [
    'coordinate',
    'grid_edges',
    'grid_widths',
    'grid_centers',
    'electron_temperature',
    'hydrogen_temperature',
    'electron_density',
    'hydrogen_density',
    'electron_pressure',
    'hydrogen_pressure',
    'velocity',
    'sound_speed',
    'electron_heat_flux',
    'hydrogen_heat_flux',
    'mass_drhobydt',
    'mass_advection',
    'momentum_drho_vbydt',
    'momentum_advection',
    'momentum_pressure_gradient',
    'momentum_gravity',
    'momentum_viscous_stress',
    'momentum_numerical_viscosity',
    'electron_dTEKEbydt',
    'electron_enthalpy',
    'electron_conduction',
    'electron_collisions',
    'electron_heating',
    'electron_radiative_loss',
    'electron_electric_field',
    'electron_ionization',
    'hydrogen_dTEKEbydt',
    'hydrogen_enthalpy',
    'hydrogen_conduction',
    'hydrogen_gravity',
    'hydrogen_collisions',
    'hydrogen_heating',
    'hydrogen_electric_field',
    'hydrogen_viscous_stress',
    'hydrogen_numerical_viscosity',
]


@pytest.fixture
def strand(hydrad):
    return Strand(hydrad)

@pytest.fixture
def strand_only_amr(hydrad):
    return Strand(hydrad,
                  read_phy=False,
                  read_ine=False,
                  read_trm=False,
                  read_hstate=False,
                  read_scl=False,
                  )


def test_parse_initial_conditions(strand):
    assert hasattr(strand, 'initial_conditions')


@pytest.mark.parametrize('quantity', VAR_NAMES)
def test_has_quantity(strand, quantity):
    for p in strand:
        assert hasattr(p, quantity)
        assert isinstance(getattr(p, quantity), u.Quantity)


@pytest.mark.parametrize('quantity', [
    'electron_temperature',
    'hydrogen_temperature',
    'electron_density',
    'hydrogen_density',
    'electron_pressure',
    'hydrogen_pressure',
    'velocity',
])
def test_no_amr_run_has_quantity(strand_only_amr, quantity):
    # check that Strand falls back to deriving thermodynamic quantities from
    # conserved quantities when we do not read the phy files
    for p in strand_only_amr:
        assert hasattr(p, quantity)
        assert isinstance(getattr(p, quantity), u.Quantity)


def test_time_arrays_same(hydrad, strand):
    """Check that reading time arrays different ways yields same answer"""
    strand2 = Strand(hydrad, read_from_cfg=True)
    print(strand.time)
    print(strand2.time)
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

def test_term_file_output(strand):
    for p in strand:
        # The electron energy equation's numerical viscosity term is always 0:
        assert u.allclose(
            p.electron_numerical_viscosity,
            np.zeros_like(p.electron_numerical_viscosity),
            rtol=0.0,
        )
        # The hydrogen energy equation's collision rate is never 0 at all positions:
        assert not u.allclose(
            p.hydrogen_collisions,
            np.zeros_like(p.hydrogen_collisions),
            rtol=0.0,
        )

def test_term_file_units(strand):
    assert strand[0].mass_advection.unit == u.Unit('g s-1 cm-3')
    assert strand[0].momentum_gravity.unit == u.Unit('dyne s-1 cm-3')
    assert strand[0].electron_viscous_stress.unit == u.Unit('erg s-1 cm-3')
    assert strand[0].hydrogen_collisions.unit == u.Unit('erg s-1 cm-3')
