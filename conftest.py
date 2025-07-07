"""
Test rendering of templates for configuration and header files
"""
import astropy.units as u
import pytest

import pydrad.configure
from pydrad.configure.util import get_clean_hydrad


def pytest_addoption(parser):
    parser.addoption('--hydrad-dir', action='store', default=None)


@pytest.fixture(scope='session')
def hydrad_clean(tmpdir_factory, request):
    # Returns a local path to a copy of the HYDRAD code
    # If a path is not passed as a command line argument to pytest,
    # a new copy is cloned from GitHub
    hydrad_dir = request.config.getoption('--hydrad-dir')
    if hydrad_dir is None:
        hydrad_dir = tmpdir_factory.mktemp('hydrad_tmp_clean')
        get_clean_hydrad(hydrad_dir, base_path=None, from_github=True, overwrite=True)
    return hydrad_dir


def get_configuration_dict():
    return {
        'general': {
            'footpoint_height': 5.e+08*u.cm,
            'loop_inclination': 0.*u.deg,
            'force_single_fluid': False,
            'heat_flux_timestep_limit':  1.e-10*u.s,
            'logging_frequency': 1000,
            'minimum_collisional_coupling_timescale':  0.01*u.s,
            'output_interval':  1.*u.s,
            'heat_flux_limiting_coefficient': 0.167,
            'use_kinetic_model': False,
            'write_file_equation_terms': True,
            'write_file_hydrogen_level_populations': True,
            'write_file_ion_populations': True,
            'write_file_physical': True,
            'write_file_timescales': True,
            'loop_length':  90.*u.Mm,
            'total_time':  2.*u.s,
        },
        'grid': {
            'adapt': True,
            'adapt_every_n_time_steps': 1000,
            'enforce_conservation': True,
            'linear_restriction': True,
            'maximum_cell_width':  0.5*u.Mm,
            'maximum_fractional_difference': 0.1,
            'maximum_refinement_level': 12,
            'initial_refinement_level': 10,
            'maximum_variation': 0.1,
            'minimum_delta_s':  1.*u.cm,
            'minimum_fractional_difference': 0.05,
            'refine_on_density': True,
            'refine_on_electron_energy': True,
            'refine_on_hydrogen_energy': True,
        },
        'heating': {
            'alfven_wave': False,
            'background': {'use_initial_conditions': True},
            'beam': False,
            'events': [],
            'electron_heating': 1.0,
        },
        'initial_conditions': {
            'footpoint_density':  1.e+12*u.cm**(-3),
            'footpoint_temperature':  10000.*u.K,
            'heating_range_fine_tuning': 10000.0,
            'heating_range_lower_bound':  1.e-08*u.erg / (u.cm**3*u.s),
            'heating_range_step_size': 0.01,
            'heating_range_upper_bound':  100.*u.erg / (u.cm**3*u.s),
            'isothermal': False,
            'use_poly_fit_gravity': False,
            'use_poly_fit_magnetic_field': False,
            'heating_location':  45.*u.Mm,
            'heating_scale_height':  1.e+300*u.cm
        },
        'radiation': {
            'abundance_dataset': 'asplund',
            'decouple_ionization_state_solver': False,
            'density_dependent_rates': False,
            'elements_equilibrium': [],
            'elements_nonequilibrium': [],
            'emissivity_dataset': 'chianti_v7',
            'nlte_chromosphere': False,
            'optically_thick_radiation': False,
            'ranges_dataset': 'ranges',
            'rates_dataset': 'chianti_v7',
            'use_power_law_radiative_losses': True
        },
        'solver': {
            'cutoff_ion_fraction': 1e-15,
            'epsilon': 0.01,
            'epsilon_d': 0.1,
            'epsilon_r': 1.8649415311920072,
            'maximum_optically_thin_density': 1.e+12*u.cm**(-3),
            'minimum_radiation_temperature':  20000.*u.K,
            'minimum_temperature':  10000.*u.K,
            'safety_advection': 1.0,
            'safety_atomic': 1.0,
            'safety_conduction': 1.0,
            'safety_radiation': 0.1,
            'safety_viscosity': 1.0,
            'timestep_increase_limit': 0.05,
            'zero_over_temperature_interval':  500.*u.K,
        }
    }


@pytest.fixture(scope='function')
def configuration():
    return pydrad.configure.Configure(get_configuration_dict(), freeze_date=True)


@pytest.fixture(scope='session')
def hydrad(tmpdir_factory, hydrad_clean):
    # NOTE: purposefully reconstructing the configuration here, rather than using the fixture
    # above in order to keep configuration at the function scope level
    configuration = pydrad.configure.Configure(get_configuration_dict(), freeze_date=True)
    # Run a HYDRAD simulation for the given configuration and return the path to the directory
    # containing the results.
    hydrad_tmp = tmpdir_factory.mktemp('hydrad_tmp')
    configuration.setup_simulation(hydrad_tmp, hydrad_clean, overwrite=True)
    pydrad.configure.util.run_shell_command(hydrad_tmp / 'HYDRAD.exe')
    return hydrad_tmp
