"""
Test rendering of templates for configuration and header files
"""
import pytest
import astropy.units as u

import pydrad.configure


def pytest_addoption(parser):
    parser.addoption('--hydrad-dir', action='store', default=None)


@pytest.fixture
def hydrad_clean(request):
    return request.config.getoption('--hydrad-dir')


@pytest.fixture
def configuration_dict():
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
            'heat_electrons': True,
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


@pytest.fixture
def configuration(configuration_dict):
    return pydrad.configure.Configure(configuration_dict, freeze_date=True)


@pytest.fixture
def hydrad(tmp_path, configuration, hydrad_clean):
    if hydrad_clean is None:
        pytest.skip('Path to HYDRAD code not specified. Skipping those tests that require HYDRAD.')
    else:
        hydrad_tmp = tmp_path / 'hydrad_tmp'
        configuration.setup_simulation(hydrad_tmp, hydrad_clean)
        pydrad.configure.util.run_shell_command(
            ['./HYDRAD.exe'],
            hydrad_tmp,
        )
        return hydrad_tmp
