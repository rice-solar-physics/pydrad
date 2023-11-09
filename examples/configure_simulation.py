"""
======================
Configure a Simulation
======================

This example demonstrates the various ways to configure a HYDRAD
simulation.
"""
import pathlib
import tempfile

import astropy.units as u

from pydrad.configure import Configure
from pydrad.configure.data import get_defaults
from pydrad.configure.util import get_clean_hydrad

#################################################################
# HYDRAD requires setting many different configuration files, making it
# difficult and tedious to configure a simulation by hand. Fortunately,
# the ``pydrad`` package makes it easy to quickly configure a new
# simulation using Python.
#
# Rather than editing each configuration file individually, an entire
# simulation can be configured using a single Python `dict`.
# Below is an example of a dictionary for configuring a simulation of a 80
# Mm loop lasting 5000 s heated by a single 200 s nanoflare
# solved on an adaptive grid.
# A complete list of configuration parameters can be found in
# the `configuration-tables`_ page.
config_dict = {
    'general': {
        'loop_length': 80*u.Mm,
        'total_time': 5000*u.s,
        'footpoint_height': 5*u.Mm,
        'loop_inclination': 0*u.deg,
        'force_single_fluid': False,
        'heat_flux_timestep_limit':  1e-10*u.s,
        'logging_frequency': 1000,
        'minimum_collisional_coupling_timescale':  0.01*u.s,
        'output_interval':  10*u.s,
        'heat_flux_limiting_coefficient': 0.167,
        'use_kinetic_model': False,
        'write_file_equation_terms': True,
        'write_file_hydrogen_level_populations': True,
        'write_file_ion_populations': True,
        'write_file_physical': True,
        'write_file_timescales': True,
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
        'electron_heating': 1.0,
        'events': [
            {'time_start': 0.*u.s,
             'rise_duration': 100*u.s,
             'decay_duration': 100*u.s,
             'total_duration': 200*u.s,
             'location': 40*u.Mm,
             'scale_height': 1e300 * u.cm,
             'rate': 0.1 * u.erg/u.s/(u.cm**3)},
        ],
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

#################################################################
# The next step is to pass it to the `~pydrad.configure.Configure`
# class which can parse it and print the needed configuration files,
c = Configure(config_dict)

#################################################################
# We can preview some of the different configuration files that setup the
# HYDRAD simulation. This is often useful for debugging.
print(c.initial_conditions_header)
print(c.initial_conditions_cfg)
print(c.hydrad_header)

#################################################################
# To actually setup a simulation, you'll first need to obtain a
# copy of HYDRAD. You can use a copy you already have locally or
# use the following convenience function to grab the most recent
# version from GitHub
tmpdir = pathlib.Path(tempfile.mkdtemp())  # Change this to wherever you want to save your clean HYDRAD copy
hydrad_clean = tmpdir / 'hydrad-clean'
get_clean_hydrad(hydrad_clean, from_github=True)

#################################################################
# The way HYDRAD works, a complete copy of the source code is required
# for every simulation run. When setting up a simulation, we pass both
# the path to our "clean" copy (which will remain unchanged) and the path
# to our new simulation. This function will write all of the needed
# aforementioned configuration files and run the hydrostatic solver which
# will provide the initial conditions for the actual simulation.
c.setup_simulation(tmpdir / 'test-run', hydrad_clean)

#################################################################
# To avoid having to repeatedly setup the configuration dictionary,
# the simulation setup can be serialized to an `ASDF
# file <https://asdf.readthedocs.io>`__, a human-readable, structured
# plain text file in the YAML format. These ASDF files are structured just
# like the config directory and can be easily read and written.
#
# To save the configuration to disk and then load it back into a `dict`,
asdf_config = tmpdir / 'test_config.asdf'
c.save_config(asdf_config)
config_from_disk = Configure.load_config(asdf_config)
print(config_from_disk)

#################################################################
# Additionally, `~pydrad` also includes a default set of configuration
# parameters. These are often a useful starting point when creating
# your own configuration.
config_default = get_defaults()
print(config_default)

#################################################################
# If you make local modifications to the HYDRAD code, you may need
# configuration options not in the included templates. To use custom
# configuration options, you can inject your own modified templates of the
# configuration files which take advantage of any custom options.
#
# To see the available templates,
print(c.templates)

#################################################################
# Say we want to add an option, ``MY_NEW_PARAM`` to the ``collisions.h``
# file. To get the current unrendered template,
print(c.get_raw_template('collisions.h'))

#################################################################
# We can then create a new template with our new value,
new_collisions = """// ****
// *
// * #defines for configuring the shortest collisional coupling timescale
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {{ date }}
// *
// ****

// **** Physics ****
#define MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE {{ general.minimum_collisional_coupling_timescale | units_filter('s') }}
{% if general.force_single_fluid -%}#define FORCE_SINGLE_FLUID{%- endif %}
// **** End of Physics ****
#define MY_NEW_PARAM {{ general.my_new_param }}"""

#################################################################
# add our new parameter to the configuration directory,
# and then pass the template to the ``Configure`` object,
config_dict['general']['my_new_param'] = 100
c_custom = Configure(config_dict, templates={'collisions.h': new_collisions})
print(c_custom.collisions_header)
