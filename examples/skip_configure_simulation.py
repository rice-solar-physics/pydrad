"""
======================
Configure a Simulation
======================

This example demonstrates the various ways to configure a HYDRAD
simulation.
"""
import astropy.units as u
from pydrad import configure

#################################################################
# HYDRAD requires setting many different configuration files, making it
# difficult and tedious to configure a simulation by hand. Fortunately,
# the ``pydrad`` package makes it easy to quickly configure a new
# simulation using Python.
# 
# Rather than editing each configuration file individually, an entire
# simulation can be configured using a single `Python
# dictionary <https://docs.python.org/3/library/stdtypes.html#dict>`__.
# Below is an example of a dictionary for configuring a simulation of a 80
# Mm loop lasting 5000 s with a single heating pulse using an adaptive
# grid. Many of the configuration options have been excluded for brevity.
# A complete list of configuration parameters can be found in the `section
# below <#configuration-parameters>`__.
config_dict = {
    'general': {
        'total_time': 5e3 * u.s,
        'loop_length': 80 * u.Mm,
        'footpoint_height': 5e8 * u.cm,
        # ... other config parameters ...
    },
    'initial_conditions': {
        'footpoint_temperature': 2e4 * u.K,
        'footpoint_density': 1e11 * u.cm**(-3),
        'isothermal': False,
        # ... other config parameters ...
    },
    'radiation': {
        'use_power_law_radiative_losses': True,
        # ... other config parameters ...
    },
    'heating': {
        'heat_electrons': True,
        'events': [
            {'time_start': 0.*u.s, 'rise_duration': 100*u.s,
            'decay_duration': 100*u.s, 'total_duration': 200*u.s,
            'location': 4e9*u.cm, 'scale_height': 1e300 * u.cm,
            'rate': 0.1 *u.erg/u.s/(u.cm**3), },
        ],
    },
    'solver': {
        'safety_radiation': 1.0,
        'safety_conduction': 1.0,
        'safety_advection': 1.0,
        'safety_atomic': 1.0,
        # ... other grid parameters ...
    },
    'grid': {
        'adapt': True,
        # ... other grid parameters ...
    }
} 

#################################################################
# The next step is to pass it to the ``pydrad`` class which can parse it
# and print the needed configuration files,
c = configure.Configure(config_dict)

#################################################################
# We can preview the different configuration files that setup the
# HYDRAD simulation,
print(c.initial_conditions_header)
print(c.initial_conditions_cfg)
print(c.hydrad_header)

#################################################################
# To print all configuration files, run the initial conditions, and copy
# all of this to a new location,
c.setup_simulation('/path/to/simulation/dir/new_hydrad_sim',
                      base_path='/path/to/clean/HYDRAD')

#################################################################
# This will create all of the needed input files from the options in
# ``config``, compile the initial conditions code, run
# ``Initial_Conditions.exe``, compile the main HYDRAD, and copy it all to
# the directory ``/path/to/simulation/dir/new_hydrad_sim``.


#################################################################
# HYDRAD requires a lot of configuration options and it can be annoying to
# have put them all in a configuration file. To avoid this, you can load a
# default configuration from an `ASDF
# file <https://asdf.readthedocs.io>`__, a human-readable, structured
# plain text file in the YAML format. These ASDF files are structured just
# like the config directory and can be easily read and written.
# 
# To load a configuration from a file and then save back to disk,
config = configure.Configure.load_config('/path/to/config/defaults.asdf')
c = configure.Configure(config)
c.save_config('/path/to/config/my_config.asdf')

#################################################################
# An example default configuration file can be found in the root of the
# `pydrad repository <https://github.com/wtbarnes/pydrad>`__.

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
config['general']['my_new_param'] = 100
c_new = Configure(config, templates={'collisions.h': new_collisions})
print(c_new.collisions_header)
