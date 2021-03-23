"""
Some tests for the configuration module
"""
import pytest
import numpy as np
import astropy.units as u

from pydrad.configure.util import MissingParameter

from . import assert_ignore_blanks


def test_collisions_header(configuration):
    collisions_header = f"""// ****
// *
// * #defines for configuring the shortest collisional coupling timescale
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Physics ****
#define MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE 0.01

// **** End of Physics ****"""
    assert_ignore_blanks(configuration.collisions_header, collisions_header)
    configuration.config['general']['force_single_fluid'] = True
    collisions_header = f"""// ****
// *
// * #defines for configuring the shortest collisional coupling timescale
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Physics ****
#define MINIMUM_COLLISIONAL_COUPLING_TIME_SCALE 0.01
#define FORCE_SINGLE_FLUID
// **** End of Physics ****"""
    assert_ignore_blanks(configuration.collisions_header, collisions_header)


def test_heating_config(configuration):
    # No events, no background
    configuration.config['heating']['background'] = {
        'location': 0*u.cm,
        'scale_height': 0*u.cm,
        'rate': 0.0*u.erg/(u.cm**3)/u.s 
    }
    heating_config = f"""0.00000000e+00 0.00000000e+00 0.0

0


Configuration file generated by pydrad on {configuration.date}"""
    assert configuration.heating_cfg == heating_config
    # No events, custom background
    configuration.config['heating']['background'] = {
        'location': 0*u.cm,
        'scale_height': 1e300*u.cm,
        'rate': 0.1*u.erg/(u.cm**3)/u.s 
    }
    heating_config = f"""0.0 1e+300 0.1

0


Configuration file generated by pydrad on {configuration.date}"""
    # Initial conditions not run
    configuration.config['heating']['background'] = {
        'use_initial_conditions': True,
    }
    with pytest.raises(AttributeError):
        configuration.heating_cfg
    # 1 event, no background
    configuration.config['heating']['background'] = {
        'location': 0*u.cm,
        'scale_height': 0*u.cm,
        'rate': 0.0*u.erg/(u.cm**3)/u.s 
    }
    configuration.config['heating']['events'] = [
        {'time_start': 0.*u.s,
         'rise_duration': 100.*u.s,
         'decay_duration': 100.*u.s,
         'total_duration': 200.*u.s,
         'location': 0*u.cm,
         'scale_height': 1e300*u.cm,
         'rate': 0.1*u.erg/(u.cm**3)/u.s}
    ]
    heating_config = f"""0.00000000e+00 0.00000000e+00 0.0

1

0.00000000e+00 1.00000000e+300 0.1 0.0 100.0 100.0 200.0

Configuration file generated by pydrad on {configuration.date}"""
    assert configuration.heating_cfg == heating_config


def test_heating_header(configuration):
    # Electron Heating
    header = f'''// ****
// *
// * #defines for configuring the heating model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****
#define ELECTRON_HEATING 1.0
#define HYDROGEN_HEATING 0.0
#define ELECTRON_HEATING_ONLY
#include "../../Radiation_Model/source/config.h"'''
    assert_ignore_blanks(configuration.heating_header, header)
    # Ion Heating
    configuration.config['heating']['electron_heating'] = 0.0
    header = f'''// ****
// *
// * #defines for configuring the heating model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****
#define ELECTRON_HEATING 0.0
#define HYDROGEN_HEATING 1.0
#define HYDROGEN_HEATING_ONLY
#include "../../Radiation_Model/source/config.h"'''
    assert_ignore_blanks(configuration.heating_header, header)
    # Mix between electron and ion Heating
    configuration.config['heating']['electron_heating'] = 0.5
    header = f'''// ****
// *
// * #defines for configuring the heating model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****
#define ELECTRON_HEATING 0.5
#define HYDROGEN_HEATING 0.5
#include "../../Radiation_Model/source/config.h"'''
    assert_ignore_blanks(configuration.heating_header, header)
    # Other Heating Models
    configuration.config['heating']['alfven_wave'] = True
    configuration.config['heating']['beam'] = True
    header = f'''// ****
// *
// * #defines for configuring the heating model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****
#define ELECTRON_HEATING 0.5
#define HYDROGEN_HEATING 0.5
#define BEAM_HEATING
#define ALFVEN_WAVE_HEATING
#include "../../Radiation_Model/source/config.h"'''
    assert_ignore_blanks(configuration.heating_header, header)


def test_hydrad_config(configuration):
    config = f"""Initial_Conditions/profiles/initial.amr
Initial_Conditions/profiles/initial.amr.gravity
2.0
1.0

Configuration file generated by pydrad on {configuration.date}"""
    assert configuration.hydrad_cfg == config
    configuration.config['general']['poly_fit_gravity'] = np.array([1, 2, 3])
    configuration.config['general']['poly_fit_magnetic_field'] = np.array([1, 2, 3])
    config = f"""Initial_Conditions/profiles/initial.amr
poly_fit.gravity
poly_fit.magnetic_field
2.0
1.0

Configuration file generated by pydrad on {configuration.date}"""
    assert configuration.hydrad_cfg == config
    configuration.config['general']['initial_amr_file'] = 'Results/profile10.amr'
    config = f"""Results/profile10.amr
poly_fit.gravity
poly_fit.magnetic_field
2.0
1.0

Configuration file generated by pydrad on {configuration.date}"""


def test_hydrad_header(configuration):
    header = f"""// ****
// *
// * #defines for configuring the hydrodynamic model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Output ****
#define WRITE_FILE_PHYSICAL
#define WRITE_FILE_ION_POPULATIONS
#define WRITE_FILE_HSTATE
#define WRITE_FILE_SCALES
#define WRITE_FILE_TERMS
#define OUTPUT_EVERY_N_TIME_STEPS 1000
// **** End of Output ****

// **** Physics ****
#include "../../Heating_Model/source/config.h"
#include "../../Radiation_Model/source/config.h"
#define HEAT_FLUX_LIMITING_COEFFICIENT 0.167
#define TIME_STEP_LIMIT 1e-10
#include "collisions.h"

// **** End of Physics ****

// **** Solver ****

#define SAFETY_RADIATION 0.1
#define SAFETY_CONDUCTION 1.0
#define SAFETY_ADVECTION 1.0
#define SAFETY_VISCOSITY 1.0
#define TIME_STEP_INCREASE_LIMIT 1.05

#define MINIMUM_RADIATION_TEMPERATURE 20000.0
#define ZERO_OVER_TEMPERATURE_INTERVAL 500.0
#define MINIMUM_TEMPERATURE 10000.0

// **** End of Solver ****

// **** Grid ****
#define MAX_REFINEMENT_LEVEL 12
#define INITIAL_REFINEMENT_LEVEL 10
#define ADAPT
#define ADAPT_EVERY_N_TIME_STEPS 1000
#define REFINE_ON_DENSITY
#define REFINE_ON_ELECTRON_ENERGY
#define REFINE_ON_HYDROGEN_ENERGY
#define MIN_FRAC_DIFF 0.05
#define MAX_FRAC_DIFF 0.1
#define LINEAR_RESTRICTION
#define ENFORCE_CONSERVATION
// **** End of Grid ****"""
    assert_ignore_blanks(configuration.hydrad_header, header)
    configuration.config['general']['poly_fit_magnetic_field'] = np.array([1, 2, 3])
    header = f"""// ****
// *
// * #defines for configuring the hydrodynamic model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Output ****
#define WRITE_FILE_PHYSICAL
#define WRITE_FILE_ION_POPULATIONS
#define WRITE_FILE_HSTATE
#define WRITE_FILE_SCALES
#define WRITE_FILE_TERMS
#define OUTPUT_EVERY_N_TIME_STEPS 1000
// **** End of Output ****

// **** Physics ****
#include "../../Heating_Model/source/config.h"
#include "../../Radiation_Model/source/config.h"
#define HEAT_FLUX_LIMITING_COEFFICIENT 0.167
#define TIME_STEP_LIMIT 1e-10
#include "collisions.h"
#define USE_POLY_FIT_TO_MAGNETIC_FIELD
// **** End of Physics ****

// **** Solver ****

#define SAFETY_RADIATION 0.1
#define SAFETY_CONDUCTION 1.0
#define SAFETY_ADVECTION 1.0
#define SAFETY_VISCOSITY 1.0
#define TIME_STEP_INCREASE_LIMIT 1.05

#define MINIMUM_RADIATION_TEMPERATURE 20000.0
#define ZERO_OVER_TEMPERATURE_INTERVAL 500.0
#define MINIMUM_TEMPERATURE 10000.0

// **** End of Solver ****

// **** Grid ****
#define MAX_REFINEMENT_LEVEL 12
#define INITIAL_REFINEMENT_LEVEL 10
#define ADAPT
#define ADAPT_EVERY_N_TIME_STEPS 1000
#define REFINE_ON_DENSITY
#define REFINE_ON_ELECTRON_ENERGY
#define REFINE_ON_HYDROGEN_ENERGY
#define MIN_FRAC_DIFF 0.05
#define MAX_FRAC_DIFF 0.1
#define LINEAR_RESTRICTION
#define ENFORCE_CONSERVATION
// **** End of Grid ****"""
    assert_ignore_blanks(configuration.hydrad_header, header)


def test_initial_conditions_config(configuration):
    configuration.config['initial_conditions']['isothermal'] = True
    # Isothermal case
    config = f"""Initial_Conditions/profiles/initial.amr

9.00000000e+09
0.0
5.00000000e+08

1.000000e+04

1.0000000e+12

4.50000000e+09
1.00000000e+300
0.0
0.0
1.0
1.0

Configuration file generated by pydrad on {configuration.date}"""
    assert configuration.initial_conditions_cfg == config
    # Hydrostatic equilibrium
    configuration.config['initial_conditions']['isothermal'] = False
    configuration.config['initial_conditions']['heating_range_lower_bound'] = 1e-8*u.erg/(u.cm**3)/u.s
    configuration.config['initial_conditions']['heating_range_upper_bound'] = 1e2*u.erg/(u.cm**3)/u.s
    configuration.config['initial_conditions']['heating_range_step_size'] = 0.01
    configuration.config['initial_conditions']['heating_range_fine_tuning'] = 10000
    config = f"""Initial_Conditions/profiles/initial.amr

9.00000000e+09
0.0
5.00000000e+08

1.000000e+04

1.0000000e+12

4.50000000e+09
1.00000000e+300
-8.0
2.0
0.01
10000

Configuration file generated by pydrad on {configuration.date}"""
    assert configuration.initial_conditions_cfg == config


def test_initial_conditions_header(configuration):
    configuration.config['initial_conditions']['isothermal'] = False
    header = f"""// ****
// *
// * #defines for configuring the hydrostatic model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Output ****
// **** End of Output ****

// **** Physics ****
#include "../../Radiation_Model/source/config.h"

// **** Solver ****
#define EPSILON 0.01

// **** Grid ****
#define ADAPT
#define MIN_CELLS 180
#define MAX_CELLS 737280
#define MAX_REFINEMENT_LEVEL 12
#define INITIAL_REFINEMENT_LEVEL 10
#define MIN_DS 1.00000000e+00
#define MAX_VARIATION 1.1"""
    assert_ignore_blanks(configuration.initial_conditions_header, header)
    configuration.config['initial_conditions']['isothermal'] = True
    configuration.config['general']['poly_fit_gravity'] = np.array([1, 2, 3])
    configuration.config['initial_conditions']['use_poly_fit_gravity'] = True
    header = f"""// ****
// *
// * #defines for configuring the hydrostatic model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Output ****
// **** End of Output ****

// **** Physics ****
#include "../../Radiation_Model/source/config.h"
#define ISOTHERMAL
#define USE_POLY_FIT_TO_GRAVITY
#define POLY_FIT_TO_GRAVITY_FILE "poly_fit.gravity"

// **** Solver ****
#define EPSILON 0.01

// **** Grid ****
#define ADAPT
#define MIN_CELLS 180
#define MAX_CELLS 737280
#define MAX_REFINEMENT_LEVEL 12
#define INITIAL_REFINEMENT_LEVEL 10
#define MIN_DS 1.00000000e+00
#define MAX_VARIATION 1.1"""
    assert_ignore_blanks(configuration.initial_conditions_header, header)


def test_radiation_header(configuration):
    # Test all true
    configuration.config['radiation']['elements_nonequilibrium'] = ['iron']
    configuration.config['radiation']['decouple_ionization_state_solver'] = True
    configuration.config['radiation']['density_dependent_rates'] = True
    configuration.config['radiation']['optically_thick_radiation'] = True
    configuration.config['radiation']['nlte_chromosphere'] = True
    configuration.config['radiation']['minimum_density_limit'] = 1e12*u.cm**(-3)
    configuration.config['solver']['cutoff_ion_fraction'] = 1e-6
    header = f"""// ****
// *
// * #defines for configuring the radiation model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Physics ****
#define NON_EQUILIBRIUM_RADIATION
#define DECOUPLE_IONISATION_STATE_SOLVER
#define USE_POWER_LAW_RADIATIVE_LOSSES
#define DENSITY_DEPENDENT_RATES
#define OPTICALLY_THICK_RADIATION
#define NLTE_CHROMOSPHERE
#define MIN_DENSITY_LIMIT 1000000000000.0
#include "../../HYDRAD/source/collisions.h"
// **** End of Physics ****

// **** Solver ****
#define MAX_OPTICALLY_THIN_DENSITY 1000000000000.0
#define SAFETY_ATOMIC 1.0
#define CUTOFF_ION_FRACTION 1e-06
#define EPSILON_D 0.1
#define EPSILON_R 1.8649415311920072
// **** End of Solver ****"""
    assert_ignore_blanks(configuration.radiation_header, header)
    # Test missing minimum density limit raises error
    del configuration.config['radiation']['minimum_density_limit']
    with pytest.raises(MissingParameter):
        configuration.radiation_header
    # Test all false
    configuration.config['radiation']['use_power_law_radiative_losses'] = False
    configuration.config['radiation']['density_dependent_rates'] = False
    configuration.config['radiation']['optically_thick_radiation'] = False
    configuration.config['radiation']['nlte_chromosphere'] = False
    header = f"""// ****
// *
// * #defines for configuring the radiation model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Physics ****

#define NON_EQUILIBRIUM_RADIATION
#define DECOUPLE_IONISATION_STATE_SOLVER




#include "../../HYDRAD/source/collisions.h"
// **** End of Physics ****

// **** Solver ****
#define MAX_OPTICALLY_THIN_DENSITY 1000000000000.0
#define SAFETY_ATOMIC 1.0
#define CUTOFF_ION_FRACTION 1e-06
#define EPSILON_D 0.1
#define EPSILON_R 1.8649415311920072
// **** End of Solver ****"""
    assert_ignore_blanks(configuration.radiation_header, header)
    # Test that removing non-equilibrium radiation removes the
    # decouple flag
    configuration.config['radiation']['elements_nonequilibrium'] = []
    header = f"""// ****
// *
// * #defines for configuring the radiation model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by pydrad on {configuration.date}
// *
// ****

// **** Physics ****





#include "../../HYDRAD/source/collisions.h"
// **** End of Physics ****

// **** Solver ****
#define MAX_OPTICALLY_THIN_DENSITY 1000000000000.0
#define SAFETY_ATOMIC 1.0
#define CUTOFF_ION_FRACTION 1e-06
#define EPSILON_D 0.1
#define EPSILON_R 1.8649415311920072
// **** End of Solver ****"""


def test_radiation_config_equilibrium(configuration):
    configuration.config['radiation']['elements_equilibrium'] = [
        1,
        'He',
        'iron'
    ]
    config = f"""ranges
chianti_v7
asplund
chianti_v7
3
h
1
he
2
fe
26

Configuration file generated by pydrad on {configuration.date}"""
    assert configuration.radiation_equilibrium_cfg == config


def test_radiation_config_nonequilibrium(configuration):
    configuration.config['radiation']['elements_nonequilibrium'] = ['iron']
    config = f"""ranges
chianti_v7
asplund
chianti_v7
1
fe
26

Configuration file generated by pydrad on {configuration.date}"""
    assert configuration.radiation_nonequilibrium_cfg == config


def test_unit_conversion_error(configuration):
    configuration.config['general']['total_time'] = 1 * u.cm
    with pytest.raises(u.UnitsError):
        configuration.hydrad_cfg