"""
===================================
Configure a Beam Heating Simulation
===================================

This example demonstrates a typical configuration of a HYDRAD
simulation with electron beam heating.
"""
import pathlib
import tempfile

import astropy.units as u

from pydrad.configure import Configure
from pydrad.configure.data import get_defaults
from pydrad.configure.util import get_clean_hydrad

#################################################################
# HYDRAD can be used for simulations of solar flares, in addition
# to active region loops.  Since flares have significantly stronger
# energy release, however, some of the parameters need to be set
# more stringently than might be done for e.g. nanoflare heating.
# This example will show how to set up a reasonable flare simulation.
#
# Let's start with the defaults:
config = get_defaults()

# Configure the basic parameters of the simulation first.  Let's assume a
# loop length of 50 Mm, and run the simulation for one hour of simulation time.
config['general']['loop_length'] = 50 * u.Mm
config['general']['total_time'] = 1 * u.h

# We typically assume that the loop is heated by electron beam heating,
# where electrons are accelerated to tens of keV in energy in the corona,
# and stream towards the chromosphere, where they deposit that energy.
#
# With beam heating, we typically want to use a slightly more detailed
# chromosphere (rather than assuming constant temperature).
# This is for two reasons: (1) it needs to be dense enough to stop the
# beam, and (2) the energy deposition will produce a more realistic
# chromospheric evaporation.
# To do this, we turn on optically thick radiation, which will calculate
# optically thick radiative losses in the chromosphere, produced by
# hydrogen, magnesium, and calcium (based on the formulation by
# Carlsson & Leenaarts 2012).
config['radiation']['optically_thick_radiation'] = True

# This will set the chromospheric temperature profile to the so-called
# VAL C model.  For consistency with that model, we also
# need to set a few boundary conditions
config['general']['footpoint_height'] = 2.26 * u.Mm
config['initial_conditions']['footpoint_temperature'] = 24000 * u.K
config['initial_conditions']['footpoint_density'] = 4.2489e9 * u.cm**(-3)
config['solver']['minimum_radiation_temperature'] = 24000 * u.K
config['solver']['minimum_temperature'] = 4170 * u.K

# There is one further option for the chromosphere.  It is not generally
# recommended, but produces more accurate electron densities.  This will
# solve an approximation to radiative transfer for hydrogen with the caveat that
# it will slow the code by well over an order of magnitude.  It is most
# useful for users who wish to synthesize chromospheric line profiles.
# Change the value to True if you would like to try using it.
config['radiation']['nlte_chromosphere'] = False

# This takes care of the chromosphere, but we should also set up the radiation
# calculation in the corona (the optically thin losses).  Let's use the 15 most
# abundant elements in the Sun for this calculation
config['radiation']['use_power_law_radiative_losses'] = False
config['radiation']['elements_equilibrium'] = ['H','He','C','N','O','Ne','Na','Mg','Al','Si','S','Ar','Ca','Fe','Ni']

# We still need to set up the electron beam itself.
# Let's assume a constant heating for 10 seconds, using a moderate
# energy flux, typical low energy cut-off of 15 keV, and typical spectral index of 5.
config['heating']['beam'] = [
            {'time_start': 10.0*u.s,
            'flux': 3e10*u.erg/(u.cm**2)/u.s,
            'cut_off': 15.0*u.keV,
            'index': 5.0},
    ]

# Finally, let's set parameters to ensure the simulation runs smoothly.
# There are a lot of numerical parameters, so we'll explain what these mean.
#
# First, let's consider the time-stepping.  HYDRAD solves its equations explicitly,
# meaning that the CFL condition must be met at all times for numerical stability.
# This means that the time-step, :math:`dt`, must always be sufficiently small to produce
# an accurate numerical solution.
#
# Let's first set "safety" parameters for the radiation and conduction time-scales,
# which simply reduces their values to make it more likely that the CFL condition is met.
config['solver']['safety_radiation'] = 0.1
config['solver']['safety_conduction'] = 0.1
# Since conduction can be extremely limiting, we also set an effective floor:
config['general']['heat_flux_timestep_limit'] = 1e-6 * u.s

# Next, we want to make sure that the grid is spatially resolved
# We are telling the code to refine the grid every time step, and that
# any grid cell can be refined up to 12 times.  For extremely strong heating events,
# you might consider increasing this to 14.
config['solver']['adapt_every_n_time_steps'] = 1
config['solver']['initial_refinement_level'] = 12
config['solver']['maximum_refinement_level'] = 12

# Finally, we change a few parameters from their defaults, only for stability
#   Do not check for numerical precision errors in the conservation of energy:
config['solver']['enforce_conservation'] = False
#   Do not refine the grid on hydrogen energy:
config['solver']['refine_on_hydrogen_energy'] = False


#################################################################
# As in the basic tutorial, we now use the configuration defined here
# to set up the simulation by writing configuration files.
c = Configure(config)
tmpdir = pathlib.Path(tempfile.mkdtemp())  # Change this to wherever you want to save the unedited HYDRAD copy
hydrad_clean = tmpdir / 'hydrad-clean'
get_clean_hydrad(hydrad_clean, from_github=True)
c.setup_simulation(tmpdir / 'test-run', hydrad_clean)  # The location of the simulation itself

# Let's finally save the configuration, in case we want to reuse it or modify it later.
asdf_config = tmpdir / 'beam_config.asdf'

# This will have set up the simulation, which can now be run!  Since
# beam heating simulations can be somewhat slower, it is recommended
# to run this by hand in a terminal.  Also, if using Mac, be sure to
# use the "caffeinate" command to keep the computer from sleeping.
