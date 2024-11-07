"""
==================================
Parse and Visualize HYDRAD Results
==================================

This example demonstrates how to configure a simple hydrad simulation and
read in and visualize the outputs as a function of time and field-aligned
spatial coordinate.
"""
import pathlib
import tempfile

import astropy.units as u
from astropy.visualization import ImageNormalize, LogStretch

from pydrad.configure import Configure
from pydrad.configure.data import get_defaults
from pydrad.configure.util import get_clean_hydrad, run_shell_command
from pydrad.parse import Strand

tmpdir = pathlib.Path(tempfile.mkdtemp())  # Change to wherever you want to save your clean HYDRAD copy


#################################################################
# HYDRAD prints all results to the ``Results`` directory in the main code
# directory. These filenames are in the format ``profile{index}.{ext}``
# where ``index`` is the timestep number and ``ext`` can be one of the
# following filetypes,
#
# -  ``phy`` – main physics results file containing temperature, density,
#    etc.
# -  ``amr`` – data related to the adaptive mesh
# -  ``trm`` – equation terms at each grid cell
#
# Parsing all of these files can be difficult. `pydrad` provides a
# convenient object, `~pydrad.parse.Strand`, for easily accessing the
# data associated with a particular HYDRAD run.
#
# First, we'll setup a simple HYDRAD run and run the simulation.
# Grab a new copy of HYDRAD
hydrad_clean = tmpdir / 'hydrad-clean'
get_clean_hydrad(hydrad_clean, from_github=True)

#################################################################
# We'll start with the default configuration and setup a simple
# simulation for a 50 Mm loop lasting 500 s in which there are no
# heating events such that the loop is maintained in static
# equilibrium.
config = get_defaults()
config['general']['total_time'] = 500 * u.s
config['general']['output_interval'] = 10 * u.s
config['general']['loop_length'] = 50 * u.Mm
config['heating']['background']['use_initial_conditions'] = True
config['heating']['events'] = []

################################################################
# To reduce the commputation time and memory footprint of our
# example simulation, we reduce the degree to which our spatial
# grid is refined.
# In general, this is something that should be done cautiously.
# See
# `Bradshaw and Cargill (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJ...770...12B>`__
# for more details on appropriately resolving gradients in the transition
# region. This is an especially important consideration in impulsive
# heating scenarios.
config['grid']['initial_refinement_level'] = 6
config['grid']['maximum_refinement_level'] = 6

#################################################################
# Next, we'll setup and run the simulation. We'll run HYDRAD
# directly from Python, but this is easily done via the command
# line as well. This can take a few minutes.
c = Configure(config)
hydrad_results = tmpdir / 'steady-run'
c.setup_simulation(hydrad_results, hydrad_clean)
run_shell_command(hydrad_results / 'HYDRAD.exe')

#################################################################
# To parse the results of a simulation, we create a `~pydard.parse.Strand`
# object which holds information about the entire simulation.
s = Strand(hydrad_results)
print(s)

#################################################################
# We can get basic information about the simulation such as the
# time and loop length
print(s.time)
print(s.loop_length)

#################################################################
# Indexing the `~pydrad.parse.Strand` object, we get back a
# `~pydrad.parse.Profile` object which holds information about
# a single snapshot of the simulation.
p = s[0]
print(p)

#################################################################
# The `~pydrad.parse.Profile` object also holds the thermodynamic
# quantities as a function of spatial coordinate at this particular
# timestep. Note that the quantity at each time step is stored as a
# separate 1D array as the spatial grid from one time step to the
# next need not be the same due to the use of the adaptively refined
# grid. Note that each profile is a `~astropy.units.Quantity` and can
# be easily converted to any compatible unit.
print(p.electron_temperature)
print(p.ion_density)
print(p.velocity)

#################################################################
# The `~pydrad.parse.Strand` object also holds an additional
# `~pydrad.parse.Profile` object for the initial conditions,
print(s.initial_conditions)

#################################################################
# The `~pydrad.parse.Strand` object can also be indexed in more complex ways
# in the same manner as a `list` or array to get the simulation only over a
# particular time range or at a particular time cadence.
# As an example, we can grab only the first 5 time steps or every
# other time step.
print(s[:5])
print(s[::2])

#################################################################
# The simplest way to visualize HYDRAD output is to look at a
# single time step as a function of field-aligned coordinate.
# The `~pydrad.parse.Profile` object provides a simple quicklook
# method for this,
p.peek()

#################################################################
# Similarly, we can call this method on a `~pydrad.parse.Strand` to
# overplot every profile at each time step, where the time step is
# mapped to the colormap.
# For additional information about what arguments can be passed in
# here, see the `~pydrad.visualize.plot_profile` documentation.
s.peek()

#################################################################
# Often we want to visualize the whole loop as a function of time
# as a "time-distance" plot. The `~pydrad.parse.Strand` object
# provides a convenience method for this visualization as well.
# Because of the differences in spatial grids between time steps,
# the quantities are reinterpolated onto a uniform grid of specified
# resolution, in this case 0.5 Mm, before plotting.
s.peek_time_distance('electron_temperature', 0.5*u.Mm)

#################################################################
# This quicklook function is flexible enough to visualize multiple
# quantities on different colormaps
s.peek_time_distance(
    ['electron_temperature', 'electron_density', 'velocity'],
    0.5*u.Mm,
    cmap={'velocity': 'RdBu_r'},
    labels={'electron_temperature': r'$T_e$',
            'electron_density': r'$n_e$'},
    norm={'electron_density': ImageNormalize(vmin=1e8, vmax=1e12, stretch=LogStretch())},
    time_unit='h',
    space_unit='Mm',
)
