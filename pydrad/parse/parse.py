"""
Interface for easily parsing HYDRAD results
"""
import os
import glob

import numpy as np
from scipy.interpolate import splev, splrep
import astropy.units as u

from pydrad.visualize import (plot_strand,
                              plot_profile,
                              animate_strand,
                              plot_time_distance)

__all__ = ['Strand', 'Profile']


def get_master_time(hydrad_root):
    amr_files = glob.glob(os.path.join(hydrad_root, 'Results/profile*.amr'))
    time = []
    for af in amr_files:
        with open(af, 'r') as f:
            time.append(float(f.readline()))
    return sorted(time) * u.s


class Strand(object):
    """
    Container for parsing HYDRAD results

    # Parameters
    hydrad_root (`str`): Path to HYDRAD simulation directory
    """

    def __init__(self, hydrad_root, **kwargs):
        self.hydrad_root = hydrad_root
        # NOTE: time is only specified when slicing a Strand. When not
        # slicing, it should be read from the results.
        self._time = kwargs.pop('time', None)
        if self._time is None:
            self._time = get_master_time(self.hydrad_root)
        # This is different than time depending on the slicing. We allow this
        # to be passed as a kwarg to avoid repeatedly reading multiple files
        # when slicing.
        self._master_time = kwargs.pop('master_time', None)
        if self._master_time is None:
            self._master_time = get_master_time(self.hydrad_root)
        self._profile_kwargs = kwargs

    def __repr__(self):
        return f"""HYDrodynamics and RADiation (HYDRAD) Code
-----------------------------------------
Results path: {self.hydrad_root}
Time interval: [{self.time[0]}, {self.time[-1]}]
Number of profiles: {len(self)}
Loop length: {self.loop_length.to(u.Mm):.3f}"""

    def __len__(self):
        return self.time.shape[0]

    def __getitem__(self, index):
        # NOTE: This will throw an index error to stop iteration
        _ = self.time[index]
        if self.time[index].shape:  # empty if time[index] is a scalar
            return Strand(self.hydrad_root,
                          time=self.time[index],
                          master_time=self._master_time,
                          **self._profile_kwargs)
        else:
            return Profile(self.hydrad_root,
                           self.time[index],
                           master_time=self._master_time,
                           **self._profile_kwargs)

    @property
    def time(self):
        """
        Simulation time
        """
        return self._time

    @property
    def loop_length(self):
        """
        Footpoint-to-footpoint loop length
        """
        with open(os.path.join(self.hydrad_root, 'Results/profile0.amr'), 'r') as f:
            tmp = f.readlines()
            loop_length = float(tmp[2])
        return loop_length * u.cm

    @property
    def initial_conditions(self):
        """
        #Profile for the solutions to the hydrostatic
        equations used as initial conditions.
        """
        return InitialProfile(self.hydrad_root,
                              0*u.s,
                              master_time=self._master_time,
                              **self._profile_kwargs)

    def peek(self, **kwargs):
        """
        Quick look at all profiles for the run on a single plot. Takes
        the same keyword arguments as #pydrad.visualize.plot_strand
        """
        plot_strand(self, **kwargs)

    @u.quantity_input
    def peek_time_distance(self, quantities, delta_s: u.cm, **kwargs):
        """
        Quick look at time-distance plots of various quantities. Takes
        the same keyword arguments as #pydrad.visualize.plot_time_distance
        """
        plot_time_distance(self, quantities, delta_s, **kwargs)

    def animate(self, **kwargs):
        """
        Simple animation of time-dependent loop profiles. Takes the same
        keyword arguments as #pydrad.visualize.animate_strand
        """
        return animate_strand(self, **kwargs)

    @u.quantity_input
    def get_uniform_grid(self, delta_s: u.cm):
        """
        Create a spatial grid with uniform spacing `delta_s`.
        
        # Parameters:
        delta_s (`astropy.units.Quantity`): Spacing between each grid point
        """
        return np.arange(
            0, self.loop_length.to(u.cm).value, delta_s.to(u.cm).value)*u.cm

    def get_unique_grid():
        all_coordinates = [p.coordinate.to(u.cm).value for p in self]
        return np.unique(np.concatenate(all_coordinates).ravel()) * u.cm

    @u.quantity_input
    def to_constant_grid(self, name, grid: u.cm):
        """
        Interpolate a given quantity onto a spatial grid that is the same at
        each time step.

        # Parameters
        name (`str`): Name of quantity
        grid (`astropy.units.Quantity`): Spatial grid to interpolate onto
        """
        q_uniform = np.zeros(self.time.shape+grid.shape)
        grid_cm = grid.to(u.cm).value
        for i, p in enumerate(self):
            q = getattr(p, name)
            tsk = splrep(p.coordinate.to(u.cm).value, q.value,)
            q_uniform[i, :] = splev(grid_cm, tsk, ext=0)

        return q_uniform * q.unit


class Profile(object):
    """
    Container for HYDRAD results at a given timestep. Typically accessed
    through #Strand

    # Parameters
    hydrad_root (`str`): Path to HYDRAD directory
    time (`int`): Timestep index
    """

    @u.quantity_input
    def __init__(self, hydrad_root, time: u.s, **kwargs):
        self.hydrad_root = hydrad_root
        if time.shape:
            raise ValueError('time must be a scalar')
        self.time = time
        self._master_time = kwargs.get('master_time')
        if self._master_time is None:
            self._master_time = get_master_time(self.hydrad_root)
        # Read results files
        self._read_phy()
        if kwargs.get('read_amr', True):
            self._read_amr()

    @property
    def _amr_filename(self):
        return os.path.join(
            self.hydrad_root,
            f'Results/profile{self._index:d}.amr')

    @property
    def _phy_filename(self):
        return os.path.join(
            self.hydrad_root,
            f'Results/profile{self._index:d}.phy')

    @property
    def _index(self):
        return np.where(self.time == self._master_time)[0][0]

    def __repr__(self):
        return f"""HYDRAD Timestep Profile
-----------------------
Filename: {self._phy_filename}
Timestep #: {self._index}"""

    def _read_phy(self):
        """
        Parse the hydrodynamics results file
        """
        self.results = np.loadtxt(self._phy_filename)

    def _read_amr(self):
        """
        Parse the adaptive mesh grid file
        """
        # TODO: Read other quantities from .amr file
        with open(self._amr_filename) as f:
            lines = f.readlines()
            self._grid_centers = np.zeros((int(lines[3]),))
            self._grid_widths = np.zeros((int(lines[3]),))
            for i, l in enumerate(lines[4:]):
                tmp = np.array(l.split(), dtype=float)
                self._grid_centers[i] = tmp[0]
                self._grid_widths[i] = tmp[1]

    @property
    def grid_centers(self):
        """
        Spatial location of the grid centers
        """
        return u.Quantity(self._grid_centers, 'cm')

    @property
    def grid_widths(self):
        """
        Spatial width of each grid cell
        """
        return u.Quantity(self._grid_widths, 'cm')

    @property
    def grid_edges(self):
        """
        Spatial location of left edge of each grid cell
        """
        return self.grid_centers - self.grid_widths/2.

    @property
    def coordinate(self):
        """
        Field-aligned loop coordinate $s$
        """
        return self.results[:, 0] * u.cm

    @property
    def electron_temperature(self):
        """
        Electron temperature $T_e$ as a function of $s$
        """
        return self.results[:, -4] * u.K

    @property
    def ion_temperature(self):
        """
        Ion temperature $T_i$ as a function of $s$
        """
        return self.results[:, -3] * u.K

    @property
    def electron_density(self):
        """
        Electron density $n_e$ as a function of $s$
        """
        return self.results[:, 3] * u.cm**(-3)

    @property
    def ion_density(self):
        """
        Ion density $n_i$ as a function of $s$
        """
        return self.results[:, 4] * u.cm**(-3)

    @property
    def electron_pressure(self):
        """
        Electron pressure $p_e$ as a function of $s$
        """
        return self.results[:, 5] * u.dyne * u.cm**(-2)

    @property
    def ion_pressure(self):
        """
        Ion pressure $p_i$ as a function of $s$
        """
        return self.results[:, 6] * u.dyne * u.cm**(-2)

    @property
    def velocity(self):
        """
        Velocity $v$ as a function of $s$
        """
        return self.results[:, 1] * u.cm / u.s

    def spatial_average(self, quantity, bounds=None):
        """
        Compute a spatial average of a specific quantity
        """
        if bounds is None:
            bounds = self.coordinate[[0, -1]]
        i_bounds = np.where(np.logical_and(self.coordinate >= bounds[0],
                                           self.coordinate <= bounds[1]))
        quantity_bounds = getattr(self, quantity)[i_bounds]
        grid_widths_bounds = self.grid_widths[i_bounds]
        return np.average(quantity_bounds, weights=grid_widths_bounds)

    def peek(self, **kwargs):
        """
        Quick look at profiles at a given timestep.
        """
        plot_profile(self, **kwargs)


class InitialProfile(Profile):

    @property
    def _amr_filename(self):
        return os.path.join(
            self.hydrad_root,
            'Initial_Conditions/profiles/initial.amr')

    @property
    def _phy_filename(self):
        return os.path.join(
            self.hydrad_root,
            'Initial_Conditions/profiles/initial.amr.phy')
