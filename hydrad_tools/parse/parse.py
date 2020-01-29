"""
Interface for easily parsing HYDRAD results
"""
import os
import glob

import numpy as np
from scipy.interpolate import splev, splrep
import astropy.units as u

from hydrad_tools.visualize import plot_strand, animate_strand

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
Number of profiles: {len(self.time)}
Loop length: {self.loop_length.to(u.Mm):.3f}"""

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

    def peek(self, start=0, stop=None, step=100, **kwargs):
        """
        Take a quick look at all profiles for the run on a single plot. Takes
        the same keyword arguments as #hydrad_tools.visualize.plot_strand
        """
        plot_strand(self, start=start, stop=stop, step=step, **kwargs)

    def animate(self, start=0, stop=None, step=100, **kwargs):
        """
        Simple animation of time-dependent loop profiles. Takes the same
        keyword arguments as #hydrad_tools.visualize.animate_strand
        """
        return animate_strand(self, start=start, stop=step, step=step, **kwargs)

    def to_uniform_grid(self, name, delta_s: u.cm):
        """
        Calculate a given quantity on a uniform spatial grid at every time step

        # Parameters
        name (`str`): Name of quantity
        delta_s (`astropy.units.Quantity`): Spatial resolution of uniform grid
        """
        s_uniform = np.arange(0, self.loop_length.to(u.cm).value, delta_s.to(u.cm).value)*u.cm
        q_uniform = np.zeros(self.time.shape+s_uniform.shape)
        # Interpolate each quantity at each timestep
        for i, p in enumerate(self):
            q = getattr(p, name)
            tsk = splrep(p.coordinate.to(u.cm).value, q.value,)
            q_uniform[i, :] = splev(s_uniform.value, tsk, ext=0)

        return s_uniform, q_uniform * q.unit


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
        self._fname = os.path.join(
            hydrad_root, 'Results/profile{index:d}.{ext}')
        # Read results files
        self._read_phy()
        if kwargs.get('read_amr', True):
            self._read_amr()

    @property
    def _index(self):
        return np.where(self.time == self._master_time)[0][0]

    def __repr__(self):
        return f"""HYDRAD Timestep Profile
--------------
Filename: {self._fname.format(index=self._index, ext='phy')}
Timestep #: {self._index}"""

    def _read_phy(self):
        """
        Parse the hydrodynamics results file
        """
        self.results = np.loadtxt(self._fname.format(index=self._index, ext='phy'))

    def _read_amr(self):
        """
        Parse the adaptive mesh grid file
        """
        # TODO: Read other quantities from .amr file
        with open(self._fname.format(index=self._index, ext='amr')) as f:
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
        Quick look at profiles at at given timestep.
        """
        plot_strand(self, start=self._index, stop=self._index+1, step=1, **kwargs)
