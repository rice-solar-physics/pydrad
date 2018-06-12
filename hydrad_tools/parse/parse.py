"""
Interface for easily parsing HYDRAD results
"""
import os
import glob

import numpy as np
import astropy.units as u

from hydrad_tools.visualize import plot_strand, animate_strand

__all__ = ['Strand', 'Profile']


class Strand(object):
    """
    Container for parsing HYDRAD results

    # Parameters
    hydrad_root (`str`): Path to HYDRAD simulation directory
    """

    def __init__(self, hydrad_root, **kwargs):
        self.hydrad_root = hydrad_root
        self._profile_kwargs = kwargs
        self._time = self._read_time()

    def _read_time(self):
        amr_files = glob.glob(os.path.join(self.hydrad_root, 'Results/profile*.amr'))
        time = []
        for af in amr_files:
            with open(af, 'r') as f:
                time.append(float(f.readline()))
        return sorted(time) * u.s

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
        if index < self.time.shape[0]:
            return Profile(self.hydrad_root, index, **self._profile_kwargs)
        else:
            raise IndexError

    def peek(self, start=0, stop=None, step=100, **kwargs):
        """
        Take a quick look at all profiles for the run on a single plot. Takes
        the same keyword arguments as #hydrad_tools.visualize.plot_strand
        """
        plot_strand(self, start=start, stop=stop, step=step, **kwargs)

    def animate(self, start=0, stop=None, step=100, **kwargs):
        """
        Simple animation of time-dependent loop profiles. Takes the same keyword
        arguments as #hydrad_tools.visualize.animate_strand
        """
        return animate_strand(self, start=start, stop=step, step=step, **kwargs)


class Profile(object):
    """
    Container for HYDRAD results at a given timestep. Typically accessed through #Strand

    # Parameters
    hydrad_root (`str`): Path to HYDRAD directory
    index (`int`): Timestep index
    """

    def __init__(self, hydrad_root, index, **kwargs):
        self.hydrad_root = hydrad_root
        self._index = index
        self._fname = os.path.join(hydrad_root, 'Results/profile{index:d}.{ext}')
        # Read results files
        self._read_phy()
        if kwargs.get('read_amr', True):
            self._read_amr()

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
