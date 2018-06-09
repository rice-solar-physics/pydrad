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

    def __init__(self, hydrad_root):
        self.hydrad_root = hydrad_root

    @property
    def time(self):
        """
        Simulation time
        """
        amr_files = glob.glob(os.path.join(self.hydrad_root, 'Results/profile*.amr'))
        time = []
        for af in amr_files:
            with open(af, 'r') as f:
                time.append(float(f.readline()))
        return sorted(time) * u.s

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
            return Profile(os.path.join(self.hydrad_root, 'Results'), index)
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
    results_dir (`str`): Path to HYDRAD results directory
    index (`int`): Timestep index
    """

    def __init__(self, results_dir, index):
        self._index = index
        self.results = np.loadtxt(os.path.join(results_dir, f'profile{index:d}.phy'))

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

    def peek(self, **kwargs):
        """
        Quick look at profiles at at given timestep.
        """
        plot_strand(self, start=self.index, stop=self.index+1, step=1, **kwargs)
