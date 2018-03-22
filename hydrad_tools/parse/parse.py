"""
Interface for easily parsing HYDRAD results
"""
import os
import glob

import numpy as np
import astropy.units as u

__all__ = ['Strand', 'Profile']


class Strand(object):

    def __init__(self, hydrad_root):
        self.hydrad_root = hydrad_root

    @property
    def time(self):
        amr_files = glob.glob(os.path.join(self.hydrad_root, 'Results/profile*.amr'))
        time = []
        for af in amr_files:
            with open(af, 'r') as f:
                time.append(float(f.readline()))
        return sorted(time) * u.s

    @property
    def loop_length(self):
        with open(os.path.join(self.hydrad_root, 'Results/profile0.amr'), 'r') as f:
            tmp = f.readlines()
            loop_length = float(tmp[2])
        return loop_length * u.cm

    def __getitem__(self, index):
        if index < self.time.shape[0]:
            return Profile(os.path.join(self.hydrad_root, 'Results'), index)
        else:
            raise IndexError


class Profile(object):

    def __init__(self, results_dir, index):
        self.results = np.loadtxt(os.path.join(results_dir, f'profile{index:d}.phy'))

    @property
    def coordinate(self):
        return self.results[:, 0] * u.cm

    @property
    def electron_temperature(self):
        return self.results[:, -4] * u.K

    @property
    def ion_temperature(self):
        return self.results[:, -3] * u.K

    @property
    def electron_density(self):
        return self.results[:, 3] * u.cm**(-3)

    @property
    def ion_density(self):
        return self.results[:, 4] * u.cm**(-3)

    @property
    def electron_pressure(self):
        return self.results[:, 5] * u.dyne * u.cm**(-2)

    @property
    def ion_pressure(self):
        return self.results[:, 6] * u.dyne * u.cm**(-2)

    @property
    def velocity(self):
        return self.results[:, 1] * u.cm / u.s
