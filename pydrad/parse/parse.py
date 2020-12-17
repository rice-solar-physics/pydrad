"""
Interface for easily parsing HYDRAD results
"""
import os
import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep
import astropy.units as u
import plasmapy.particles
import h5py

from pydrad import log
from pydrad.visualize import (plot_strand,
                              plot_profile,
                              animate_strand,
                              plot_time_distance,
                              plot_histogram)

__all__ = ['Strand', 'Profile', 'InitialProfile']


def get_master_time(hydrad_root, read_from_cfg=False):
    amr_files = glob.glob(os.path.join(hydrad_root, 'Results/profile*.amr'))
    if read_from_cfg:
        log.debug('Creating master time array from config files')
        with open(os.path.join(hydrad_root, 'HYDRAD/config/hydrad.cfg'), 'r') as f:
            lines = f.readlines()
        cadence = float(lines[3])
        with open(amr_files[0]) as f:
            start_time = float(f.readline())
        time = start_time + np.arange(len(amr_files)) * cadence
    else:
        log.debug('Reading master time array from all AMR files')
        time = np.zeros((len(amr_files),))
        for i, af in enumerate(amr_files):
            with open(af, 'r') as f:
                time[i] = f.readline()
    return sorted(time) * u.s


class Strand(object):
    """
    Container for parsing HYDRAD results

    # Parameters
    hydrad_root (`str`): Path to HYDRAD simulation directory
    """

    def __init__(self, hydrad_root, **kwargs):
        self.hydrad_root = hydrad_root
        # This is different than time depending on the slicing. We allow this
        # to be passed as a kwarg to avoid repeatedly reading multiple files
        # when slicing.
        self._master_time = kwargs.pop('master_time', None)
        if self._master_time is None:
            self._master_time = get_master_time(self.hydrad_root,
                                                read_from_cfg=kwargs.get('read_from_cfg', False))
        # NOTE: time is only specified when slicing a Strand. When not
        # slicing, it should be read from the results.
        self._time = kwargs.pop('time', self._master_time)
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

    def to_hdf5(self, filename, *variables):
        """
        Save variables to an HDF5 file

        # Parameters
        filename (`str` or path-like): path to HDF file
        variables (`str`): Names of variables to save to file
        """
        with h5py.File(filename, 'w') as hf:
            ds = hf.create_dataset('time', data=self.time.value)
            ds.attrs['unit'] = self.time.unit.to_string()
            for p in self:
                grp = hf.create_group(f'index{p._index}')
                for v in variables:
                    data = getattr(p, v)
                    ds = grp.create_dataset(v, data=data.value)
                    ds.attrs['unit'] = data.unit.to_string()

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
                              master_time=[0]*u.s,
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

    def get_unique_grid(self):
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
    time (`astropy.units.Quantity`): 
    """

    @u.quantity_input
    def __init__(self, hydrad_root, time: u.s, **kwargs):
        self.hydrad_root = hydrad_root
        if time.shape:
            raise ValueError('time must be a scalar')
        self.time = time
        self._master_time = kwargs.get('master_time')
        if self._master_time is None:
            self._master_time = get_master_time(self.hydrad_root,
                                                read_from_cfg=kwargs.get('read_from_cfg', False))
        # Read results files
        self._read_phy()
        self._read_amr()
        if kwargs.get('read_hstate', True):
            self._read_hstate()
        if kwargs.get('read_ine', True):
            self._read_ine()
        if kwargs.get('read_trm', True):
            self._read_trm()

    @property
    def _amr_filename(self):
        return os.path.join(self.hydrad_root,
                            f'Results/profile{self._index:d}.amr')

    @property
    def _phy_filename(self):
        return os.path.join(self.hydrad_root,
                            f'Results/profile{self._index:d}.phy')

    @property
    def _trm_filename(self):
        return os.path.join(self.hydrad_root,
                            f'Results/profile{self._index:d}.trm')
    
    @property
    def _ine_filename(self):
        return os.path.join(self.hydrad_root,
                            f'Results/profile{self._index:d}.ine')
    
    @property
    def _hstate_filename(self):
        return os.path.join(self.hydrad_root,
                            f'Results/profile{self._index:d}.Hstate')

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
        self._phy_data = np.loadtxt(self._phy_filename)

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

    def _read_trm(self):
        """
        Parse the equation terms file
        """
        try:
            with open(self._trm_filename, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            log.debug(f'{self._trm_filename} not found')
            return
        
        n_elements = int(len(lines)/5)
        self._trm_data = np.zeros([n_elements, 3])

        # The files come in sets of 5 rows
        #   -- Loop coordinate, and at each one:
        #   -- Terms of mass equation
        #   -- Terms of momentum equation
        #   -- Terms of electron energy equation
        #   -- Terms of hydrogen energy equation
        # Right now, we only read 3 values from this:
        #  the electron heating, hydrogen heating,
        #  and bolometric radiative losses
        for i in range(len(lines)):
            j = int(i/5)
            line = lines[i].strip().split()
            # Electron heating and radiative loss terms from the 
            # electron energy equation
            if i % 5 == 3:
                self._trm_data[j, 0] = float(line[5])
                self._trm_data[j, 1] = float(line[6])
            # Hydrogen heating term from the hydrogen energy
            # equation
            if i % 5 == 4:
                self._trm_data[j, 2] = float(line[5])
        
        properties = [('electron_heating_term', '_trm_data', 0, 'erg cm-3 s-1'),
                      ('hydrogen_heating_term', '_trm_data', 2, 'erg cm-3 s-1'),
                      ('radiative_loss_term', '_trm_data', 1, 'erg cm-3 s-1'),]
                     
        for p in properties:
            add_property(*p)
        
    def _read_ine(self):
        """
        Parse non-equilibrium ionization population fraction files
        and set attributes for relevant quantities
        """
        # TODO: clean this up somehow? I've purposefully included
        # a lot of comments because the format of this file makes
        # the parsing code quite opaque
        try:
            with open(self._ine_filename, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            log.debug(f'{self._ine_filename} not found')
            return
        # First parse all of the population fraction arrays
        n_s = self.coordinate.shape[0]
        # NOTE: Have to calculate the number of elements we have
        # computed population fractions for as we do not necessarily
        # know this ahead of time
        n_e = int(len(lines)/n_s - 1)
        # The file is arranged in n_s groups of n_e+1 lines each where the first
        # line is the coordinate and the subsequent lines are the population fraction
        # for each element, with each column corresponding to an ion of that element
        # First, separate by coordinate
        pop_lists = [lines[i*(n_e+1)+1:(i+1)*(n_e+1)] for i in range(n_s)]
        # Convert each row of each group into a floating point array
        pop_lists = [[np.array(l.split(), dtype=float) for l in p] for p in pop_lists]
        # NOTE: each row has Z+2 entries as the first entry is the atomic number Z
        # Get these from just the first group as the number of elements is the same
        # for each
        Z = np.array([p[0] for p in pop_lists[0]], dtype=int)
        pop_arrays = [np.zeros((n_s, z+1))for z in Z]
        for i, p in enumerate(pop_lists):
            for j, l in enumerate(p):
                pop_arrays[j][i, :] = l[1:]  # Skip first entry, it is the atomic number

        # Then set attributes for each ion of each element
        for z, p in zip(Z, pop_arrays):
            name = plasmapy.particles.element_name(z)
            attr = f'_population_fraction_{name}'
            setattr(self, attr, p)
            for p in [(f'{attr[1:]}_{i+1}', attr, i, '') for i in range(z+1)]:
                add_property(*p)

    def _read_hstate(self):
        """
        Parse the hydrogen energy level populations file
        """
        try:
            self._hstate_data = np.loadtxt(self._hstate_filename)
        except OSError:
            log.debug(f'{self._hstate_filename} not found')
            self._hstate_data = None

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

    @u.quantity_input
    def column_emission_measure(self, bins:u.K=None):
        """
        Computes the column emission measure, where it is assumed that the loop is
        confined to a single pixel and oriented along the LOS
        
        # Parameters
        bins (`astropy.units.Quantity`): temperature bin edges, including rightmost
        edge. If None (default), the bins will be equally-spaced in $\log{T}$, with
        a left edge at $\log{T}=3$, a right edge at $\log{T}=8$, and a bin width of
        $0.05$.
        
        # Returns
        em (astropy.units.Quantity): the column emission measure in each bin
        
        bins (astropy.units.Quantity): temperature bin edges. Note that `len(bins)=len(em)+1`.
        """
        if bins is None:
            bins = 10.0**(np.arange(3.0, 8.0, 0.05)) * u.K
        weights = self.electron_density * self.ion_density * self.grid_widths
        H, _, _ = np.histogram2d(self.grid_centers, self.electron_temperature,
                                 bins=(self.grid_edges, bins), weights=weights)
        return H.sum(axis=0), bins

    def peek(self, **kwargs):
        """
        Quick look at profiles at a given timestep.
        """
        plot_profile(self, **kwargs)

    def peek_emission_measure(self, **kwargs):
        """
        Quick look at the column emission measure
        """
        bins = kwargs.pop('bins', None)
        if 'color' not in kwargs:
            kwargs['color'] = 'C0'
        em, bins = self.column_emission_measure(bins=bins)
        ax = plot_histogram(em.to('cm-5').value, bins.to('K').value, **kwargs)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(kwargs.get('xlim'))
        ax.set_ylim(kwargs.get('ylim'))
        ax.set_xlabel(r'$T$ [K]')
        ax.set_ylabel(r'EM [cm$^{-5}$]')
        plt.show()


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


def add_property(name, attr, index, unit):
    """
    Auto-generate properties for various pieces of data
    """
    def property_template(self):
        data = getattr(self, attr)
        if data is None:
            raise ValueError(f'No data available for {name}')
        return u.Quantity(data[:, index], unit)
    property_template.__doc__ = f'{" ".join(name.split("_")).capitalize()} as a function of $s$'
    property_template.__name__ = name
    setattr(Profile, property_template.__name__, property(property_template))


properties = [
    ('coordinate', '_phy_data', 0, 'cm'),
    ('velocity', '_phy_data', 1, 'cm / s'),
    ('sound_speed', '_phy_data', 2, 'cm / s'),
    ('electron_density', '_phy_data', 3, 'cm-3'),
    ('ion_density', '_phy_data', 4, 'cm-3'),
    ('electron_pressure', '_phy_data', 5, 'dyne cm-2'),
    ('ion_pressure', '_phy_data', 6, 'dyne cm-2'),
    ('electron_temperature', '_phy_data', 7, 'K'),
    ('ion_temperature', '_phy_data', 8, 'K'),
    ('electron_conduction', '_phy_data', 9, 'erg s-1 cm-2'),
    ('ion_conduction', '_phy_data', 10, 'erg s-1 cm-2'),
]
properties += [(f'level_population_hydrogen_{i}', '_hstate_data', i, '') for i in range(1, 7)]
for p in properties:
    add_property(*p)
