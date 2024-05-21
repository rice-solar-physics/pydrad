"""
Interface for easily parsing HYDRAD results
"""
import glob
import os
import pathlib

import astropy.units as u
import h5py
import matplotlib.pyplot as plt
import numpy as np
import plasmapy.particles
from pandas import read_csv
from scipy.interpolate import splev, splrep

from pydrad import log
from pydrad.configure import Configure
from pydrad.visualize import (animate_strand, plot_histogram, plot_profile,
                              plot_strand, plot_time_distance, plot_time_mesh)

__all__ = ['Strand', 'Profile', 'InitialProfile']


def get_master_time(hydrad_root, read_from_cfg=False):
    """
    Get array of times that correspond to each timestep for the entire simulation.

    Parameters
    ----------
    hydrad_root : `str` or path-like
    read_from_cfg : `bool`, optional
        If True, create the time array from the cadence as specified in
        `HYDRAD/config/hydrad.cfg` and the start time as given in the first
        AMR file. Note that this is substantially faster than reading the time
        from every AMR file, but there may be small differences between these
        approximate time steps and the exact time steps listed in the AMR files.

    Returns
    -------
    : `~astropy.units.Quantity`
    """
    hydrad_root = pathlib.Path(hydrad_root)
    amr_files = sorted((hydrad_root / 'Results').glob('profile*.amr'))
    if read_from_cfg:
        log.debug('Creating master time array from config files')
        # NOTE: Sometimes this file is capitalized and some OSes are sensitive to this
        cfg_file = hydrad_root / 'HYDRAD' / 'config' / 'hydrad.cfg'
        if not cfg_file.is_file():
            log.debug('hydrad.cfg not found; trying HYDRAD.cfg')
            cfg_file = hydrad_root / 'HYDRAD' / 'config' / 'HYDRAD.cfg'
        with cfg_file.open() as f:
            lines = f.readlines()
        cadence = float(lines[3])
        with amr_files[0].open() as f:
            start_time = float(f.readline())
        time = start_time + np.arange(len(amr_files)) * cadence
    else:
        log.debug('Reading master time array from all AMR files')
        time = np.zeros((len(amr_files),))
        for i, af in enumerate(amr_files):
            with af.open() as f:
                time[i] = f.readline()
    return sorted(time) * u.s


class Strand(object):
    """
    Container for parsing HYDRAD results

    Parameters
    ----------
    hydrad_root : path-like
        Path to HYDRAD simulation directory
    read_from_cfg : `bool`, optional
        If True, create the time array from the cadence as specified in
        `HYDRAD/config/hydrad.cfg` and the start time as given in the first
        AMR file. Note that this is substantially faster than reading the time
        from every AMR file, but there may be small differences between these
        approximate time steps and the exact time steps listed in the AMR files.
    """

    def __init__(self, hydrad_root, **kwargs):
        self.hydrad_root = pathlib.Path(hydrad_root)
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

        Parameters
        ----------
        filename : `str` or path-like
            Path to HDF file
        variables : `list`
            Names of variables to save to file
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
    def config(self):
        """
        Configuration options. This will only work if the simuation was also
        configured by pydrad.
        """
        return Configure.load_config(self.hydrad_root / 'pydrad_config.asdf')

    @property
    @u.quantity_input
    def time(self) -> u.s:
        """
        Simulation time
        """
        return self._time

    @property
    @u.quantity_input
    def loop_length(self) -> u.cm:
        """
        Footpoint-to-footpoint loop length
        """
        with (self.hydrad_root / 'Results' / 'profile0.amr').open() as f:
            tmp = f.readlines()
            loop_length = float(tmp[2])
        return loop_length * u.cm

    @property
    def initial_conditions(self):
        """
        `~pydrad.parse.Profile` for the solutions to the hydrostatic
        equations used as initial conditions.
        """
        return InitialProfile(self.hydrad_root,
                              0*u.s,
                              master_time=[0]*u.s,
                              **self._profile_kwargs)

    def peek(self, **kwargs):
        """
        Quick look at all profiles for the run on a single plot. Takes
        the same keyword arguments as `~pydrad.visualize.plot_strand`
        """
        plot_strand(self, **kwargs)
        plt.show()

    @u.quantity_input
    def peek_time_distance(self, quantities=None, delta_s: u.cm = 1*u.Mm, **kwargs):
        """
        Quick look at time-distance plots of various quantities. Takes
        the same keyword arguments as `~pydrad.visualize.plot_time_distance`
        """
        if quantities is None:
            quantities = ['electron_temperature', 'electron_density', 'velocity']
        _ = plot_time_distance(self, quantities, delta_s, **kwargs)
        plt.show()

    def animate(self, **kwargs):
        """
        Simple animation of time-dependent loop profiles. Takes the same
        keyword arguments as `~pydrad.visualize.animate_strand`
        """
        return animate_strand(self, **kwargs)

    @u.quantity_input
    def get_uniform_grid(self, delta_s: u.cm) -> u.cm:
        """
        Create a spatial grid with uniform spacing `delta_s`.

        Parameters
        ----------
        delta_s : `astropy.units.Quantity`
            Spacing between each grid point
        """
        return np.arange(
            0, self.loop_length.to(u.cm).value, delta_s.to(u.cm).value)*u.cm

    @u.quantity_input
    def get_unique_grid(self) -> u.cm:
        all_coordinates = [p.coordinate.to(u.cm).value for p in self]
        return np.unique(np.concatenate(all_coordinates).ravel()) * u.cm

    @u.quantity_input
    def to_constant_grid(self, name, grid: u.cm, order=1):
        """
        Interpolate a given quantity onto a spatial grid that is the same at
        each time step.

        Parameters
        ----------
        name : `str`
        grid : `~astropy.units.Quantity`
            Spatial grid to interpolate onto
        order : `int`
            Order of the spline interpolation. Default is 1.
        """
        q_uniform = np.zeros(self.time.shape+grid.shape)
        grid_cm = grid.to(u.cm).value
        for i, p in enumerate(self):
            q = getattr(p, name)
            tsk = splrep(p.coordinate.to(u.cm).value, q.value, k=order)
            q_uniform[i, :] = splev(grid_cm, tsk, ext=0)

        return q_uniform * q.unit

    def spatial_average(self, quantity, bounds=None):
        """
        Compute a spatial average of a specific quantity or quantities
        """
        return u.Quantity([p.spatial_average(quantity, bounds=bounds) for p in self])

    def column_emission_measure(self, bins=None, bounds=None):
        """
        Column emission measure as a function of time

        See Also
        --------
        Profile.column_emission_measure
        """
        _, bins = self[0].column_emission_measure(bins=bins, bounds=bounds)
        em = np.stack([p.column_emission_measure(bins=bins, bounds=bounds)[0] for p in self])
        return em, bins

    def peek_emission_measure(self, bins=None, bounds=None, **kwargs):
        em, bins = self.column_emission_measure(bins=bins, bounds=bounds)
        bin_centers = (bins[1:] + bins[:-1])/2
        # Make API consistent with plot_time_mesh
        for k in ['cmap', 'norm', 'units', 'labels']:
            if k in kwargs:
                kwargs[k] = {'EM': kwargs[k]}
        _ = plot_time_mesh(self, [('EM', em)], bin_centers, r'$T$', yscale='log', **kwargs)
        plt.show()


class Profile(object):
    """
    Container for HYDRAD results at a given timestep. Typically accessed
    through `pydrad.parse.Strand`

    Parameters
    ----------
    hydrad_root : `str`
        Path to HYDRAD directory
    time : `~astropy.units.Quantity`
        Timestep corresponding to the profile of interest
    """

    @u.quantity_input
    def __init__(self, hydrad_root, time: u.s, **kwargs):
        self.hydrad_root = pathlib.Path(hydrad_root)
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
        return self.hydrad_root / 'Results' / f'profile{self._index:d}.amr'

    @property
    def _phy_filename(self):
        return self.hydrad_root / 'Results' / f'profile{self._index:d}.phy'

    @property
    def _trm_filename(self):
        return self.hydrad_root / 'Results' / f'profile{self._index:d}.trm'

    @property
    def _ine_filename(self):
        return self.hydrad_root / 'Results' / f'profile{self._index:d}.ine'

    @property
    def _hstate_filename(self):
        return self.hydrad_root / 'Results' / f'profile{self._index:d}.Hstate'

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
        with self._amr_filename.open() as f:
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

        The files come in sets of 5 rows with variable number of columns:
            -- Loop coordinate (1 column), and at each position:
            -- Terms of mass equation (2 columns)
            -- Terms of momentum equation (6 columns)
            -- Terms of electron energy equation (11 columns)
            -- Terms of hydrogen energy equation (11 columns)
        """
        mass_columns = ['mass_drhobydt','mass_advection']
        momentum_columns = ['momentum_drho_vbydt','momentum_advection','momentum_pressure_gradient',
        'momentum_gravity','momentum_viscous_stress','momentum_numerical_viscosity']
        electron_columns = ['electron_dTEKEbydt','electron_enthalpy','electron_conduction',
        'electron_gravity','electron_collisions','electron_heating','electron_radiative_loss',
        'electron_electric_field','electron_viscous_stress','electron_numerical_viscosity',
        'electron_ionization']
        hydrogen_columns = ['hydrogen_dTEKEbydt','hydrogen_enthalpy','hydrogen_conduction',
        'hydrogen_gravity','hydrogen_collisions','hydrogen_heating','hydrogen_radiative_loss',
        'hydrogen_electric_field','hydrogen_viscous_stress','hydrogen_numerical_viscosity',
        'hydrogen_ionization']

        n_mass = len(mass_columns)
        n_momentum = len(momentum_columns)
        n_electron = len(electron_columns)
        n_hydrogen = len(hydrogen_columns)

        # Terms from the mass equation:
        mass_terms = read_csv(self._trm_filename, sep='\t', header=None, names=mass_columns,
                                skiprows=lambda x: x % 5 != 1 )
        # Terms from the momentum equation:
        momentum_terms = read_csv(self._trm_filename, sep='\t', header=None, names=momentum_columns,
                                skiprows=lambda x: x % 5 != 2 )
        # Terms from the electron energy equation:
        electron_terms = read_csv(self._trm_filename, sep='\t', header=None, names=electron_columns,
                                skiprows=lambda x: x % 5 != 3 )

        # Terms from the hydrogen energy equation:
        hydrogen_terms = read_csv(self._trm_filename, sep='\t', header=None, names=hydrogen_columns,
                                skiprows=lambda x: x % 5 != 4 )

        offsets = [n_mass,
                   n_mass+n_momentum,
                   n_mass+n_momentum+n_electron
                  ]

        n_elements = len(mass_terms)
        self._trm_data = np.zeros([n_elements, offsets[2]+n_hydrogen])

        for i in range(n_elements):
            for j in range(n_mass):
                self._trm_data[i,j] = mass_terms[mass_columns[j]][i]
            for j in range(n_momentum):
                self._trm_data[i,j-offsets[0]] = momentum_terms[momentum_columns[j]][i]
            for j in range(n_electron):
                self._trm_data[i,j-offsets[1]] = electron_terms[electron_columns[j]][i]
            for j in range(n_hydrogen):
                self._trm_data[i,j-offsets[2]] = hydrogen_terms[hydrogen_columns[j]][i]

        properties = []
        for i in range(n_mass):
            properties += [(mass_columns[i], '_trm_data', i, 'g cm-3 s-1')]
        for i in range(n_momentum):
            properties += [(momentum_columns[i], '_trm_data', i+offsets[0], 'dyne cm-3 s-1')]
        for i in range(n_electron):
            properties += [(electron_columns[i], '_trm_data', i+offsets[1], 'erg cm-3 s-1')]
        for i in range(n_hydrogen):
            properties += [(hydrogen_columns[i], '_trm_data', i+offsets[2], 'erg cm-3 s-1')]

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
            with self._ine_filename.open() as f:
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
    @u.quantity_input
    def grid_centers(self) -> u.cm:
        """
        Spatial location of the grid centers
        """
        return u.Quantity(self._grid_centers, 'cm')

    @property
    @u.quantity_input
    def grid_widths(self) -> u.cm:
        """
        Spatial width of each grid cell
        """
        return u.Quantity(self._grid_widths, 'cm')

    @property
    @u.quantity_input
    def grid_edges(self) -> u.cm:
        """
        Spatial location of the edges of each grid cell,
        including the rightmost edge
        """
        left = self.grid_centers - self.grid_widths/2.
        return np.append(left, left[-1:] + self.grid_widths[-1])

    @property
    @u.quantity_input
    def grid_edges_left(self) -> u.cm:
        """
        Spatial location of the left edge of each grid cell
        """
        return self.grid_edges[:-1]

    @property
    @u.quantity_input
    def grid_edges_right(self) -> u.cm:
        """
        Spatial location of the right edge of each grid cell
        """
        return self.grid_edges[1:]

    def spatial_average(self, quantity, bounds=None):
        """
        Compute a spatial average of a specific quantity

        Parameters
        ----------
        quantity : `str`
            Name of the desired quantity to average
        bounds : `~astropy.units.Quantity`, optional
            Array of length 2 specifying the range over which
            to take the spatial average.
        """
        if bounds is None:
            bounds = self.coordinate[[0, -1]]
        i_bounds = np.where(np.logical_and(self.coordinate >= bounds[0],
                                           self.coordinate <= bounds[1]))
        quantity_bounds = getattr(self, quantity)[i_bounds]
        grid_widths_bounds = self.grid_widths[i_bounds]
        return np.average(quantity_bounds, weights=grid_widths_bounds)

    @u.quantity_input
    def column_emission_measure(self, bins: u.K = None, bounds: u.cm = None):
        r"""
        Computes the column emission measure, where it is assumed that the loop
        is confined to a single pixel and oriented along the line of sight.

        Parameters
        ----------
        bins : `~astropy.units.Quantity`, optional
            Temperature bin edges, including rightmost edge. If None (default),
            the bins will be equally-spaced in :math:`\log{T}`, with a left
            edge at :math:`\log{T}=3`, a right edge at :math:`\log{T}=8`, and a
            bin width of :math:`0.05`.

        Returns
        -------
        em : `~astropy.units.Quantity`
            The column emission measure in each bin
        bins : `~astropy.units.Quantity`
            Temperature bin edges. Note that ``len(bins)=len(em)+1``.
        """
        if bins is None:
            bins = 10.0**(np.arange(3.0, 8.0, 0.05)) * u.K
        if bounds is None:
            bounds = self.grid_edges[[0, -1]]
        weights = self.electron_density * self.ion_density * self.grid_widths
        H, _, _ = np.histogram2d(self.grid_centers, self.electron_temperature,
                                 bins=(bounds, bins), weights=weights)
        return H.squeeze(), bins

    def peek(self, **kwargs):
        """
        Quick look at profiles at a given timestep.

        Takes the same keyword arguments as `~pydrad.visualize.plot_profile`.
        """
        plot_profile(self, **kwargs)
        plt.show()

    def peek_emission_measure(self, **kwargs):
        """
        Quick look at the column emission measure.

        Takes the same keyword arguments as `column_emission_measure` and
        `~pydrad.visualize.plot_histogram`.
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
        return self.hydrad_root / 'Initial_Conditions' / 'profiles' / 'initial.amr'

    @property
    def _phy_filename(self):
        return self.hydrad_root / 'Initial_Conditions' / 'profiles' / 'initial.amr.phy'


def add_property(name, attr, index, unit):
    """
    Auto-generate properties for various pieces of data
    """
    def property_template(self):
        data = getattr(self, attr)
        if data is None:
            raise ValueError(f'No data available for {name}')
        return u.Quantity(data[:, index], unit)
    property_template.__doc__ = f'{" ".join(name.split("_")).capitalize()} as a function of :math:`s`'
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
    ('electron_heat_flux', '_phy_data', 9, 'erg s-1 cm-2'),
    ('ion_heat_flux', '_phy_data', 10, 'erg s-1 cm-2'),
]
properties += [(f'level_population_hydrogen_{i}', '_hstate_data', i, '') for i in range(1, 7)]
for p in properties:
    add_property(*p)
