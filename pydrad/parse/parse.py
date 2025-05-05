"""
Interface for easily parsing HYDRAD results
"""
import pathlib

import astropy.constants as const
import astropy.units as u
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splev, splrep

import pydrad.util.constants
from pydrad import log
from pydrad.configure import Configure
from pydrad.parse.util import (read_amr_file, read_hstate_file, read_ine_file,
                               read_phy_file, read_trm_file)
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
        profile_kwargs = self._profile_kwargs.copy()
        # NOTE: Ensures that all results are only derived from the AMR file
        # as the quantities in the phy file for the initial conditions can potentially
        # have a different shape compared to those from the AMR file.
        profile_kwargs['read_phy'] = False
        # NOTE: These files do not exist for the initial conditions
        profile_kwargs['read_hstate'] = False
        profile_kwargs['read_ine'] = False
        profile_kwargs['read_trm'] = False
        return InitialProfile(self.hydrad_root,
                              0*u.s,
                              master_time=[0]*u.s,
                              **profile_kwargs)

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
    hydrad_root : path-like
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
        self._read_amr()
        if kwargs.get('read_phy', True):
            self._read_phy()
        if kwargs.get('read_ine', True):
            self._read_ine()
        if kwargs.get('read_hstate', True):
            self._read_hstate()
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

    def _read_amr(self):
        """
        Parse the adaptive mesh refinement (``.amr``) file
        """
        self._amr_data = read_amr_file(self._amr_filename)

    def _read_phy(self):
        """
        Parse the physical variables (``.phy``) file and set attributes
        for relevant quantities.
        """
        if not self._phy_filename.is_file():
            log.warning(f'{self._phy_filename} not found. Skipping parsing of .phy files. Set read_phy=False to suppress this warning.')
            return
        self._phy_data = read_phy_file(self._phy_filename)
        # NOTE: only adding three columns in this manner as the remaining columns are dealt with explicitly
        # below because of the overlap with the derived quantities from the .amr files.
        for column in ['sound_speed', 'electron_heat_flux', 'hydrogen_heat_flux']:
            _add_property(column, '_phy_data')

    def _read_trm(self):
        """
        Parse the equation terms (``.trm``) file and set attributes
        for relevant quantities.
        """
        if not self._trm_filename.is_file():
            log.warning(f'{self._trm_filename} not found. Skipping parsing of .trm files. Set read_trm=False to suppress this warning.')
            return
        self._trm_data = read_trm_file(self._trm_filename)
        for col in self._trm_data.colnames:
            _add_property(col, '_trm_data')

    def _read_ine(self):
        """
        Parse non-equilibrium ionization population fraction (``.ine``) file
        and set attributes for relevant quantities
        """
        if not self._ine_filename.is_file():
            log.warning(f'{self._ine_filename} not found. Skipping parsing of .ine files. Set read_ine=False to suppress this warning.')
            return
        self._ine_data = read_ine_file(self._ine_filename, self.coordinate.shape[0])
        for col in self._ine_data.colnames:
            _add_property(col, '_ine_data')

    def _read_hstate(self):
        """
        Parse the hydrogen energy level populations (``.hstate``) file and set
        the relevant attributes.
        """
        if not self._hstate_filename.is_file():
            log.warning(f'{self._hstate_filename} not found. Skipping parsing of .hstate files. Set read_hstate=False to suppress this warning.')
            return
        self._hstate_data = read_hstate_file(self._hstate_filename)
        for col in self._hstate_data.colnames:
            _add_property(col, '_hstate_data')

    @property
    @u.quantity_input
    def grid_centers(self) -> u.cm:
        """
        Spatial location of the grid centers
        """
        return self._amr_data['grid_centers']

    @property
    @u.quantity_input
    def grid_widths(self) -> u.cm:
        """
        Spatial width of each grid cell
        """
        return self._amr_data['grid_widths']

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

    @property
    @u.quantity_input
    def coordinate(self) -> u.cm:
        """
        Spatial coordinate, :math:`s`, in the field-aligned direction.

        An alias for `~pydrad.parse.Profile.grid_centers`.
        """
        return self.grid_centers

    @property
    @u.quantity_input
    def electron_mass_density(self) -> u.Unit('g cm-3'):
        # TODO: Account for possible presence of both electron
        # and ion mass densities. If the electron mass density
        # is present in this file, it will shift all of the
        # indices in the .amr file by one.
        if 'electron_mass_density' in self._amr_data.colnames:
            return self._amr_data['electron_mass_density']
        return self.mass_density

    @property
    @u.quantity_input
    def mass_density(self) -> u.Unit('g cm-3'):
        r"""
        Mass density, :math:`\rho`, as a function of :math:`s`.

        .. note:: This is a conserved quantity in HYDRAD.
        """
        return self._amr_data['mass_density']

    @property
    @u.quantity_input
    def momentum_density(self) -> u.Unit('g cm-2 s-1'):
        r"""
        Momentum density, :math:`\rho v`, as a function of :math:`s`.

        .. note:: This is a conserved quantity in HYDRAD.
        """
        return self._amr_data['momentum_density']

    @property
    @u.quantity_input
    def electron_energy_density(self) -> u.Unit('erg cm-3'):
        r"""
        Electron energy density, :math:`E_e`, as a function of :math:`s`.

        .. note:: This is a conserved quantity in HYDRAD.
        """
        return self._amr_data['electron_energy_density']

    @property
    @u.quantity_input
    def hydrogen_energy_density(self) -> u.Unit('erg cm-3'):
        r"""
        Hydrogen energy density, :math:`E_H`, as a function of :math:`s`.

        .. note:: This is a conserved quantity in HYDRAD.
        """
        return self._amr_data['hydrogen_energy_density']

    @property
    @u.quantity_input
    def velocity(self) -> u.cm/u.s:
        r"""
        Bulk velocity, :math:`v`, as a function of :math:`s`.
        """
        if hasattr(self, '_phy_data'):
            return self._phy_data['velocity']
        return self.momentum_density / self.mass_density

    @property
    @u.quantity_input
    def hydrogen_density(self) -> u.Unit('cm-3'):
        r"""
        Hydrogen density, :math:`n_H=\rho/\bar{m_i}`, as a function of :math:`s`,
        where :math:`\bar{m_i} is the average ion mass of a H-He plasma.
        """
        if hasattr(self, '_phy_data'):
            return self._phy_data['hydrogen_density']
        return self.mass_density / pydrad.util.constants.m_avg_ion

    @property
    @u.quantity_input
    def electron_density(self) -> u.Unit('cm-3'):
        r"""
        Electron density, :math:`n_e`, as a function of :math:`s`.
        In nearly all cases, this is equal to `~pydrad.parse.Profile.hydrogen_density`.
        """
        if hasattr(self, '_phy_data'):
            return self._phy_data['electron_density']
        # FIXME: If this exists as a separate column in the .amr file then
        # choose that. Otherwise, the electron and ion densities are assumed
        # to be the same.
        return self.hydrogen_density

    @property
    @u.quantity_input
    def electron_pressure(self) -> u.Unit('dyne cm-2'):
        r"""
        Electron pressure, :math:`P_e=(\gamma - 1)E_e`, as a function of :math:`s`.
        """
        if hasattr(self, '_phy_data'):
            return self._phy_data['electron_pressure']
        return self.electron_energy_density / (pydrad.util.constants.gamma - 1)

    @property
    @u.quantity_input
    def hydrogen_pressure(self) -> u.Unit('dyne cm-2'):
        r"""
        Hydrogen pressure, :math:`P_H = (\gamma - 1)(E_H - \rho v^2/2)`, as a function of :math:`s`.
        """
        if hasattr(self, '_phy_data'):
            return self._phy_data['hydrogen_pressure']
        return (
            self.hydrogen_energy_density
            - self.momentum_density**2/(2*self.mass_density)
        )*(pydrad.util.constants.gamma - 1)

    @property
    @u.quantity_input
    def total_pressure(self) -> u.Unit('dyne cm-2'):
        r"""
        Total pressure, :math:`P = P_e + P_H`, as a function of :math:`s`.
        """
        return self.electron_pressure + self.hydrogen_pressure

    @property
    @u.quantity_input
    def electron_temperature(self) -> u.Unit('K'):
        r"""
        Electron temperature, :math:`T_e = P_e / (k_B n_e)`, as a function of :math:`s`.
        """
        if hasattr(self, '_phy_data'):
            return self._phy_data['electron_temperature']
        return self.electron_pressure / (const.k_B*self.electron_density)

    @property
    @u.quantity_input
    def hydrogen_temperature(self) -> u.Unit('K'):
        r"""
        Hydrogen temperature, :math:`T_H = P_H / (k_B n_H)`, as a function of :math:`s`.
        """
        if hasattr(self, '_phy_data'):
            return self._phy_data['hydrogen_temperature']
        return self.hydrogen_pressure / (const.k_B*self.hydrogen_density)

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
        weights = self.electron_density * self.hydrogen_density * self.grid_widths
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


def _add_property(name, attr,):
    """
    Auto-generate properties for various pieces of data
    """
    def property_template(self):
        data = getattr(self, attr)
        if data is None:
            raise ValueError(f'No data available for {name}')
        return u.Quantity(data[name])
    property_template.__doc__ = f'{" ".join(name.split("_")).capitalize()} as a function of :math:`s`'
    property_template.__name__ = name
    setattr(Profile, property_template.__name__, property(property_template))
