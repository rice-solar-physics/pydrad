"""
Utilities related to parsing HYDRAD results
"""
import astropy.table
import astropy.units as u
import numpy as np
import pathlib
import plasmapy.particles

from pandas import read_csv

from pydrad import log

__all__ = [
    'read_master_time',
    'read_amr_file',
    'read_phy_file',
    'read_ine_file',
    'read_trm_file',
    'read_hstate_file',
    'read_scl_file',
]


# Do this here as calling this each time adds significant overhead
# when parsing a file.
ELEMENT_NAME_MAPPING = {
    z: plasmapy.particles.element_name(z) for z in range(1,31)
}


def read_master_time(hydrad_root, read_from_cfg=False):
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


def read_amr_file(filename):
    """
    Parse the adaptive mesh refinement ``.amr`` files containing grid parameters
    and conserved quantities as a function of position.

    .. note:: This is not intended for direct use. Use the `pydrad.parse.Strand`
              object instead.
    """
    columns = [
        'grid_centers',
        'grid_widths',
        'electron_mass_density',
        'mass_density',
        'momentum_density',
        'electron_energy_density',
        'hydrogen_energy_density',
    ]
    units = {
        'grid_centers': 'cm',
        'grid_widths': 'cm',
        'electron_mass_density': 'g cm-3',
        'mass_density': 'g cm-3',
        'momentum_density': 'g cm-2 s-1',
        'electron_energy_density': 'erg cm-3',
        'hydrogen_energy_density': 'erg cm-3',
    }
    # NOTE: Purposefully using pandas explicitly as it seems to be faster
    # than astropy.io.ascii.read for tables with this particular delimiter.
    # I am not completely sure why this is the case but the difference is
    # almost an order of magnitude.
    table = read_csv(
        filename,
        skiprows=4,
        sep=r'\s+',
        header=None,
        engine='c',
    )
    # NOTE: The columns we care about are doubles in HYDRAD, while the
    # other columns are integers with information about the
    # refinement level of the grid cell. As a result, if electron
    # mass density is not present in the .amr file, then the
    # seventh column is an integer.
    if table.dtypes[len(columns)-1] == np.int64:
        columns.remove('electron_mass_density')
        del units['electron_mass_density']
    # NOTE: This is done after creating the table because the
    # remaining number of columns can be variable and thus we
    # cannot assign all of the column names at once.
    table = table.truncate(
        after=len(columns)-1,
        axis='columns'
    )
    table.rename(columns={i:name for i, name in enumerate(columns)}, inplace=True)
    table = astropy.table.QTable.from_pandas(table, units=units)
    return table


def read_phy_file(filename):
    """
    Parse physical variables ``.phy`` files containing thermodynamic
    quantities as a function of position.

    .. note:: This is not intended for direct use. Use the `pydrad.parse.Strand`
              object instead.
    """
    columns = [
        'coordinate',
        'velocity',
        'sound_speed',
        'electron_density',
        'hydrogen_density',
        'electron_pressure',
        'hydrogen_pressure',
        'electron_temperature',
        'hydrogen_temperature',
        'electron_heat_flux',
        'hydrogen_heat_flux',
    ]
    units = {
        'coordinate': 'cm',
        'velocity': 'cm / s',
        'sound_speed': 'cm / s',
        'electron_density': 'cm-3',
        'hydrogen_density': 'cm-3',
        'electron_pressure': 'dyne cm-2',
        'hydrogen_pressure': 'dyne cm-2',
        'electron_temperature': 'K',
        'hydrogen_temperature': 'K',
        'electron_heat_flux': 'erg s-1 cm-2',
        'hydrogen_heat_flux': 'erg s-1 cm-2',
    }
    return astropy.table.QTable.read(
        filename,
        format='ascii',
        names=columns,
        units=units,
    )


def read_ine_file(filename, n_s):
    """
    Parse the ionization non-equilibrium ``.ine`` files containing
    ionization fraction as a function of position.

    .. note:: This is not intended for direct use. Use the `pydrad.parse.Strand`
              object instead.

    Parameters
    ----------
    filename: path-like
    n_s: `int`
        The number of grid cells in the snapshot corresponding to this file.
    """
    # This file is grouped into n_s groups each of length n_el + 1 (because the first entry is
    # the spatial coordinate) such that the total number of lines is n_s*(n_el + 1).
    # Each line in the group (except the first line) has Z+2 entries corresponding to Z followed
    # by the ionization fraction of the Z+1 ionization stages of element Z at the spatial
    # coordinate specified in the first line of the group.
    # Because of the complexity of the structure of this file, we need to parse it line by line.
    with filename.open(mode='r') as f:
        lines = f.readlines()
    n_el = int(len(lines)/n_s - 1)
    # The innermost loop parses the ionization fraction for all ionization stages of a given element Z
    # at all spatial coordinates and casts it to an array. This innermost array has dimensions (n_s,Z+1).
    # The outermost array iterates over all elements. The result is a list of length n_el where each entry
    # contains the ionization fractions at all ionization stages of a given element at all spatial coordinates.
    data = [
        np.asarray(
            [lines[(1+n_el)*i_s+1+i_z].split()[1:] for i_s in range(n_s)],
            dtype=np.float64
        )
        for i_z in range(n_el)
    ]
    Z = [x.shape[1]-1 for x in data]
    colnames = []
    for z in Z:
        # A precomputed mapping between Z and element name is used as calling plasmapy.particles.element_name
        # each time leads to significant overhead.
        colnames += [f'{ELEMENT_NAME_MAPPING[z]}_{i}' for i in range(1, z+2)]
    return astropy.table.Table(data=np.hstack(data), names=colnames, copy=False)


def read_trm_file(filename):
    """
    Parse ``.trm`` files with hydrodynamic equation terms as a function of position.

    The files come in sets of 5 rows with variable number of columns:

    * Loop coordinate (1 column), and at each position:
    * Terms of mass equation (2 columns)
    * Terms of momentum equation (6 columns)
    * Terms of electron energy equation (11 columns)
    * Terms of hydrogen energy equation (11 columns)
    """
    units = {
        'mass': 'g cm-3 s-1',
        'momentum': 'dyne cm-3 s-1',
        'electron': 'erg cm-3 s-1',
        'hydrogen': 'erg cm-3 s-1',
    }
    columns = {
        'mass': [
            'mass_drhobydt',
            'mass_advection',
        ],
        'momentum': [
            'momentum_drho_vbydt',
            'momentum_advection',
            'momentum_pressure_gradient',
            'momentum_gravity',
            'momentum_viscous_stress',
            'momentum_numerical_viscosity',
        ],
        'electron': [
            'electron_dTEKEbydt',
            'electron_enthalpy',
            'electron_conduction',
            'electron_gravity',
            'electron_collisions',
            'electron_heating',
            'electron_radiative_loss',
            'electron_electric_field',
            'electron_viscous_stress',
            'electron_numerical_viscosity',
            'electron_ionization'
        ],
        'hydrogen': [
            'hydrogen_dTEKEbydt',
            'hydrogen_enthalpy',
            'hydrogen_conduction',
            'hydrogen_gravity',
            'hydrogen_collisions',
            'hydrogen_heating',
            'hydrogen_radiative_loss',
            'hydrogen_electric_field',
            'hydrogen_viscous_stress',
            'hydrogen_numerical_viscosity',
            'hydrogen_ionization',
        ]
    }
    tables = []
    for i,(k,v) in enumerate(columns.items()):
        tables.append(
            astropy.table.QTable.from_pandas(
                read_csv(
                    filename,
                    sep='\t',
                    header=None,
                    names=v,
                    skiprows=lambda x: x % 5 != (i+1)
                ),
                units={c: units[k] for c in v}
            )
        )
    return astropy.table.hstack(tables)


def read_hstate_file(filename):
    """
    Parse the ``.hstate`` files containing the hydrogen level populations as a function of position
    """
    columns = ['coordinate'] + [f'hydrogen_I_level_{i}' for i in range(1,6)] + ['hydrogen_II_fraction']
    table = astropy.table.QTable.read(
        filename,
        format='ascii.no_header',
        delimiter='\t',
        names=columns,
    )
    # The coordinate is already stored from the .amr file, so we don't need to save it.
    # However, we need to parse the correct number of columns.
    table.remove_column('coordinate')
    return table


def read_scl_file(filename):
    """
    Parse the ``.scl`` files containing the time-scales as a function of position
    """
    columns = [
        'coordinate',
        'grid_widths',
        'advective_timescale',
        'electron_conductive_timescale',
        'ion_conductive_timescale',
        'viscous_timescale',
        'collisional_timescale',
        'radiative_timescale',
        'free_bound_timescale',
    ]
    units = {
        'coordinate': 'cm',
        'grid_widths': 'cm',
        'advective_timescale': 's',
        'electron_conductive_timescale': 's',
        'ion_conductive_timescale': 's',
        'viscous_timescale': 's',
        'collisional_timescale': 's',
        'radiative_timescale': 's',
        'free_bound_timescale': 's',
    }
    table = astropy.table.QTable.read(
        filename,
        format='ascii',
    )
    # NOTE: This is done after creating the table because there
    # is an extra column (the free-bound or bound-free timescale)
    # when non-equilibrium ionization is switched on.
    n_columns = len(table.colnames)
    table.rename_columns(
        table.colnames[:n_columns],
        columns[:n_columns],
    )
    for column in columns[:n_columns]:
        table[column].unit = units[column]
    return table
