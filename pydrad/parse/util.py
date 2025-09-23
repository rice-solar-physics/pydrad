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
    # NOTE: Using pandas here because it is much faster than using the
    # I/O functionality in astropy for these types of tables
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
    """
    # TODO: clean this up somehow? I've purposefully included
    # a lot of comments because the format of this file makes
    # the parsing code quite opaque
    with pathlib.Path(filename).open() as f:
        lines = f.readlines()
    # First parse all of the population fraction arrays
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
    pop_arrays = [np.zeros((n_s, z+1)) for z in Z]
    for i, p in enumerate(pop_lists):
        for j, line in enumerate(p):
            pop_arrays[j][i, :] = line[1:]  # Skip first entry, it is the atomic number
    columns = []
    for z in Z:
        el_name = plasmapy.particles.element_name(z)
        columns += [f'{el_name}_{i+1}' for i in range(z+1)]
    data = np.hstack([p for p in pop_arrays])
    return astropy.table.QTable(data=data, names=columns)


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
