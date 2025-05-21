"""
Utilities related to parsing HYDRAD results
"""
import pathlib

import astropy.table
import numpy as np
import plasmapy.particles
from pandas import read_csv

__all__ = [
    'read_amr_file',
    'read_phy_file',
    'read_ine_file',
    'read_trm_file',
    'read_hstate_file',
    'read_scl_file',
]


def read_amr_file(filename):
    """
    Parse the adaptive mesh refinement ``.amr`` files containing grid parameters
    and conserved quantities as a function of position.

    .. note:: This is not intended for direct use. Use the `pydrad.parse.Strand`
              object instead.
    """
    # TODO: Account for possible presence of both electron
    # and ion mass densities. If the electron mass density
    # is present in this file, it will shift all of the
    # columns from the index=2 position to the right.
    columns = [
        'grid_centers',
        'grid_widths',
        'mass_density',
        'momentum_density',
        'electron_energy_density',
        'hydrogen_energy_density',
    ]
    units = {
        'grid_centers': 'cm',
        'grid_widths': 'cm',
        'mass_density': 'g cm-3',
        'momentum_density': 'g cm-2 s-1',
        'electron_energy_density': 'erg cm-3',
        'hydrogen_energy_density': 'erg cm-3',
    }
    table = astropy.table.QTable.read(
        filename,
        format='ascii',
        data_start=4,
    )
    # NOTE: This is done after creating the table because the
    # remaining number of columns can be variable and thus we
    # cannot assign all of the column names at once.
    table.rename_columns(
        table.colnames[:len(columns)],
        columns,
    )
    for column in columns:
        table[column].unit = units[column]
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
        -- Loop coordinate (1 column), and at each position:
        -- Terms of mass equation (2 columns)
        -- Terms of momentum equation (6 columns)
        -- Terms of electron energy equation (11 columns)
        -- Terms of hydrogen energy equation (11 columns)
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
    columns = [f'hydrogen_I_level_{i}' for i in range(1,6)] + ['hydrogen_II_fraction']
    return astropy.table.QTable.read(
        filename,
        format='ascii',
        names=columns,
    )

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
