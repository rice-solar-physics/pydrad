"""
Some basic smoke tests for HYDRAD runs for different configurations
"""
import astropy.units as u
import platform
import pytest

from pydrad.configure import Configure
from pydrad.configure.util import run_shell_command
from pydrad.parse import Strand


def test_beam_heating_run(tmpdir_factory, configuration, hydrad_clean):
    configuration.config['heating']['beam'] = [
        {
            'time_start': 10*u.s,
            'flux': 1e9 * u.Unit('erg cm-2 s-1'),
            'cut_off': 15 * u.keV,
            'index': 7,
        },
    ]
    configuration.config['general']['total_time'] = 5 * u.s
    hydrad_tmp = tmpdir_factory.mktemp('hydrad_tmp')
    configuration.setup_simulation(hydrad_tmp, hydrad_clean, overwrite=True)
    run_shell_command(hydrad_tmp / 'HYDRAD.exe')
    strand = Strand(hydrad_tmp)
    assert len(strand) == 6


@pytest.mark.skipif('darwin' in platform.system().lower(),
                    reason='OpenMP not included by default on macOS')
def test_use_openmp(tmpdir_factory, configuration, hydrad_clean):
    omp_configuration = Configure(configuration.config)
    omp_configuration.config['general']['use_openmp'] = True
    omp_configuration.config['general']['grid_cells_per_thread'] = 10
    hydrad_tmp = tmpdir_factory.mktemp('hydrad_tmp')
    omp_configuration.setup_simulation(hydrad_tmp, hydrad_clean, overwrite=True)
    run_shell_command(hydrad_tmp / 'HYDRAD.exe')
    strand = Strand(hydrad_tmp)
    assert len(strand) == 6
