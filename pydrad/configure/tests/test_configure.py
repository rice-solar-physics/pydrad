"""
Test HYDRAD simulation setup
"""
import pathlib
import pytest

from pydrad.configure.data import get_defaults


@pytest.mark.parametrize(
    'filename',
    ['Initial_Conditions.exe',
     'HYDRAD.exe',
     pathlib.Path('Initial_Conditions') / 'profiles' / 'initial.amr',
     pathlib.Path('Initial_Conditions') / 'profiles' / 'initial.amr.gravity',
     pathlib.Path('Initial_Conditions') / 'profiles' / 'initial.amr.phy',
     pathlib.Path('Initial_Conditions') / 'profiles' / 'initial.amr.sol'])
def test_generated_files_exist(hydrad, filename):
    generated_file = pathlib.Path(hydrad / filename)
    assert generated_file.is_file()

def test_default_config():
    default_config = get_defaults()
    assert default_config['radiation']['abundance_dataset'] == 'asplund'
    assert not default_config['initial_conditions']['isothermal']
    assert default_config['general']['write_file_physical']
