"""
Test HYDRAD simulation setup
"""
import pathlib

import pytest


@pytest.mark.parametrize(
    'filename',
    ['Initial_Conditions.exe',
     'HYDRAD.exe',
     pathlib.Path('Initial_Conditions') / 'profiles' / 'initial.amr',
     pathlib.Path('Initial_Conditions') / 'profiles' / 'initial.amr.gravity',
     pathlib.Path('Initial_Conditions') / 'profiles' / 'initial.amr.phy',
     pathlib.Path('Initial_Conditions') / 'profiles' / 'initial.amr.sol'])
def test_generated_files_exist(hydrad, filename):
    assert (hydrad / filename).is_file()
