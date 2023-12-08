import pytest

from pydrad.configure.util import get_clean_hydrad


def test_clean_hydrad_from_local(tmpdir_factory, hydrad):
    new_hydrad_copy = tmpdir_factory.mktemp('new_hydrad')
    get_clean_hydrad(new_hydrad_copy, base_path=hydrad, from_github=False, overwrite=True)


def test_clean_hydrad_no_args_raises_error():
    with pytest.raises(ValueError, match='Specify local path to HYDRAD or clone from GitHub'):
        get_clean_hydrad('foo', base_path=None, from_github=False)
