"""
Utilities for HYDRAD configuration
"""
import os
import pathlib
import platform
import shutil
import stat
import subprocess
import tempfile

import astropy.units as u

from pydrad import log

__all__ = ['MissingParameter',
           'HYDRADError',
           'run_shell_command',
           'get_clean_hydrad',
           'get_equilibrium_heating_rate']


class MissingParameter(Exception):
    """
    An error to raise if a parameter is missing in the configuration
    """
    pass


class HYDRADError(Exception):
    """
    An error to raise if something has gone wrong when running the HYDRAD code
    """
    pass


def run_shell_command(path, shell=True, **kwargs):
    """
    Wrapper function for running shell commands.

    This is essentially a light wrapper around `subprocess.run`
    and includes some error handling and logging specific to
    running the shell commands needed to run and compile HYDRAD.

    Parameters
    ----------
    cmd : `list`
        Command to run
    cwd : `str` or path-like
        Directory to run the command in
    shell : `bool`, optional
        See `~subprocess.run`

    Raises
    ------
    HYDRADError
        This error is raised if HYDRAD or initial conditions code
        fails to compile or if there is a runtime error in the
        initial conditions code.
    """
    path = pathlib.Path(path)
    on_windows = platform.system().lower() == 'windows'
    cmd = subprocess.run(
        path.name if on_windows else f'./{path.name}',
        cwd=path.parent,
        shell=shell,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=os.environ,
        **kwargs,
    )
    stdout = f"{cmd.stdout.decode('utf-8')}"
    stderr = f"{cmd.stderr.decode('utf-8')}"
    if stdout:
        log.info(stdout)
    if stderr:
        log.warning(stderr)
    hydrad_error_messages = ['segmentation fault', 'abort', 'error', 'trace trap']
    if any([e in s.lower() for s in [stderr, stdout] for e in hydrad_error_messages]):
        raise HYDRADError(f'{stderr}\n{stdout}')


def get_clean_hydrad(output_path, base_path=None, from_github=False, overwrite=False):
    """
    Create a clean copy of HYDRAD with only the files necessary to run the code.
    May be useful when making many copies.

    Parameters
    ----------
    output_path : pathlike
        Path to the new stripped down version
    base_path : pathlike, optional
        Path to the original copy. This will not be modified.
    from_github : `bool`, optional
        If True, grab the latest copy of HYDRAD from GitHub. In this case,
        `base_path` is ignored. Note that this requires the GitPython package.
    overwrite : `bool`, optional
        If True, overwrite the directory at `output_path` if it exists. You may need
        to set this to true of the path you are writing your clean copy to HYDRAD to
        already exists, but is empty.
    """
    # NOTE: This function is needed for handling permissions errors that occur on Windows
    # when trying to remove files. See this GH issue comment for more information:
    # https://github.com/python/cpython/issues/87823#issuecomment-1093908280
    def _handle_readonly(func, path, exc_info):
        if func not in (os.unlink, os.rmdir) or exc_info[1].winerror != 5:
            raise exc_info[1]
        os.chmod(path, stat.S_IWRITE)
        func(path)

    # NOTE: this is all done in a temp directory and then copied over
    # so that if something fails, all the files are cleaned up
    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)
        if from_github:
            import git
            git.Repo.clone_from('https://github.com/rice-solar-physics/HYDRAD',
                                tmpdir)
        elif base_path:
            shutil.copytree(base_path, tmpdir, dirs_exist_ok=True)
        else:
            raise ValueError('Specify local path to HYDRAD or clone from GitHub')
        rm_dirs = [
            'Forward_Model',
            'HYDRAD_GUI',
            'Visualisation',
            '.git',
        ]
        rm_files = [
            'HYDRAD_GUI.jar',
            'LICENSE',
            'README.md',
            '.gitignore',
            'hydrad-logo.png',
        ]
        # NOTE: Use glob because name of user guide PDF may change
        rm_files += [f.name for f in tmpdir.glob('HYDRAD_User_Guide*.pdf')]
        for d in rm_dirs:
            try:
                shutil.rmtree(pathlib.Path(tmpdir) / d, onerror=_handle_readonly)
            except FileNotFoundError:
                log.warning(f'Cannot remove {d}. Directory not found.')
        for f in rm_files:
            (tmpdir / f).unlink(missing_ok=True)
        shutil.copytree(tmpdir, output_path, dirs_exist_ok=overwrite)


def get_equilibrium_heating_rate(root_dir):
    """
    Read equilibrium heating rate from initial conditions results

    Parameters
    ----------
    root_dir : `str` or pathlike
        Path to HYDRAD directory
    """
    filename = pathlib.Path(root_dir) / 'Initial_Conditions' / 'profiles' / 'initial.amr.sol'
    with filename.open() as f:
        equilibrium_heating_rate = float(f.readline()) * u.Unit('erg cm-3 s-1')
    return equilibrium_heating_rate
