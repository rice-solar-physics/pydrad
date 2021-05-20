"""
Utilities for HYDRAD configuration
"""
from distutils.dir_util import copy_tree
import glob
import os
import platform
import shutil
import subprocess
import tempfile

import astropy.units as u

from pydrad import log

__all__ = ['MissingParameter',
           'HYDRADError',
           'run_shell_command',
           'on_windows',
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


def run_shell_command(cmd, cwd, shell=True):
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
    # Remove "./" from commands if the user is working on Windows
    if on_windows() and cmd[0][0:2] == './':
        cmd[0] = cmd[0][2:]
    cmd = subprocess.run(
        cmd,
        cwd=cwd,
        shell=shell,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout = f"{cmd.stdout.decode('utf-8')}"
    stderr = f"{cmd.stderr.decode('utf-8')}"
    if stdout:
        log.info(stdout)
    if stderr:
        log.warn(stderr)
    hydrad_error_messages = ['segmentation fault', 'abort', 'error']
    if any([e in s.lower() for s in [stderr, stdout] for e in hydrad_error_messages]):
        raise HYDRADError(f'{stderr}\n{stdout}')


def on_windows():
    """
    Determine whether the user's operating system is Windows
    """
    return platform.system().lower() == 'windows'


def get_clean_hydrad(output_path, base_path=None, from_github=False):
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
    """
    # NOTE: this is all done in a temp directory and then copied over
    # so that if something fails, all the files are cleaned up
    with tempfile.TemporaryDirectory() as tmpdir:
        if from_github:
            import git
            git.Repo.clone_from('https://github.com/rice-solar-physics/HYDRAD',
                                tmpdir)
        elif base_path:
            copy_tree(base_path, tmpdir)
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
        # The name of the user guide PDF may change
        user_guide = list(map(os.path.basename, glob.glob(
            os.path.join(tmpdir, 'HYDRAD_User_Guide*.pdf'))))
        rm_files += user_guide
        for d in rm_dirs:
            try:
                shutil.rmtree(os.path.join(tmpdir, d))
            except FileNotFoundError:
                log.warn(f'Cannot remove {d}. Directory not found.')
        for f in rm_files:
            try:
                os.remove(os.path.join(tmpdir, f))
            except FileNotFoundError:
                log.warn(f'Cannot remove {f}. File not found.')
        shutil.copytree(tmpdir, output_path)


def get_equilibrium_heating_rate(root_dir):
    """
    Read equilibrium heating rate from initial conditions results

    Parameters
    ----------
    root_dir : `str` or pathlike
        Path to HYDRAD directory
    """
    filename = os.path.join(root_dir, 'Initial_Conditions/profiles/initial.amr.sol')
    with open(filename, 'r') as f:
        equilibrium_heating_rate = float(f.readline()) * u.Unit('erg cm-3 s-1')
    return equilibrium_heating_rate
