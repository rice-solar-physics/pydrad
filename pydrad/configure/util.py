"""
Utilities for HYDRAD configuration
"""
import subprocess
import platform

from pydrad import log


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
    hydrad_error_messages = ['segmentation fault', 'abort', 'error:']
    if any([e in s.lower() for s in [stderr, stdout] for e in hydrad_error_messages]):
        raise HYDRADError(f'{stderr}\n{stdout}')


def on_windows():
    """
    Determine whether the user's operating system is Windows
    """
    return platform.system().lower() == 'windows'
