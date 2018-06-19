"""
Configure HYDRAD simulations
"""
import os
import copy
import warnings
import datetime
import tempfile
import shutil
import subprocess
from distutils.dir_util import copy_tree

import numpy as np
import astropy.units as u
from jinja2 import Environment, PackageLoader
import asdf
try:
    import git
except ImportError:
    warnings.warn('GitPython not installed. Cannot retrieve base copy from GitHub.')

from . import filters

REMOTE_REPO = 'https://github.com/rice-solar-physics/HYDRAD'

__all__ = ['Configure']


class Configure(object):
    """
    Configure HYDRAD simulations from a single Python `dict`

    # Parameters
    config (`dict`): All input parameters for configuring simulation
    """

    def __init__(self, config, **kwargs):
        self.config = copy.deepcopy(config)
        self.env = Environment(loader=PackageLoader('hydrad_tools', 'configure/templates'))
        self.env.filters['units_filter'] = filters.units_filter
        self.env.filters['log10_filter'] = filters.log10_filter
        self.env.filters['get_atomic_symbol'] = filters.get_atomic_symbol
        self.env.filters['get_atomic_number'] = filters.get_atomic_number
        self.env.filters['sort_elements'] = filters.sort_elements
        self.env.filters['is_required'] = filters.is_required
        # Freeze the date at instantiation; testing purposes only
        if kwargs.get('freeze_date', False):
            self._date = self.date
            self._freeze_date = True

    @staticmethod
    def load_config(filename):
        """
        Load a base configuration from as an ASDF file

        # Parameters
        filename (`str`): Path to ASDF configuration file
        """
        with asdf.open(filename) as af:
            config = af.tree
        return config

    def save_config(self, filename):
        """
        Save the simulation configuration as an ASDF file

        # Parameters
        filename (`str`): Path to YAML configuration file
        """
        asdf.AsdfFile(self.config).write_to(filename)
    
    def setup_simulation(self, output_path, base_path=None, name=None, verbose=True):
        """
        Setup a HYDRAD simulation with desired outputs from a clean copy

        # Parameters
        output_path (`str`):
        base_path (`str`): If None (default), clone a new copy from GitHub 
        (appropriate permissions required)
        name (`str`): Name of the output directory. If None (default), use timestamp
        verbose (`bool`):
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            if base_path is None:
                git.Repo.clone_from(REMOTE_REPO, tmpdir)
            else:
                copy_tree(base_path, tmpdir)
            self.setup_initial_conditions(tmpdir, execute=True, verbose=verbose)
            self.setup_hydrad(tmpdir, verbose=verbose)
            self.save_config(os.path.join(tmpdir, 'hydrad_tools_config.asdf'))
            if name is None:
                name = f'hydrad_{self.date}'
            output_dir = os.path.join(output_path, name)
            shutil.copytree(tmpdir, output_dir)

    def setup_initial_conditions(self, root_dir, execute=True, verbose=True):
        """
        Compile and execute code to get the initial loop profile

        # Parameters
        root_dir (`str`):
        execute (`bool`): If True (default), run the initial conditions code after compiling
        verbose (`bool`):
        """
        files = [
            ('Initial_Conditions/source/config.h', self.initial_conditions_header),
            ('Initial_Conditions/config/initial_conditions.cfg', self.intial_conditions_cfg),
            ('Radiation_Model/source/config.h', self.radiation_header),
            ('Radiation_Model/config/elements_eq.cfg', self.radiation_equilibrium_cfg),
            ('Radiation_Model/config/elements_neq.cfg', self.radiation_nonequilibrium_cfg),
        ]
        if self.config['initial_conditions']['use_tabulated_gravity']:
            self.config['general']['tabulated_gravity_file'] = 'tabulated.gravity'
            files += [('tabulated.gravity', self.tabulated_gravity)]
        for filename, filestring in files:
            with open(os.path.join(root_dir, filename), 'w') as f:
                f.write(filestring)
        cmd = subprocess.run(['chmod', 'u+x', 'build_initial_conditions.bat'],
                             cwd=os.path.join(root_dir, 'Initial_Conditions/build_scripts'),
                             shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,)
        if verbose:
            print(f"{cmd.stdout.decode('utf-8')}\n{cmd.stderr.decode('utf-8')}")
        cmd = subprocess.run(['./build_initial_conditions.bat'],
                             cwd=os.path.join(root_dir, 'Initial_Conditions/build_scripts'),
                             shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,)
        if verbose:
            print(f"{cmd.stdout.decode('utf-8')}\n{cmd.stderr.decode('utf-8')}")
        if not os.path.exists(os.path.join(root_dir, 'Initial_Conditions/profiles')):
            os.mkdir(os.path.join(root_dir, 'Initial_Conditions/profiles'))
        if execute:
            cmd = subprocess.run(['./Initial_Conditions.exe'], cwd=root_dir, shell=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,)
            if verbose:
                print(f"{cmd.stdout.decode('utf-8')}\n{cmd.stderr.decode('utf-8')}")
            if self.config['heating']['background_heating']:
                self.equilibrium_heating_rate = self.get_equilibrium_heating_rate(root_dir)

    def get_equilibrium_heating_rate(self, root_dir):
        """
        Read equilibrium heating rate from initial conditions results

        # Parameters
        root_dir (`str`): Path to HYDRAD directory
        """
        filename = os.path.join(root_dir, 'Initial_Conditions/profiles/initial.amr.sol')
        with open(filename, 'r') as f:
            equilibrium_heating_rate = float(f.readline()) * u.erg / u.s / (u.cm**3)
        return equilibrium_heating_rate

    def setup_hydrad(self, root_dir, verbose=True):
        """
        Compile HYDRAD code with appropriate header and config files

        # Parameters
        root_dir (`str`):
        verbose (`bool`):
        """
        files = [
            ('Heating_Model/source/config.h', self.heating_header),
            ('Heating_Model/config/heating_model.cfg', self.heating_cfg),
            ('HYDRAD/source/config.h', self.hydrad_header),
            ('HYDRAD/source/collisions.h', self.collisions_header),
            ('HYDRAD/config/HYDRAD.cfg', self.hydrad_cfg),
        ]
        if 'tabulated_gravity_profile' in self.config['general']:
            self.config['general']['tabulated_gravity_file'] = 'tabulated.gravity'
            files += [('tabulated.gravity', self.tabulated_gravity)]
        if 'tabulated_cross_section_profile' in self.config['general']:
            self.config['general']['tabulated_cross_section_file'] = 'tabulated.cross_section'
            files += [('tabulated.cross_section', self.tabulated_cross_section)]
        for filename, filestring in files:
            with open(os.path.join(root_dir, filename), 'w') as f:
                f.write(filestring)
        cmd = subprocess.run(['chmod', 'u+x', 'build_HYDRAD.bat'],
                             cwd=os.path.join(root_dir, 'HYDRAD/build_scripts'),
                             shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,)
        if verbose:
            print(f"{cmd.stdout.decode('utf-8')}\n{cmd.stderr.decode('utf-8')}")
        cmd = subprocess.run(['./build_HYDRAD.bat'],
                             cwd=os.path.join(root_dir, 'HYDRAD/build_scripts'),
                             shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,)
        if verbose:
            print(f"{cmd.stdout.decode('utf-8')}\n{cmd.stderr.decode('utf-8')}")
        if not os.path.exists(os.path.join(root_dir, 'Results')):
            os.mkdir(os.path.join(root_dir, 'Results'))
    
    @property
    def date(self):
        """
        Return the current date
        """
        # NOTE: freeze date at instantiation for testing purposes
        if hasattr(self, '_freeze_date') and self._freeze_date:
            return self._date
        else:
            return datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')

    @property
    def intial_conditions_cfg(self):
        """
        Initial conditions configuration file, `Initial_Conditions/config/initial_conditions.cfg`
        """
        return self.env.get_template('initial_conditions.cfg').render(date=self.date, **self.config)

    @property
    def initial_conditions_header(self):
        """
        Initial conditions header file, `Initial_Conditions/source/config.h`
        """
        return self.env.get_template('initial_conditions.config.h').render(
                    date=self.date, **self.config)

    @property
    def hydrad_cfg(self):
        """
        HYDRAD configuration file, `HYDRAD/config/hydrad.cfg`
        """
        return self.env.get_template('hydrad.cfg').render(date=self.date, **self.config)

    @property
    def hydrad_header(self):
        """
        HYDRAD header file, `HYDRAD/source/config.h`
        """
        return self.env.get_template('hydrad.config.h').render(date=self.date, **self.config)

    @property
    def heating_cfg(self):
        """
        Heating model configuration file, `Heating_Model/config/heating.cfg`. If background
        heating is enabled, you must run the initial conditions and set the
        `equilibrium_heating_rate` attribute first.
        """
        if self.config['heating']['background_heating']:
            if not hasattr(self, 'equilibrium_heating_rate'):
                raise AttributeError('No background heating found')
            background_heating_rate = self.equilibrium_heating_rate
        else:
            background_heating_rate = None
        return self.env.get_template('heating.cfg').render(
            date=self.date, background_heating_rate=background_heating_rate, **self.config)

    @property
    def heating_header(self):
        """
        Heating model header file, `Heating_Model/source/config.h`
        """
        return self.env.get_template('heating.config.h').render(date=self.date, **self.config)

    @property
    def radiation_equilibrium_cfg(self):
        """
        Equilibrium elements configuration file, `Radiation_Model/config/elements_eq.cfg`
        """
        elements = self.config['radiation'].get('elements_equilibrium', [])
        return self.env.get_template('radiation.elements.cfg').render(
                    date=self.date, elements=elements, **self.config)

    @property
    def radiation_nonequilibrium_cfg(self):
        """
        Non-equilibrium elements configuration file, `Radiation_Model/config/elements_neq.cfg`
        """
        elements = self.config['radiation'].get('elements_nonequilibrium', [])
        return self.env.get_template('radiation.elements.cfg').render(
                    date=self.date, elements=elements, **self.config)

    @property
    def radiation_header(self):
        """
        Radiation model header file, `Radiation_Model/source/config.h`
        """
        return self.env.get_template('radiation.config.h').render(date=self.date, **self.config)

    @property
    def collisions_header(self):
        """
        Collisions header file, `HYDRAD/source/collisions.h`
        """
        return self.env.get_template('collisions.h').render(date=self.date, **self.config)

    @property
    def tabulated_cross_section(self):
        """
        Sixth-order polynomial fit coefficients for computing flux tube expansion
        """
        return self.env.get_template('coefficients.cfg').render(
            date=self.date,
            coefficients=self.config['general']['tabulated_cross_section_profile'])

    @property
    def tabulated_gravity(self):
        """
        Sixth-order polynomial fit coefficients for computing gravitational acceleration
        """
        return self.env.get_template('coefficients.cfg').render(
            date=self.date,
            coefficients=self.config['general']['tabulated_gravity_profile'])
