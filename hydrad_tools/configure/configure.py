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
from jinja2 import Environment, PackageLoader, ChoiceLoader, DictLoader
import asdf
try:
    import git
except ImportError:
    warnings.warn('GitPython not installed. Cannot clone HYDRAD from GitHub.')

from . import filters

REMOTE_REPO = 'https://github.com/rice-solar-physics/HYDRAD'

__all__ = ['Configure']


class Configure(object):
    """
    Configure HYDRAD simulations from a single Python `dict`

    # Parameters
    config (`dict`): All input parameters for configuring simulation
    templates (`dict`): Templates to override defaults, optional
    """

    def __init__(self, config, **kwargs):
        self.config = copy.deepcopy(config)
        loader = ChoiceLoader([
            DictLoader(kwargs.get('templates', {})),
            PackageLoader('hydrad_tools', 'configure/templates')
        ])
        self.env = Environment(loader=loader)
        self.env.filters['units_filter'] = filters.units_filter
        self.env.filters['log10_filter'] = filters.log10_filter
        self.env.filters['get_atomic_symbol'] = filters.get_atomic_symbol
        self.env.filters['get_atomic_number'] = filters.get_atomic_number
        self.env.filters['sort_elements'] = filters.sort_elements
        self.env.filters['is_required'] = filters.is_required
        # NOTE: Freeze the date at instantiation so that files can be compared
        # exactly for testing
        if kwargs.get('freeze_date', False):
            self._date = self.date
            self._freeze_date = True

    def _run_shell_command(self, cmd, cwd, shell=True, verbose=False):
        cmd = subprocess.run(
            cmd,
            cwd=cwd,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if verbose:
            print(f"{cmd.stdout.decode('utf-8')}\n{cmd.stderr.decode('utf-8')}")

    @staticmethod
    def load_config(filename):
        """
        Load a base configuration from an ASDF file

        # Parameters
        filename (`str`): Path to ASDF configuration file
        """
        with asdf.open(filename) as af:
            config = copy.deepcopy(dict(af.tree))
        return config

    def save_config(self, filename):
        """
        Save the simulation configuration as an ASDF file

        # Parameters
        filename (`str`): Path to YAML configuration file
        """
        asdf.AsdfFile(self.config).write_to(filename)

    def setup_simulation(self, output_path, base_path=None, name=None,
                         verbose=True, run_initial_conditions=True):
        """
        Setup a HYDRAD simulation with desired outputs from a clean copy

        # Parameters
        output_path (`str`): Directory in which to put new copy of HYDRAD
        base_path (`str`): Path to existing HYDRAD. If None (default), clone a
        new copy from GitHub (appropriate permissions required)
        name (`str`): Name of the output directory. If None (default), use
        timestamp
        verbose (`bool`):
        run_initial_conditions (`bool`): If True, compile and run the initial
        conditions code
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            if base_path is None:
                git.Repo.clone_from(REMOTE_REPO, tmpdir)
            else:
                copy_tree(base_path, tmpdir)
            if run_initial_conditions:
                self.setup_initial_conditions(tmpdir, execute=True,
                                              verbose=verbose)
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
        execute (`bool`): If True (default), compute initial conditions
        verbose (`bool`):
        """
        files = [
            ('Initial_Conditions/source/config.h',
             self.initial_conditions_header),
            ('Initial_Conditions/config/initial_conditions.cfg',
             self.initial_conditions_cfg),
            ('Radiation_Model/source/config.h',
             self.radiation_header),
            ('Radiation_Model/config/elements_eq.cfg',
             self.radiation_equilibrium_cfg),
            ('Radiation_Model/config/elements_neq.cfg',
             self.radiation_nonequilibrium_cfg),
        ]
        # NOTE: there are two options here so that the gravitational and
        # magnetic field polynomial fits can be applied just to the
        # hydrodynamic step and not to the initial conditions. Sometimes an
        # initial condition cannot be found if the gravitational and/or magnetic
        # field profile is a bit strange.
        if (self.config['initial_conditions']['use_poly_fit_gravity']
                and 'poly_fit_gravity' in self.config['general']):
            files += [('poly_fit.gravity', self.poly_fit_gravity)]
        if (self.config['initial_conditions']['use_poly_fit_magnetic_field']
                and 'poly_fit_magnetic_field' in self.config['general']):
            files += [('poly_fit.magnetic_field', self.poly_fit_magnetic_field)]
        for filename, filestring in files:
            with open(os.path.join(root_dir, filename), 'w') as f:
                f.write(filestring)
        # NOTE: make sure we have needed permissions to run compile script
        self._run_shell_command(
            ['chmod', 'u+x', 'build_initial_conditions.bat'],
            os.path.join(root_dir, 'Initial_Conditions/build_scripts'),
            shell=False,
            verbose=verbose
        )
        self._run_shell_command(
            ['./build_initial_conditions.bat'],
            os.path.join(root_dir, 'Initial_Conditions/build_scripts'),
            verbose=verbose
        )
        if not os.path.exists(os.path.join(root_dir,
                                           'Initial_Conditions/profiles')):
            os.mkdir(os.path.join(root_dir, 'Initial_Conditions/profiles'))
        if execute:
            self._run_shell_command(
                ['./Initial_Conditions.exe'],
                root_dir,
                verbose=verbose,
            )
            if self.config['heating']['background'].get('use_initial_conditions', False):
                self.equilibrium_heating_rate = self.get_equilibrium_heating_rate(root_dir)

    def get_equilibrium_heating_rate(self, root_dir):
        """
        Read equilibrium heating rate from initial conditions results

        # Parameters
        root_dir (`str`): Path to HYDRAD directory
        """
        filename = os.path.join(root_dir,
                                'Initial_Conditions/profiles/initial.amr.sol')
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
            ('Radiation_Model/source/config.h',
             self.radiation_header),
            ('Radiation_Model/config/elements_eq.cfg',
             self.radiation_equilibrium_cfg),
            ('Radiation_Model/config/elements_neq.cfg',
             self.radiation_nonequilibrium_cfg),
            ('Heating_Model/source/config.h',
             self.heating_header),
            ('Heating_Model/config/heating_model.cfg',
             self.heating_cfg),
            ('HYDRAD/source/config.h',
             self.hydrad_header),
            ('HYDRAD/source/collisions.h',
             self.collisions_header),
            ('HYDRAD/config/HYDRAD.cfg',
             self.hydrad_cfg),
        ]
        if 'poly_fit_gravity' in self.config['general']:
            files += [('poly_fit.gravity', self.poly_fit_gravity)]
        if 'poly_fit_magnetic_field' in self.config['general']:
            files += [('poly_fit.magnetic_field', self.poly_fit_magnetic_field)]
        for filename, filestring in files:
            with open(os.path.join(root_dir, filename), 'w') as f:
                f.write(filestring)
        # NOTE: using OpenMP requires an alternate compile script
        if self.config['general'].get('use_openmp', False):
            build_script = 'build_HYDRAD_OPENMP.bat'
        else:
            build_script = 'build_HYDRAD.bat'
        # NOTE: make sure we have needed permissions to run compile script
        self._run_shell_command(
            ['chmod', 'u+x', build_script],
            os.path.join(root_dir, 'HYDRAD/build_scripts'),
            shell=False,
            verbose=verbose,
        )
        self._run_shell_command(
            [f'./{build_script}'],
            os.path.join(root_dir, 'HYDRAD/build_scripts'),
            verbose=verbose,
        )
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
            return datetime.datetime.utcnow().strftime('%Y-%m-%d_%H.%M.%S UTC')

    @property
    def templates(self,):
        """
        List of available templates
        """
        return self.env.list_templates()

    def get_raw_template(self, name):
        """
        Return the unrendered template.
        """
        with open(self.env.get_template(name).filename, 'r') as f:
            return f.read()

    @property
    def initial_conditions_cfg(self):
        """
        Initial conditions configuration file,
        `Initial_Conditions/config/initial_conditions.cfg`
        """
        return self.env.get_template('initial_conditions.cfg').render(
            date=self.date,
            **self.config
        )

    @property
    def initial_conditions_header(self):
        """
        Initial conditions header file, `Initial_Conditions/source/config.h`
        """
        return self.env.get_template('initial_conditions.config.h').render(
            date=self.date,
            maximum_cells=self.maximum_cells,
            minimum_cells=self.minimum_cells,
            **self.config
        )

    @property
    def hydrad_cfg(self):
        """
        HYDRAD configuration file, `HYDRAD/config/hydrad.cfg`
        """
        return self.env.get_template('hydrad.cfg').render(
            date=self.date,
            **self.config
        )

    @property
    def hydrad_header(self):
        """
        HYDRAD header file, `HYDRAD/source/config.h`
        """
        return self.env.get_template('hydrad.config.h').render(
            date=self.date,
            **self.config
        )

    @property
    def heating_cfg(self):
        """
        Heating model configuration file, `Heating_Model/config/heating.cfg`.
        If background heating is enabled and you want to use the equilibrium
        values, you must run the initial conditions and set the
        `equilibrium_heating_rate` attribute first.
        """
        if self.config['heating'].get('background', False):
            bg = self.config['heating']['background']
            if bg.get('use_initial_conditions', False):
                background = {
                    'rate': self.equilibrium_heating_rate,
                    'location': self.config['initial_conditions']['heating_location'],
                    'scale_height': self.config['initial_conditions']['heating_scale_height'],
                }
            elif all((k in bg for k in ('rate', 'location', 'scale_height'))):
                background = self.config['heating']['background']
            else:
                raise ValueError(
                    'Set use_initial_conditions to True or set parameters '
                    'explicitly in order to use background heating.')
        else:
            background = {
                'rate': 0*u.erg/(u.cm**3 * u.s),
                'location': 0*u.cm,
                'scale_height': 0*u.cm
            }
        return self.env.get_template('heating.cfg').render(
            date=self.date,
            background=background,
            **self.config
        )

    @property
    def heating_header(self):
        """
        Heating model header file, `Heating_Model/source/config.h`
        """
        return self.env.get_template('heating.config.h').render(
            date=self.date,
            **self.config
        )

    @property
    def radiation_equilibrium_cfg(self):
        """
        Equilibrium elements configuration file,
        `Radiation_Model/config/elements_eq.cfg`
        """
        elements = self.config['radiation'].get('elements_equilibrium', [])
        return self.env.get_template('radiation.elements.cfg').render(
            date=self.date,
            elements=elements,
            **self.config
        )

    @property
    def radiation_nonequilibrium_cfg(self):
        """
        Non-equilibrium elements configuration file,
        `Radiation_Model/config/elements_neq.cfg`
        """
        elements = self.config['radiation'].get('elements_nonequilibrium', [])
        return self.env.get_template('radiation.elements.cfg').render(
            date=self.date,
            elements=elements,
            **self.config
        )

    @property
    def radiation_header(self):
        """
        Radiation model header file, `Radiation_Model/source/config.h`
        """
        return self.env.get_template('radiation.config.h').render(
            date=self.date,
            **self.config
        )

    @property
    def collisions_header(self):
        """
        Collisions header file, `HYDRAD/source/collisions.h`
        """
        return self.env.get_template('collisions.h').render(
            date=self.date,
            **self.config
        )

    @property
    def poly_fit_magnetic_field(self):
        """
        Sixth-order polynomial fit coefficients for computing flux tube
        expansion
        """
        return self.env.get_template('coefficients.cfg').render(
            date=self.date,
            coefficients=self.config['general']['poly_fit_magnetic_field']
        )

    @property
    def poly_fit_gravity(self):
        """
        Sixth-order polynomial fit coefficients for computing gravitational
        acceleration
        """
        return self.env.get_template('coefficients.cfg').render(
            date=self.date,
            coefficients=self.config['general']['poly_fit_gravity']
        )

    @property
    def minimum_cells(self):
        """
        Minimum allowed number of grid cells,
        $n_{min}=\lceil L/\Delta s_{max}\\rceil$, where $L$ is the loop
        length and $\Delta s_{max}$ is the maximum allowed grid cell width.
        """
        n_min = self.config['general']['loop_length'] / self.config['grid']['maximum_cell_width']
        if n_min.decompose().unit != u.dimensionless_unscaled:
            raise u.UnitConversionError(
                f'''Maximum cell width must be able to be converted to 
                {self.config['general']['loop_length'].unit}''')
        return int(np.ceil(n_min.decompose()))

    @property
    def maximum_cells(self):
        """
        Maximum allowed number of grid cells,
        $n_{max}=\lfloor 2^{L_R}/n_{min}\\rfloor$, where $L_R$ is the maximum
        refinement level and $n_{min}$ is the minimum allowed number of
        grid cells.
        """
        n_min = self.config['general']['loop_length'] / self.config['grid']['maximum_cell_width']
        if n_min.decompose().unit != u.dimensionless_unscaled:
            raise u.UnitConversionError(
                f'''Maximum cell width must be able to be converted to
                {self.config['general']['loop_length'].unit}''')
        return int(np.floor(
            2**self.config['grid']['maximum_refinement_level'] * n_min))
