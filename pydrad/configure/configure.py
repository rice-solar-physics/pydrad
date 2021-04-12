"""
Configure HYDRAD simulations
"""
import os
import copy
import datetime
import tempfile
import shutil
from distutils.dir_util import copy_tree

import numpy as np
import astropy.units as u
from jinja2 import Environment, PackageLoader, ChoiceLoader, DictLoader
import asdf

from . import filters
from .util import run_shell_command, on_windows, get_equilibrium_heating_rate

__all__ = ['Configure']


class Configure(object):
    """
    Configure HYDRAD simulations from a single Python `dict`

    Parameters
    ----------
    config : `dict`
        All input parameters for configuring simulation
    templates : `dict`
        Templates to override defaults, optional
    """

    def __init__(self, config, **kwargs):
        self.config = copy.deepcopy(config)
        loader = ChoiceLoader([
            DictLoader(kwargs.get('templates', {})),
            PackageLoader('pydrad', 'configure/templates')
        ])
        self.env = Environment(loader=loader)
        self.env.filters['units_filter'] = filters.units_filter
        self.env.filters['log10_filter'] = filters.log10_filter
        self.env.filters['get_atomic_symbol'] = filters.get_atomic_symbol
        self.env.filters['get_atomic_number'] = filters.get_atomic_number
        self.env.filters['sort_elements'] = filters.sort_elements
        self.env.filters['is_required'] = filters.is_required
        self.env.filters['sci_notation'] = filters.sci_notation
        # NOTE: Freeze the date at instantiation so that files can be compared
        # exactly for testing
        if kwargs.get('freeze_date', False):
            self._date = self.date
            self._freeze_date = True

    @staticmethod
    def load_config(filename):
        """
        Load a base configuration from an ASDF file

        Parameters
        ----------
        filename : `str`
        """
        with asdf.open(filename) as af:
            config = copy.deepcopy(dict(af.tree))
        return config

    def save_config(self, filename):
        """
        Save the simulation configuration as an ASDF file

        Parameters
        ----------
        filename : `str`
        """
        asdf.AsdfFile(self.config).write_to(filename)

    def setup_simulation(self,
                         output_path,
                         base_path,
                         run_initial_conditions=True,
                         overwrite=False,
                         **kwargs):
        """
        Setup a HYDRAD simulation with desired outputs from a clean copy

        Parameters
        ----------
        output_path : `str`
            Path to new copy of HYDRAD
        base_path : `str`
            Path to existing HYDRAD
        run_initial_conditions : `bool`
            If True, compile and run the initial conditions code
        overwrite : `bool`
            If True, overwrite an existing HYDRAD instance at ``output_path``
            if it exists.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # NOTE: this is all done in a temp directory and then copied over
            # so that if something fails, all the files are cleaned up
            copy_tree(base_path, tmpdir)
            if run_initial_conditions:
                execute = kwargs.get('execute', True)
                self.setup_initial_conditions(tmpdir, execute=execute)
            self.setup_hydrad(tmpdir)
            self.save_config(os.path.join(tmpdir, 'pydrad_config.asdf'))
            if overwrite:
                if os.path.exists(output_path):
                    shutil.rmtree(output_path)
            shutil.copytree(tmpdir, output_path)

    def setup_initial_conditions(self, root_dir, execute=True):
        """
        Compile and run the initial conditions code to get the initial
        loop profile

        Parameters
        ----------
        root_dir : `str`
            Path to new HYDRAD copy
        execute : `bool`
            If True (default), compute initial conditions. Otherwise, they
            are only compiled. This is useful for debugging.
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
        # Only do this on Unix-based systems
        if not on_windows():
            run_shell_command(
                ['chmod', 'u+x', 'build_initial_conditions.bat'],
                os.path.join(root_dir, 'Initial_Conditions/build_scripts'),
                shell=False
            )
        run_shell_command(
            ['./build_initial_conditions.bat'],
            os.path.join(root_dir, 'Initial_Conditions/build_scripts')
        )
        if not os.path.exists(os.path.join(root_dir,
                                           'Initial_Conditions/profiles')):
            os.mkdir(os.path.join(root_dir, 'Initial_Conditions/profiles'))
        if execute:
            run_shell_command(
                ['./Initial_Conditions.exe'],
                root_dir,
            )
            if self.config['heating']['background'].get('use_initial_conditions', False):
                self.equilibrium_heating_rate = get_equilibrium_heating_rate(root_dir)

    def setup_hydrad(self, root_dir):
        """
        Compile HYDRAD code with appropriate header and config files

        Parameters
        -----------
        root_dir : `str`
            Path to new HYDRAD copy
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
            ('HYDRAD/config/hydrad.cfg',
             self.hydrad_cfg),
        ]
        if 'poly_fit_gravity' in self.config['general']:
            files += [('poly_fit.gravity', self.poly_fit_gravity)]
        if 'poly_fit_magnetic_field' in self.config['general']:
            files += [('poly_fit.magnetic_field', self.poly_fit_magnetic_field)]
        if self.config['heating'].get('beam', False):
            files += [('Heating_Model/config/beam_heating_model.cfg',
                       self.beam_heating_cfg)]
        for filename, filestring in files:
            with open(os.path.join(root_dir, filename), 'w') as f:
                f.write(filestring)
        # NOTE: using OpenMP requires an alternate compile script
        if self.config['general'].get('use_openmp', False):
            build_script = 'build_HYDRAD_OPENMP.bat'
        else:
            build_script = 'build_HYDRAD.bat'
        # NOTE: make sure we have needed permissions to run compile script
        # Only do this on Unix-based systems (chmod not valid on Windows)
        if not on_windows():
            run_shell_command(
                ['chmod', 'u+x', build_script],
                os.path.join(root_dir, 'HYDRAD/build_scripts'),
                shell=False,
            )
        run_shell_command(
            [f'./{build_script}'],
            os.path.join(root_dir, 'HYDRAD/build_scripts'),
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
        background = copy.deepcopy(self.config['heating']['background'])
        if background.get('use_initial_conditions', False):
            background['rate'] = self.equilibrium_heating_rate
            background['location'] = self.config['initial_conditions']['heating_location']
            background['scale_height'] = self.config['initial_conditions']['heating_scale_height']
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
    def beam_heating_cfg(self):
        """
        Beam heating model configuration file.
        """
        return self.env.get_template('heating.beam.cfg').render(
            date=self.date,
            **self.config,
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
        fit = self._fit_poly_domains('poly_fit_magnetic_field', 'G')
        return self.env.get_template('coefficients.cfg').render(
            date=self.date,
            **fit,
        )

    @property
    def poly_fit_gravity(self):
        """
        Sixth-order polynomial fit coefficients for computing gravitational
        acceleration
        """
        fit = self._fit_poly_domains('poly_fit_gravity', 'cm s-2')
        return self.env.get_template('coefficients.cfg').render(
            date=self.date,
            **fit,
        )

    def _fit_poly_domains(self, name, unit):
        """
        Perform polynomial fit to quantity as a function of field aligned coordinate
        over multiple domains and return fitting coefficients.
        """
        # TODO: refactor to be independent of dictionary
        fit = copy.deepcopy(self.config['general'][name])
        x = (fit['x'] / self.config['general']['loop_length']).decompose().to(u.dimensionless_unscaled).value
        y = fit['y'].to(unit).value
        coefficients = []
        minmax = []
        for i in range(len(fit['domains'])-1):
            i_d = np.where(np.logical_and(x>=fit['domains'][i], x<=fit['domains'][i+1]))
            coefficients.append(np.polyfit(x[i_d], y[i_d], fit['order'])[::-1])
            minmax.append([y[i_d].min(), y[i_d].max()])
        fit['minmax'] = minmax
        fit['coefficients'] = coefficients
        return fit

    @property
    def minimum_cells(self):
        """
        Minimum allowed number of grid cells,
        :math:`n_{min}=\lceil L/\Delta s_{max}\\rceil`, where :math:`L` is the loop
        length and :math:`\Delta s_{max}` is the maximum allowed grid cell width.
        Optionally, if the minimum number of cells is specified
        in ``config['grid']['minimum_cells']``, this value will take
        precedence.
        """
        if 'minimum_cells' in self.config['grid']:
            return int(self.config['grid']['minimum_cells'])
        n_min = self.config['general']['loop_length'] / self.config['grid']['maximum_cell_width']
        if n_min.decompose().unit != u.dimensionless_unscaled:
            raise u.UnitConversionError(
                f'''Maximum cell width must be able to be converted to 
                {self.config['general']['loop_length'].unit}''')
        return int(np.round(n_min.decompose()))

    @property
    def maximum_cells(self):
        """
        Maximum allowed number of grid cells,
        :math:`n_{max}=\lfloor 2^{L_R}n_{min}\\rfloor`, where :math:`L_R` is the maximum
        refinement level and :math:`n_{min}` is the minimum allowed number of
        grid cells. Optionally, if the maximum number of cells is specified
        in ``config['grid']['maximum_cells']``, this value will take
        precedence.
        """
        if 'maximum_cells' in self.config['grid']:
            return int(self.config['grid']['maximum_cells'])
        n_min = self.config['general']['loop_length'] / self.config['grid']['maximum_cell_width']
        if n_min.decompose().unit != u.dimensionless_unscaled:
            raise u.UnitConversionError(
                f'''Maximum cell width must be able to be converted to
                {self.config['general']['loop_length'].unit}''')
        return int(np.floor(
            2**self.config['grid']['maximum_refinement_level'] * n_min))
