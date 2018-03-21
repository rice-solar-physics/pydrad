"""
Configure HYDRAD simulations
"""
import os
import datetime
import tempfile
import shutil
from distutils.dir_util import copy_tree

import numpy as np
from jinja2 import Environment, PackageLoader
import yaml
import git

from . import filters

REMOTE_REPO = 'https://github.com/rice-solar-physics/HYDRAD'

__all__ = ['Configure']


class Configure(object):

    def __init__(self, config, use_default_options=True):
        if use_default_options:
            with open(os.path.join(os.environ['HOME'], '.hydrad_tools', 'defaults.yml')) as f:
                self.config = yaml.load(f)
            for k in self.config:
                if k in config:
                    self.config[k].update(config[k])
        else:
            self.config = config
        self.env = Environment(loader=PackageLoader('hydrad_tools', 'configure/templates'))
        self.env.filters['units_filter'] = filters.units_filter
        self.env.filters['log10_filter'] = filters.log10_filter
        self.env.filters['get_atomic_symbol'] = filters.get_atomic_symbol
        self.env.filters['get_atomic_number'] = filters.get_atomic_number
        self.env.filters['sort_elements'] = filters.sort_elements

    def setup_simulation(self, output_path, base_path=None, name=None):
        """
        Setup a HYDRAD simulation with desired outputs from a clean copy

        Parameters
        ----------
        output_path : `str`
        base_path : `str`, optional
            If None (default), clone a new copy from GitHub (appropriate permissions required)
        name : `str`, optional
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Get clean copy
            if base_path is None:
                git.Repo.clone_from(REMOTE_REPO, tmpdir)
            else:
                copy_tree(base_path, tmpdir)
            # Generate configuration files and copy them to the right locations
            # Compile executables
            # Copy to output
            if name is None:
                name = f'hydrad_{self.date}'
            output_dir = os.path.join(output_path, name)
            shutil.copytree(tmpdir, output_dir)

    @property
    def date(self):
        return datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')

    @property
    def intial_conditions_cfg(self):
        return self.env.get_template('initial_conditions.cfg').render(date=self.date, **self.config)

    @property
    def initial_conditions_header(self):
        return self.env.get_template('initial_conditions.config.h').render(
                    date=self.date, **self.config)

    @property
    def hydrad_cfg(self):
        return self.env.get_template('hydrad.cfg').render(date=self.date, **self.config)

    @property
    def hydrad_header(self):
        return self.env.get_template('hydrad.config.h').render(date=self.date, **self.config)

    @property
    def heating_cfg(self):
        return self.env.get_template('heating.cfg').render(date=self.date, **self.config)

    @property
    def heating_header(self):
        return self.env.get_template('heating.config.h').render(date=self.date, **self.config)

    @property
    def radiation_equilibrium_cfg(self):
        elements = self.config['radiation'].get('elements_equilibrium', [])
        return self.env.get_template('radiation.elements.cfg').render(
                    date=self.date, elements=elements, **self.config)

    @property
    def radiation_nonequilibrium_cfg(self):
        elements = self.config['radiation'].get('elements_nonequilibrium', [])
        return self.env.get_template('radiation.elements.cfg').render(
                    date=self.date, elements=elements, **self.config)

    @property
    def radiation_header(self):
        return self.env.get_template('radiation.config.h').render(date=self.date, **self.config)

    @property
    def collisions_header(self):
        return self.env.get_template('collisions.h').render(date=self.date, **self.config)
