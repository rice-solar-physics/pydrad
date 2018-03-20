"""
Configure HYDRAD simulations
"""
import os
import datetime

import numpy as np
from jinja2 import Environment, PackageLoader
import yaml

from . import filters

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
        # Setup paths for base simulation
        # Setup paths for output simulation
        self.env = Environment(loader=PackageLoader('hydrad_tools', 'configure/templates'))
        self.env.filters['units_filter'] = filters.units_filter
        self.env.filters['log10_filter'] = filters.log10_filter
        self.env.filters['get_atomic_symbol'] = filters.get_atomic_symbol
        self.env.filters['get_atomic_number'] = filters.get_atomic_number
        self.env.filters['sort_elements'] = filters.sort_elements

    @property
    def date(self):
        return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

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
