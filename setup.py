import os
from distutils.core import setup

setup(
    name='hydrad_tools',
    license='MIT',
    version='0.1',
    author='Will Barnes',
    url='https://github.com/rice-solar-physics/hydrad_tools',
    package_data={'hydrad_tools': ['configure/templates/*']},
    packages=['hydrad_tools'],
    author_email='will.t.barnes@gmail.com',
    description='Tools for configuring and parsing HYDRAD simulations'
)
