# pydrad

[![Build Status](https://travis-ci.org/rice-solar-physics/pydrad.svg?branch=master)](https://travis-ci.org/rice-solar-physics/pydrad)
[![codecov](https://codecov.io/gh/rice-solar-physics/pydrad/branch/master/graph/badge.svg)](https://codecov.io/gh/rice-solar-physics/pydrad)

Some Python tools to configure and parse output from the HYDrodynamics and RADiation (HYDRAD) code for field-aligned coronal loop physics.

## Install

To install the package and the needed dependencies,
```shell
$ git clone https://github.com/rice-solar-physics/pydrad.git
$ cd pydrad
$ pip install -r requirements/requirements.txt
$ python setup.py install
```

If you'd like to run the tests and confirm that everything is working alright,
```shell
$ pip install -r requirements/requirements-dev.txt
$ pytest
``` 

See the [docs](https://rice-solar-physics.github.io/pydrad/) for more info. Additionally, **you will need access to the HYDRAD source code.**

## Help
Create an issue if you run into any problems. Submit a PR if you would like to add any functionality.
