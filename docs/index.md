# pydrad
This package contains various helpers for configuring, parsing, and plotting HYDRAD simulations. The HYDrodynamics and RADiation model is a code for modeling the field-aligned dynamics of coronal loops.

## Installing this Package
Python 3.6+ and the following dependencies are required to run this package:

* numpy
* jinja2
* asdf (pip only)
* astropy
* matplotlib
* plasmapy (pip only)
* GitPython (optional)

First, clone the repository,
```shell
$ git clone https://github.com/rice-solar-physics/pydrad.git
$ cd pydrad
```

To install all of the needed dependencies,
```shell
$ pip install -r requirements/requirements/txt
```

Finally, install the pydrad package,
```shell
$ python setup.py install
```

## Testing
If you'd like to run the tests, you can install the additional development dependencies,
```shell
$ pip install -r requirements/requirements-dev.txt
```
and run the tests,
```
$ pytest
```
Note that this step is not necessary to run the code.

## Additional Resources
Below are a few papers describing the HYDRAD code,

* [Bradshaw and Mason (2003a)](http://adsabs.harvard.edu/abs/2003A%26A...401..699B)
* [Bradshaw and Mason (2003b)](http://adsabs.harvard.edu/abs/2003A%26A...407.1127B)
* [Bradshaw and Klimchuk (2011)](http://adsabs.harvard.edu/abs/2011ApJS..194...26B)
* [Bradshaw and Cargill (2013)](http://adsabs.harvard.edu/abs/2013ApJ...770...12B)
