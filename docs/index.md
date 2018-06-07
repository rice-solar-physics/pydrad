# hydrad_tools
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

The recommended method for installing Python and the associated dependencies is the [Anaconda Python distribution](https://www.anaconda.com/download/). With Anaconda installed, all of the above dependencies can be installed with,
```shell
$ conda install {package_name}
```
or with pip,
```shell
$ pip install {package_name}
```
Alternatively, you can create a conda environment with all of the needed dependencies using the conda environment file in the repository,
```shell
$ conda create -f conda_environment.yml
$ source activate hydrad_tools
```
Finally, install the hydrad_tools package,
```shell
$ git clone https://github.com/rice-solar-physics/hydrad_tools.git
$ cd hydrad_tools
$ python setup.py install
```

## Additional Resources
Below are a few papers describing the HYDRAD code,

* [Bradshaw and Mason (2003a)](http://adsabs.harvard.edu/abs/2003A%26A...401..699B)
* [Bradshaw and Mason (2003b)](http://adsabs.harvard.edu/abs/2003A%26A...407.1127B)
* [Bradshaw and Klimchuk (2011)](http://adsabs.harvard.edu/abs/2011ApJS..194...26B)
* [Bradshaw and Cargill (2013)](http://adsabs.harvard.edu/abs/2013ApJ...770...12B)
