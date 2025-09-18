.. _pydrad-getting-started:

Getting Started
===============

Installing this Package
-----------------------

First, clone the repository,

.. code:: shell

   $ git clone https://github.com/rice-solar-physics/pydrad.git
   $ cd pydrad

To install pydrad and all of the needed dependencies,

.. code:: shell

   $ pip install .

Testing
-------

If you’d like to run the tests, you can install the additional
development dependencies,

.. code:: shell

   $ pip install -e .[tests,docs]

and run the tests,

.. code:: shell

   $ pytest pydrad

If you would like to run the tests using a local copy of HYDRAD,

.. code:: shell

   $ pytest --hydrad-dir=/path/to/hydrad

If you do not specify a path to a particular version of HYDRAD, it will be automatically downloaded from GitHub in order to run the tests.
Note that running the tests is not necessary if you just want to use pydrad.

Additional Resources
--------------------

A description of the HYDRAD code is available on its Github repository,

-  `HYDRAD user's guide
   <https://github.com/rice-solar-physics/HYDRAD/blob/6344b8e3e14ba7d3d470f9a5d57b0adc16421731/HYDRAD_User_Guide(03_20_2021).pdf>`__

Below are a few papers describing the HYDRAD code,

- :cite:t:`bradshaw_self-consistent_2003`
- :cite:t:`bradshaw_radiative_2003`
- :cite:t:`bradshaw_what_2011`
- :cite:t:`bradshaw_influence_2013`
- :cite:t:`reep_efficient_2019`
