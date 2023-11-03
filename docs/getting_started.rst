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

Below are a few papers describing the HYDRAD code,

-  `Bradshaw and Mason
   (2003a) <http://adsabs.harvard.edu/abs/2003A%26A...401..699B>`__
-  `Bradshaw and Mason
   (2003b) <http://adsabs.harvard.edu/abs/2003A%26A...407.1127B>`__
-  `Bradshaw and Klimchuk
   (2011) <http://adsabs.harvard.edu/abs/2011ApJS..194...26B>`__
-  `Bradshaw and Cargill
   (2013) <http://adsabs.harvard.edu/abs/2013ApJ...770...12B>`__
