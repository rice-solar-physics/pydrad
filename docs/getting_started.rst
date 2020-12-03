Getting Started
===============

Installing this Package
-----------------------

First, clone the repository,

.. code:: shell

   $ git clone https://github.com/rice-solar-physics/pydrad.git
   $ cd pydrad

To install all of the needed dependencies,

.. code:: shell

   $ pip install -r requirements/requirements/txt

Finally, install the pydrad package,

.. code:: shell

   $ python setup.py install

Testing
-------

If you’d like to run the tests, you can install the additional
development dependencies,

.. code:: shell

   $ pip install -r requirements/requirements-dev.txt

and run the tests,

.. code:: shell

   $ pytest

If you’d like to also run the tests that require a local copy of HYDRAD,

.. code:: shell

   $ pytest --hydrad-dir=/path/to/hydrad

Note that this step is not necessary to run the code.

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
