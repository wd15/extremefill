==============
 Extreme Fill
==============

This Python package generates the data required to reproduce the
figures in `Modeling Extreme Bottom-Up Filling of through Silicon
Vias`_ by Josell, Wheeler and Moffat.

Requirements
============

The main requirements are a working version of FiPy_ to solve the PDEs
and PyTables_ to store the data in an HDF5 file. See
`requirements.txt`_ for specific versions of dependencies used to run
this package on the maintainer's system. The `requirements.txt`_ file
is auto-generated so most of the packages listed are not necessarily
required. If you are having issues with installation and missing
packages this would be a good place to start looking. If the plots are
not generated correctly for some reason or latex errors are occurring,
it might be work checking the maintainer's matplotlibrc_ file for
clues on the correct latex packages to install and the correct
matplotlib configuration.

Installation
============

The following should get you most of the way to an install.

::

$ pip install numpy
$ pip install tables
$ pip install pysparse
$ pip install fipy
$ pip install matplotlib

and then clone the git repository

::

$ git clone git://github.com/wd15/extremefill.git

See the `Github project page`_ for further details. After cloning,
install with

::

$ cd extremefill
$ python setup.py install

Documentation
=============

Use::

$ cd doc
$ make html

to generate this documentation.

.. _Modeling Extreme Bottom-Up Filling of through Silicon Vias: http://dx.doi.org/10.1149/2.009210jes
.. _requirements.txt: https://github.com/wd15/extremefill/blob/master/requirements.txt
.. _FiPy: http://www.ctcms.nist.gov/fipy/
.. _Github project page: https://github.com/wd15/extremefill
.. _PyTables: http://www.pytables.org/moin
.. _matplotlibrc: https://github.com/wd15/env/blob/021e67f5acf1344a727f3b9eb012d9f615856f23/matplotlibrc
