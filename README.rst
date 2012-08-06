========
Overview
========

This python package reproduces the data required for the figures in
"Modeling Extreme Bottom-Up Filling of through Silicon Vias" by
Josell, Wheeler and Moffat.

------------
Requirements
------------

The main requirements are a working version of FiPy. See
requirements.txt for the specific versions of dependencies. Generally
these won't matter, but the correct version of FiPy may be required.

------------
Installation
------------

Use::

$ python setup.py install

-----
Usage
-----

To generate all the data required to generate the figures.

>>> import extremefill
>>> extremefill.generateData()
 
This will run all the simulations for the paper and takes about a
day. The data is written to an hdf5 file ``data.hdf5``, which will
appear in the base directory. To generate all the figures from the
paper run.

>>> extremefill.generateFigures()

This will generate both pdf and png images of the figures in the paper
using the hdf5 file or generating data that is missing from the hdf5
file. To just generate a single figure.

>>> extremefill.generateFigure(5)

To run a simulation with varying parameters.

>>> extremefill.run(**parameters)

where a variety of parameters can be passed. See the documentation for
further details.


