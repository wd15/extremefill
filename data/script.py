import sys
import os
import argparse

import tables
from extremefill.simulation2D import Simulation2D


## Parse command line arguments.
class Bare(object):
    pass


args = Bare()
parser = argparse.ArgumentParser(description="Run Extremefill 2D example.")
parser.add_argument('parameterfile', default='default.param')
parser.parse_args(namespace=args)

## Generate parameters.
parameters = Bare()
f = open(args.parameterfile, 'r')
for line in f.readlines():
    print line
    exec('parameters.' + line)

## Generate data file path.
datadir = 'Data'
datafile = 'data.h5'
if hasattr(parameters, 'sumatra_label'):
    datadir = os.path.join(datadir, parameters.sumatra_label)    
datapath = os.path.join(datadir, datafile)

simulation = Simulation2D()
simulation.run(view=False,
               totalSteps=parameters.steps,
               sweeps=30,
               dt=0.01,
               tol=1e-1,
               Nx=parameters.Nx,
               CFL=parameters.CFL,
               PRINT=True,
               areaRatio=2 * 0.093,
               dtMax=100.,
               dataFile=datapath,
               totalTime=5000.,
               data_frequency=10,
               NxBase=1200)


