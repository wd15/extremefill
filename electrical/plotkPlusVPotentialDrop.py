import numpy as np
import parameters
from fipy import Grid1D, dump
import pylab
import matplotlib

kMinus = 2.47707636e+07
kPluses = 10**np.linspace(0, 3, 100)
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['legend.fontsize'] = 11

def getX(featureDepth):
    L = parameters.delta + featureDepth
    N = 1000
    dx = L / N 
    mesh = Grid1D(nx=N, dx=dx) - [[featureDepth]]
    x = 0
    ID = np.argmin(abs(mesh.x - x))
    if mesh.x[ID] < x:
        ID = ID + 1
    X = mesh.x[:ID + 1].value
    return X, ID, dx

X, ID, dx = getX(parameters.featureDepth)

def potentialDrop(filename):
    data = dump.read(filename)
    return data['potential'][ID]

potentials = np.zeros(len(kPluses), 'd')

for i, kPlus in enumerate(kPluses):
    filename='tmp/base-kPlus-%1.2e-kMinus-%1.2e' % (kPlus, kMinus)
    potentials[i] = potentialDrop(filename)

pylab.semilogx(kPluses, potentials, 'k', lw=2)
pylab.xlabel(r'$k^+$ (m$^3$ / mol s)', fontsize=16)
pylab.ylabel(r'$\eta$ (V)', fontsize=16)
pylab.savefig('kPlusVPotential.png')