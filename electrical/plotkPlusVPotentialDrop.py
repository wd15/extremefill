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

def plotkPlusVPotential():
    pylab.semilogx(kPluses, potentials, 'k', lw=1)
    pylab.semilogx((1, 1000), (0.25, 0.25), 'k--', lw=1)
    pylab.xlabel(r'$k^+$', fontsize=10, labelpad=-3)
    pylab.xticks((1, 10, 100, 1000), (r'$1$', r'$10$', r'$100$', r'$1000$'), fontsize=8)
    pylab.yticks((0.1, 0.2), (r'$-0.1$', r'$-0.2$'), fontsize=8)
    pylab.ylim(0.04, 0.27)

    pylab.ylabel(r'$\eta$', fontsize=10, rotation='horizontal', labelpad=-8)
##    pylab.title(r'(e)', fontsize=10)


if __name__ == '__main__':
    plotkPlusVPotential()
    pylab.savefig('kPlusVPotential.png')
    pylab.show()
