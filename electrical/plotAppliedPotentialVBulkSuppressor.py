import pylab
from fipy import dump, Grid1D
import parameters
import numpy

import matplotlib
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
rc('text', usetex=True)

matplotlib.rcParams['lines.linewidth'] = 2
#font = {'family' : 'normal',
#        'weight' : 'normal',
#        'size'   : 12}
#matplotlib.rc('font', **font)
matplotlib.rcParams['legend.fontsize'] = 11

def getX(featureDepth):
    L = parameters.delta + featureDepth
    N = 1000
    dx = L / N 
    mesh = Grid1D(nx=N, dx=dx) - [[featureDepth]]
    x = 0
    ID = numpy.argmin(abs(mesh.x - x))
    if mesh.x[ID] < x:
        ID = ID + 1
    X = mesh.x[:ID + 1].value
    return X, ID, dx

def E(potential):
    return numpy.exp(-parameters.alpha * parameters.Fbar * potential) \
        - numpy.exp((2 - parameters.alpha) * parameters.Fbar * potential)

V = parameters.omega * parameters.i0 * E(parameters.appliedPotential) / parameters.charge / parameters.faradaysConstant

X, ID, dx = getX(parameters.featureDepth)
def figureOfMerit(filename):
    data = dump.read(filename)
    theta = data['theta'][:ID + 1]
    cupric = data['cupric'][:ID + 1]
    potential = data['potential'][:ID + 1]
    cupric = data['cupric'][:ID + 1]
    I0 = parameters.i0 + parameters.i1 * theta
    cbar = cupric / parameters.bulkCupric
    current = cbar * I0 * E(-potential)
    depositionRate = parameters.omega * current / parameters.charge / parameters.faradaysConstant
    derivative = (depositionRate[2:] - depositionRate[:-2]) / dx / 2
    depositionRate = depositionRate[1:-1]
    mask = (depositionRate > (V / 10000)) & (depositionRate > (max(depositionRate) / 100))
    
    if mask.any():
        if mask.all():
            localized = 0
        else:
            localized = 1
        return max(derivative[mask]), min(derivative[mask]), max(depositionRate), localized
    else:
        return -1e-10, -1e-10, max(depositionRate), 0

def getDepositionRate(filename):
    data = dump.read(filename)
    theta = data['theta'][:ID + 1]
    cupric = data['cupric'][:ID + 1]
    potential = data['potential'][:ID + 1]
    I0 = parameters.i0 + parameters.i1 * theta
    cbar = cupric / parameters.bulkCupric
    current = cbar * I0 * E(-potential)
    return X, parameters.omega * current / parameters.charge / parameters.faradaysConstant

##kPlus=150
##kMinus=1e8
##filename='tmp/base-kPlus-' + str(kPlus) + '-kMinus-' + str(kMinus) + '.gz'
##print figureOfMerit(filename) 
##raw_input('stopped')

# ax = pylab.subplot(111)
# ax.set_xscale('log')
# ax.set_yscale('log')
# def log_10_product(x, pos):
#     """The two args are the value and tick position.
#     Label ticks with the product of the exponentiation"""
#     return '%1i' % (x)
# formatter = pylab.FuncFormatter(log_10_product)
# ax.xaxis.set_major_formatter(formatter)
# ax.yaxis.set_major_formatter(formatter)

from  potentialVsuppressorStudy import appliedPotentials, bulkSuppressors
figureOfMerits = numpy.zeros((len(appliedPotentials), len(bulkSuppressors), 4), 'd')
#print bulkSuppressors
#print appliedPotentials
#raw_input()
for i, appliedPotential in enumerate(appliedPotentials):
    if appliedPotential < -9.9e-1:
        break
    for j, bulkSuppressor in enumerate(bulkSuppressors):
        filename='tmp/base-appliedPotential-%1.3e-bulkSuppressor-%1.3e' % (appliedPotential, bulkSuppressor)
        figureOfMerits[i, j, :] = figureOfMerit(filename)

a = pylab.contourf(bulkSuppressors, -appliedPotentials, figureOfMerits[...,0], (-1e+10, 1e+10), colors=('red',))
a = pylab.contourf(bulkSuppressors, -appliedPotentials, figureOfMerits[...,1], (-1e+10, 0), colors=('orange',))
a = pylab.contourf(bulkSuppressors, -appliedPotentials, figureOfMerits[...,0], (-1e+10, 0), colors=('green',))
a = pylab.contourf(bulkSuppressors, -appliedPotentials, figureOfMerits[...,2], (-1e+10, V / 10000), colors=('blue',))
a = pylab.contour(bulkSuppressors, -appliedPotentials, figureOfMerits[...,3], (-1e+10, 0.5, 1e+10), linestyles='solid', linewidths=4, colors='k')

a.ax.set_xscale('log')
a.ax.set_yscale('log')

##pylab.loglog((0.02,), (0.25,), 'ko', lw=3) 

pylab.xlabel(r'$C_S$ (mol / m$^3$)', fontsize=14)
pylab.ylabel(r'$-\eta$ (V)', fontsize=14)

align = False

def makeBackGroudPlot(appliedPotential, bulkSuppressor, fP, fM):
    xpos = (numpy.log10(bulkSuppressor) + 4.) / 4. *  fM
    ypos = (numpy.log10(-appliedPotential) + 2) / 2. * fP   
    pylab.axes(a.ax)
    if align:
        pylab.loglog((bulkSuppressor,), (-appliedPotential,), 'ko', lw=3) 
    la = pylab.axes((xpos, ypos, 0.08, 0.08), frame_on=True, axisbg='w')
    la.patch.set_alpha(0.5)
    X, depositionRate = getDepositionRate('tmp/base-appliedPotential-%1.3e-bulkSuppressor-%1.3e' % (appliedPotential, bulkSuppressor))

    pylab.plot(X, depositionRate, 'k', lw=2)
    pylab.setp(la, ylim=(0,1e-8), xlim=(-parameters.featureDepth, 0), xticks=[], yticks=[], alpha=0.8)
    if align:
        pylab.plot((-parameters.featureDepth / 2,), (1e-8 / 2,), 'ks', lw=3) 

for fP, appliedPotential in ((0.885, -2.477e-1), (1.1, -0.02535364)):
    for fM, bulkSuppressor in ((1.8, 2.10490414e-04), (1.04, 1.96304065e-03), (0.925, 2.009e-02), (0.88, 2.05651231e-01)):
        makeBackGroudPlot(appliedPotential, bulkSuppressor, fP, fM)

makeBackGroudPlot(-0.097701, 6.57933225e-03, 0.927, 0.952)
makeBackGroudPlot(-0.39442061, 6.73415066e-02, 0.875, 0.895)
pylab.savefig('appliedPotentialVBulkSuppressor.png')
pylab.show()
