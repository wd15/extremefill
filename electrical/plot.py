import pylab
from fipy import dump, Grid1D
import parameters
import numpy
import matplotlib
import math
#from matplotlib import rc

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
#rc('text', usetex=True)

matplotlib.rcParams['lines.linewidth'] = 2
#font = {'family' : 'normal',
#        'weight' : 'normal',
#        'size'   : 12}
#matplotlib.rc('font', **font)
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['legend.labelspacing'] = 0.1
matplotlib.rcParams['figure.subplot.wspace'] = 0.3
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
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
    return X, ID

def E(potential):
    return numpy.exp(-parameters.alpha * parameters.Fbar * potential) \
        - numpy.exp((2 - parameters.alpha) * parameters.Fbar * potential)

V = parameters.omega * parameters.i0 * E(parameters.appliedPotential) / parameters.charge / parameters.faradaysConstant
trenchWidth = 2 * parameters.areaRatio / parameters.perimeterRatio

def me(n):
    s = '%1.1e' % n
    m, e = s.split('e')
    return float(m), int(e)

def plotDeposition(variables, fileprefix, label, figprefix, mulFactor=1, legend=1, loc='upper left', maxSuppressor=parameters.bulkSuppressor):

    pylab.figure()
    figDeposition = pylab.subplot(221)
    figTheta = pylab.subplot(222)
    figSuppressor = pylab.subplot(223)
    figCupric = pylab.subplot(224)
    maxFeatureDepth = 0
    print '**************'
    print 'varying ' + figprefix

    xlabel = r'$z$ ($\mu$m)'
    scaleFactor = 1000000

    for variable in variables:
        filename = fileprefix + str(variable) + '.gz'
        data = dump.read(filename)
        featureDepth = data['featureDepth']
        X, ID = getX(featureDepth)
        maxFeatureDepth = max(featureDepth, maxFeatureDepth)
        theta = data['theta'][:ID + 1]
        cupric = data['cupric'][:ID + 1]
        potential = data['potential'][:ID + 1]
        cupric = data['cupric'][:ID + 1]
        suppressor = data['suppressor'][:ID + 1]
        I0 = parameters.i0 + parameters.i1 * theta
        cbar = cupric / parameters.bulkCupric
        current = cbar * I0 * E(-potential)
        depositionRate = parameters.omega * current / parameters.charge / parameters.faradaysConstant   

        pylab.axes(figDeposition)
        pylab.plot(X * scaleFactor, depositionRate, label=label % me(variable * mulFactor))

        pylab.axes(figTheta)
        pylab.plot(X * scaleFactor, theta, label=label % me(variable * mulFactor))

        pylab.axes(figSuppressor)
        pylab.plot(X * scaleFactor, suppressor, label=label % me(variable * mulFactor))

        pylab.axes(figCupric)
        pylab.plot(X * scaleFactor, cupric, label=label % me(variable * mulFactor))
        
        print '------------'
        print figprefix + ' value: ' + str(data[figprefix])
        print 'voltage drop',-data['appliedPotential'] - potential[ID]

    ax = pylab.axes(figDeposition)
    if legend == 1:
        l = pylab.legend()
    pylab.ylabel(r'$v$ (m / s)', rotation='vertical', fontsize=14)
    pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
    pylab.xlim(xmax=0)
    pylab.ylim(ymax=8e-9)
    pylab.ylim(ymin=0.0)
    pylab.title('(a)')

    pylab.axes(figTheta)
    if legend == 2:
        l = pylab.legend(loc=loc)
    pylab.ylabel(r'$\theta$', rotation='vertical', fontsize=12)
    pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
    pylab.xlim(xmax=0)
    pylab.ylim(ymax=1* 1.05)
    pylab.ylim(ymin=0)
    pylab.title('(b)')
    
    pylab.axes(figSuppressor)
    if legend == 3:
        l = pylab.legend(loc=loc)
    pylab.xlabel(xlabel, fontsize=14)
    pylab.ylabel(r'$C_{\text{Supp}}$ (mol / m$^3$)', rotation='vertical', fontsize=14)
    pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
    pylab.xlim(xmax=0)
    pylab.ylim(ymin=0)
    pylab.ylim(ymax=maxSuppressor* 1.05)
    pylab.title('(c)')

    pylab.axes(figCupric)
    if legend == 4:
        l = pylab.legend(loc=loc)
    pylab.xlabel(xlabel, fontsize=14)
    pylab.ylabel(r'$C_{\text{Cu}}$ (mol / m$^3$)', rotation='vertical', fontsize=14)
    pylab.ylim(ymax=parameters.bulkCupric * 1.05)
    pylab.ylim(ymin=0)
    pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
    pylab.xlim(xmax=0)
    pylab.title('(d)')

    l.labelspacing = 0
    l.columnspacing = 0.1

    pylab.savefig(figprefix + '.png')


plotDeposition((0.01, 5., 25., 50., 100., 150., 1000.),
               'tmp/base-kPlus-',
               r'$k^+=%1.1f\times 10^{%i}$ (m$^3$ / mol s)',
               'kPlus')

plotDeposition((1e7, 1.5e7, 2e7, 2.5e7, 3e7),
               'tmp/base-kMinus-',
               r'$k^-=%1.1f\times 10^{%i}$ (1 / m)',
               'kMinus')

plotDeposition((0.001, 0.005, 0.01, 0.02, 0.03, 0.04),
               'tmp/base-deltaRef-',
               r'$\delta_{\text{Ref}}=%1.1f\times 10^{%i}$ (m)',
               'deltaRef',
               legend=2, 
               loc='upper right')

plotDeposition((0.005, 0.01, 0.02, 0.04),
               'tmp/base-bulkSuppressor-',
               r'$C_{\text{Supp}}^{\infty}=%1.1f\times 10^{%i}$ (mol / m$^3$)',
               'bulkSuppressor',
               maxSuppressor=0.04)

plotDeposition((-0.200, -0.250, -0.300),
               'tmp/base-appliedPotential-',
               r'$E_{\text{Applied}}=%1.1f\times 10^{%i}$ (V)',
               'appliedPotential')

plotDeposition((15e-6, 25e-6, 35e-6, 45e-6, 55e-6, 65e-6, 75e-6, 85e-6),
               'tmp/base-featureDepth-',
               r'$h=%1.1f\times 10^{%i}$ ($\mu$m)',
               'featureDepth',
               mulFactor=1000000,
               legend=3)
