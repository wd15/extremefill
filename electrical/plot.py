import pylab
from fipy import dump, Grid1D
import parameters
import numpy
import matplotlib

matplotlib.rcParams['lines.linewidth'] = 2
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
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
    return X, ID

def E(potential):
    return numpy.exp(-parameters.alpha * parameters.Fbar * potential) \
        - numpy.exp((2 - parameters.alpha) * parameters.Fbar * potential)

V = parameters.omega * parameters.i0 * E(parameters.appliedPotential) / parameters.charge / parameters.faradaysConstant
trenchWidth = 2 * parameters.areaRatio / parameters.perimeterRatio

def plotDeposition(variables, fileprefix, label, figprefix, mulFactor=1):

    figDeposition = pylab.figure()
    figTheta = pylab.figure()
    figSuppressor = pylab.figure()
    figCupric = pylab.figure()
    maxFeatureDepth = 0
    print '**************'
    print 'varying ' + figprefix

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

        pylab.figure(figDeposition.number)
        pylab.plot(X * 1000, depositionRate, label=label % (variable * mulFactor))

        pylab.figure(figTheta.number)
        pylab.plot(X * 1000, theta, label=label % (variable * mulFactor))

        pylab.figure(figSuppressor.number)
        pylab.plot(X * 1000, suppressor, label=label % (variable * mulFactor))

        pylab.figure(figCupric.number)
        pylab.plot(X * 1000, cupric, label=label % (variable * mulFactor))
        
        print '------------'
        print figprefix + ' value: ' + str(data[figprefix])
        print 'voltage drop',-data['appliedPotential'] - potential[ID]

    pylab.figure(figDeposition.number)
    pylab.xlabel(r'z (mm)', fontsize=16)
    pylab.ylabel(r'v (m / s)', rotation='vertical', fontsize=16)
    pylab.legend()
    pylab.xlim(xmin=-maxFeatureDepth * 1000)
    pylab.xlim(xmax=0)
    pylab.savefig(figprefix + 'Deposition.png')

    pylab.figure(figTheta.number)
    pylab.xlabel(r'z (mm)', fontsize=16)
    pylab.ylabel(r'$\theta$', rotation='horizontal', fontsize=16)
    pylab.legend(loc='lower right')
    pylab.xlim(xmin=-maxFeatureDepth * 1000)
    pylab.xlim(xmax=0)
    pylab.ylim(ymax=1)
    pylab.ylim(ymin=0)
    pylab.savefig(figprefix + 'Theta.png')

    pylab.figure(figSuppressor.number)
    pylab.xlabel(r'z (mm)', fontsize=16)
    pylab.ylabel(r'$C_{S}$', rotation='horizontal', fontsize=16)
    pylab.legend(loc='upper left')
    pylab.xlim(xmin=-maxFeatureDepth * 1000)
    pylab.xlim(xmax=0)
    pylab.ylim(ymin=0)
    pylab.savefig(figprefix + 'Suppressor.png')

    pylab.figure(figCupric.number)
    pylab.xlabel(r'z (mm)', fontsize=16)
    pylab.ylabel(r'$C_{Cu}$', rotation='horizontal', fontsize=16)
    pylab.ylim(ymax=parameters.bulkCupric)
    pylab.ylim(ymin=0)
    pylab.legend(loc='lower right')
    pylab.xlim(xmin=-maxFeatureDepth * 1000)
    pylab.xlim(xmax=0)
    pylab.savefig(figprefix + 'Cupric.png')


plotDeposition((0.01, 5., 25., 50., 100., 150., 1000.),
               'tmp/base-kPlus-',
               r'$k^+=%1.1e$ (m$^3$ / mol s)',
               'kPlus')

plotDeposition((1e7, 1.5e7, 2e7, 2.5e7, 3e7),
               'tmp/base-kMinus-',
               r'$k^-=%1.1e$ (1 / m)',
               'kMinus')

plotDeposition((0.001, 0.01, 0.02, 0.03, 0.04),
               'tmp/base-deltaRef-',
               r'$\delta_{ref}=%1.1e$ (m)',
               'deltaRef')

plotDeposition((0.005, 0.01, 0.02, 0.04, 0.08),
               'tmp/base-bulkSuppressor-',
               r'$C_{S}^{\infty}=%1.1e$ (mol / m$^3$)',
               'bulkSuppressor')

plotDeposition((-0.200, -0.250, -0.300),
                'tmp/base-appliedPotential-',
                r'$\eta=%1.1e$ (V)',
                'appliedPotential')

plotDeposition((15e-6, 25e-6, 35e-6, 45e-6, 55e-6, 65e-6, 75e-6, 85e-6),
                'tmp/base-featureDepth-',
                r'$h=%1.1e$ (mm)',
                'featureDepth',
               mulFactor=1000)
