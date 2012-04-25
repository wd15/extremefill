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

L = parameters.delta + parameters.featureDepth
N = 1000
dx = L / N 
mesh = Grid1D(nx=N, dx=dx) - [[parameters.featureDepth]]
x = 0
ID = numpy.argmin(abs(mesh.x - x))
if mesh.x[ID] < x:
    ID = ID + 1
X = mesh.x[:ID + 1].value

def E(potential):
    return numpy.exp(-parameters.alpha * parameters.Fbar * potential) \
        - numpy.exp((2 - parameters.alpha) * parameters.Fbar * potential)

V = parameters.omega * parameters.i0 * E(parameters.appliedPotential) / parameters.charge / parameters.faradaysConstant
trenchWidth = 2 * parameters.areaRatio / parameters.perimeterRatio

def plotDeposition(variables, fileprefix, label, figprefix):

    figDeposition = pylab.figure()
    figTheta = pylab.figure()
    figSuppressor = pylab.figure()
    figCupric = pylab.figure()

    for variable in variables:
        filename = fileprefix + str(variable) + '.gz'
        data = dump.read(filename)
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
        pylab.plot(X, depositionRate, label=label % variable)

        pylab.figure(figTheta.number)
        pylab.plot(X, theta, label=label % variable)

        pylab.figure(figSuppressor.number)
        pylab.plot(X, suppressor, label=label % variable)

        pylab.figure(figCupric.number)
        pylab.plot(X, cupric, label=label % variable)
        
    pylab.figure(figDeposition.number)
    pylab.xlabel(r'z (m)', fontsize=16)
    pylab.ylabel(r'v (m / s)', rotation='vertical', fontsize=16)
    pylab.legend()
    pylab.xlim(xmin=-parameters.featureDepth)
    pylab.xlim(xmax=0)
    pylab.savefig(figprefix + 'Depostion.png')

    pylab.figure(figTheta.number)
    pylab.xlabel(r'z (m)', fontsize=16)
    pylab.ylabel(r'$\theta$', rotation='horizontal', fontsize=16)
    pylab.legend()
    pylab.xlim(xmin=-parameters.featureDepth)
    pylab.xlim(xmax=0)
    pylab.ylim(ymax=1)
    pylab.ylim(ymin=0)
    pylab.savefig(figprefix + 'Theta.png')

    pylab.figure(figSuppressor.number)
    pylab.xlabel(r'z (m)', fontsize=16)
    pylab.ylabel(r'$c_{\theta}$', rotation='horizontal', fontsize=16)
    pylab.legend()
    pylab.xlim(xmin=-parameters.featureDepth)
    pylab.xlim(xmax=0)
    pylab.ylim(ymin=0)
    pylab.savefig(figprefix + 'Suppressor.png')

    pylab.figure(figCupric.number)
    pylab.xlabel(r'z (m)', fontsize=16)
    pylab.ylabel(r'$c_{cu}$', rotation='horizontal', fontsize=16)
    pylab.ylim(ymax=parameters.bulkCupric)
    pylab.ylim(ymin=0)
    pylab.legend()
    pylab.xlim(xmin=-parameters.featureDepth)
    pylab.xlim(xmax=0)
    pylab.savefig(figprefix + 'Cupric.png')


plotDeposition((0.01, 5., 25., 50., 75., 100., 125., 150., 300., 1000.),
               'tmp/base-kPlus-',
               r'$k^+=%1.1e$',
               'kPlus')

plotDeposition((1e7, 1.5e7, 2e7, 2.5e7, 3e7),
               'tmp/base-kMinus-',
               r'$k^-=%1.1e$',
               'kMinus')

plotDeposition((0.001, 0.01, 0.02, 0.03, 0.04),
               'tmp/base-deltaRef-',
               r'$\delta_{ref}=%1.1e$',
               'deltaRef')

plotDeposition((0.005, 0.01, 0.02, 0.04, 0.08),
               'tmp/base-bulkSuppressor-',
               r'$c_{\theta}^{\infty}=%1.1e$',
               'bulkSuppressor')

plotDeposition((0.200, 0.250, 0.300),
                'tmp/base-appliedPotential-',
                r'$\eta=-%1.1e$',
                'appliedPotentialDeposition.png')

# plotDeposition((0e-6, 5e-6, 15e-6, 25e-6, 35e-6, 45e-6, 55e-6, 65e-6, 75e-6, 85e-6),
#                'tmp/base-featureDepth-',
#                r'$h=%1.1e$ (m)',
#                'featureDepthDeposition.png')
