import pylab
from fipy import dump, Grid1D
import parameters
import numpy
import matplotlib

matplotlib.rcParams['lines.linewidth'] = 2
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}
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
fig1 = pylab.figure()
fig2 = pylab.figure()

trenchWidth = 2 * parameters.areaRatio / parameters.perimeterRatio
foms = []
kPluses = (0.01, 0.1, 1., 5., 10., 25., 50., 75., 100., 125., 150., 300., 600., 1000.)

for kPlus in kPluses:
    K = parameters.gamma * kPlus * parameters.delta / parameters.diffusionSuppressor
    filename='tmp/base-kPlus-' + str(kPlus) + '.gz'
    data = dump.read(filename)
    theta = data['theta'][:ID + 1]
    cupric = data['cupric'][:ID + 1]
    potential = data['potential'][:ID + 1]
    cupric = data['cupric'][:ID + 1]
    I0 = parameters.i0 + parameters.i1 * theta
    cbar = cupric / parameters.bulkCupric
    current = cbar * I0 * E(-potential)
    depositionRate = parameters.omega * current / parameters.charge / parameters.faradaysConstant   
    figureOfMerit = 1 - 2 * (parameters.featureDepth + X) * depositionRate / depositionRate[0] / trenchWidth
    foms += [min(figureOfMerit)]
    pylab.figure(fig1.number)
    pylab.plot(X / parameters.featureDepth, depositionRate / V, label=r'$k^+=%1.1e$' % kPlus)

    pylab.figure(fig2.number)
    pylab.plot(X / parameters.featureDepth, figureOfMerit, label=r'$k^+=%1.1e$' % kPlus)
    
pylab.figure(fig1.number)
pylab.xlabel(r'$z / h$', fontsize=16)
pylab.ylabel(r'$v / v_0$', rotation='horizontal', fontsize=16)
pylab.legend()
pylab.xlim(xmax=0)
pylab.savefig('kPlusDepositionRate.png')

pylab.figure(fig2.number)
pylab.xlabel(r'$z / h$', fontsize=16)
pylab.ylabel(r'$1 - \frac{2 \left(h + z\right) v}{w v_0}$', rotation='vertical', fontsize=16)
pylab.legend(loc='lower left')
pylab.xlim(xmax=0)
pylab.savefig('kPlusFigureOfMerit.png')

pylab.figure()
pylab.semilogx(kPluses, foms)
pylab.xlabel(r'$k^+$ (m$^3$ / mol s)', fontsize=16)
pylab.ylabel(r'$\min\left(1 - \frac{2 \left(h + z\right) v}{w v_0}\right)$', rotation='vertical', fontsize=16)
pylab.savefig('kPlusVFigureOfMerit.png')
