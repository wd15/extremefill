import pylab
from fipy import dump, Grid1D
import parameters
import numpy
import matplotlib

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
        return max(derivative[mask])
    else:
        return +1e-10

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

kPluses = 10**numpy.linspace(0, 3, 100)
kMinuses = 10**numpy.linspace(6, 9, 100)
figureOfMerits = numpy.zeros((len(kPluses), len(kMinuses)), 'd')
mx = 4e-5
for i, kPlus in enumerate(kPluses):
    print 'kPlus',kPlus
    for j, kMinus in enumerate(kMinuses):
        filename='tmp/base-kPlus-%1.2e-kMinus-%1.2e' % (kPlus, kMinus)
        fom = max(figureOfMerit(filename), -mx)
        figureOfMerits[i, j] = min(fom, mx)

#pylab.title(r'$\max \left(\frac{\partial v}{\partial z}\right),\;v > 0$', fontsize=18)
a = pylab.contourf(kMinuses, kPluses, figureOfMerits, (-4e-5, -2e-5, 1.9999e-5, 2e-5, 4e-5))
a = pylab.contour(kMinuses, kPluses, figureOfMerits, (-1e+10, 0, 1e+10), linestyles='solid', linewidths=4, colors='k')
a.ax.set_xscale('log')
a.ax.set_yscale('log')
pylab.xlabel(r'$k^-$ (1 / m)', fontsize=14)
pylab.ylabel(r'$k^+$ (m$^3$ / mol s)', fontsize=14)

def makeBackGroudPlot(kPlus, kMinus, xpos, ypos, axisbg=None):

    a = pylab.axes((xpos, ypos, 0.08, 0.08))
    X, depositionRate = getDepositionRate('tmp/base-kPlus-%1.2e-kMinus-%1.2e' % (kPlus, kMinus))
    pylab.plot(X, depositionRate, 'k', lw=2)
    pylab.setp(a, ylim=(0,1e-8), xlim=(-parameters.featureDepth, 0), xticks=[], yticks=[])

for fP, kPlus in ((1.1,3.51e+00), (1.0, 9.33e+00), (0.9, 3.51e+01), (0.89, 9.33e+01), (0.85, 3.51e+02)):
    for fM, kMinus in ((1.4, 2.48e+06), (1.1, 7.05e+6), (0.95, 2.48e+07), (0.91, 7.05e+7), (0.88, 2.48e+08)):
        xpos = (numpy.log10(kMinus) - 6.) / 3. *  fM
        ypos = numpy.log10(kPlus) / 3. * fP        
        makeBackGroudPlot(kPlus, kMinus, xpos, ypos, axisbg=None)

pylab.savefig('kPlusVkMinus.png')
pylab.show()
