import pylab
from fipy import dump, Grid1D
import parameters
import numpy

fontsize=16

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

def makeBackGroudPlot(kPlus, kMinus, xpos, ypos, axisbg=None, axeslabel=False):
    
    a = pylab.axes((xpos, ypos, 0.08, 0.08), frame_on=True, axisbg='w')
    a.patch.set_alpha(0.5)
    X, depositionRate = getDepositionRate('tmp/base-kPlus-%1.2e-kMinus-%1.2e' % (kPlus, kMinus))
    if axeslabel is True:
        pylab.xlabel(r'$z$', fontsize=fontsize, labelpad=2.4)
        pylab.ylabel(r'$v$', fontsize=fontsize, rotation='horizontal', labelpad=-.2)
    pylab.plot(X, depositionRate, 'k', lw=2)
    pylab.setp(a, ylim=(0,1e-8), xlim=(-parameters.featureDepth, 0), xticks=[], yticks=[], alpha=0.8)

def plot(filesuffix='.png'):

    pylab.figure()
    kPluses = 10**numpy.linspace(0, 3, 100)
    kMinuses = 10**numpy.linspace(6, 9, 100)
    figureOfMerits = numpy.zeros((len(kPluses), len(kMinuses), 4), 'd')

    for i, kPlus in enumerate(kPluses):
        print 'kPlus',kPlus
        for j, kMinus in enumerate(kMinuses):
            filename='tmp/base-kPlus-%1.2e-kMinus-%1.2e' % (kPlus, kMinus)
            figureOfMerits[i, j, :] = figureOfMerit(filename)


    #pylab.title(r'$\max \left(\frac{\partial v}{\partial z}\right),\;v > 0$', fontsize=18)
    #a = pylab.contourf(kMinuses, kPluses, figureOfMerits, (-4e-5, -2e-5, 1.9999e-5, 2e-5, 4e-5))
    #a = pylab.contour(kMinuses, kPluses, figureOfMerits, (-1e+10, 0, 1e+10), linestyles='solid', linewidths=4, colors='k')
    a = pylab.contourf(kMinuses, kPluses, figureOfMerits[...,0], (-1e+10, 1e+10), colors=('red',))
    a = pylab.contourf(kMinuses, kPluses, figureOfMerits[...,1], (-1e+10, 0), colors=('orange',))
    a = pylab.contourf(kMinuses, kPluses, figureOfMerits[...,0], (-1e+10, 0), colors=('green',))
    a = pylab.contourf(kMinuses, kPluses, figureOfMerits[...,2], (-1e+10, V / 10000), colors=('blue',))

    a = pylab.contour(kMinuses, kPluses, figureOfMerits[...,3], (-1e+10, 0.5, 1e+10), linestyles='solid', linewidths=4, colors='k')
    #a = pylab.contour(kMinuses, kPluses, figureOfMerits[...,1], (-1e+10, 0, 1e+10), linestyles='solid', linewidths=4, colors='g')
    #a = pylab.contour(kMinuses, kPluses, figureOfMerits[...,2], (-1e+10, V / 10000, 1e+10), linestyles='solid', linewidths=4, colors='y')
    a.ax.set_xscale('log')
    a.ax.set_yscale('log')
    pylab.xticks((10**6, 10**7, 10**8, 10**9), fontsize=fontsize)
    pylab.yticks((10**0, 10**1, 10**2, 10**3), fontsize=fontsize)
    pylab.xlabel(r'$k^-$ (1 / m)', fontsize=fontsize)
    pylab.ylabel(r'$k^+$ (m$^3$ / mol s)', fontsize=fontsize)

    pylab.text(2 * 10**6, 7 * 10**2, r'I', fontsize=fontsize)
    pylab.text(3 * 10**7, 7 * 10**2, r'II', fontsize=fontsize)
    pylab.text(6 * 10**8, 7 * 10**2, r'III', fontsize=fontsize)
    pylab.text(6 * 10**8, 7 * 10**1, r'IV', fontsize=fontsize)



    for fP, kPlus, paxeslabel in ((1.1,3.51e+00, False), (1.0, 9.33e+00, False), (0.9, 3.51e+01, False), (0.89, 9.33e+01, True
    ), (0.85, 3.51e+02, False)):
        for fM, kMinus, maxeslabel in ((1.4, 2.48e+06, False), (1.1, 7.05e+6, False), (0.95, 2.48e+07, True), (0.91, 7.05e+7, False), (0.88, 2.48e+08, False)):
            xpos = (numpy.log10(kMinus) - 6.) / 3. *  fM
            ypos = numpy.log10(kPlus) / 3. * fP        
            makeBackGroudPlot(kPlus, kMinus, xpos, ypos, axisbg=None, axeslabel=paxeslabel and maxeslabel)

    pylab.savefig('kPlusVkMinus' + filesuffix)

if __name__ == '__main__':
    plot()
    pylab.show()
