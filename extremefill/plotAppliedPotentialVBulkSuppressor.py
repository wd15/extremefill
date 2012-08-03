import pylab
from fipy import dump, Grid1D
import parameters
import numpy

import matplotlib
#from matplotlib import rc

##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
##rc('text', usetex=True)

matplotlib.rcParams['lines.linewidth'] = 2
#font = {'family' : 'normal',
#        'weight' : 'normal',
#        'size'   : 12}
#matplotlib.rc('font', **font)
fontsize=16
matplotlib.rcParams['legend.fontsize'] = 11
matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize
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

def makeBackGroudPlot(appliedPotential, bulkSuppressor, fP, fM, showPosition=False, axislabel=False, 
align=False, mainaxes=None):
    alpha = 0.5
    xpos = (numpy.log10(bulkSuppressor) + 4.) / 4. *  fM
    ypos = (numpy.log10(-appliedPotential) + 2) / 2. * fP   
    la = pylab.axes((xpos, ypos, 0.08, 0.08), frame_on=True, axisbg='w')
    la.patch.set_alpha(0.5)
    X, depositionRate = getDepositionRate('tmp/base-appliedPotential-%1.3e-bulkSuppressor-%1.3e' % (appliedPotential, bulkSuppressor))
    
    ylim=(0,1e-8)
    xlim=(-parameters.featureDepth, 0)

    if axislabel is True:
        pylab.xlabel(r'$\bm z$', fontsize=fontsize, labelpad=3.5, alpha=1)
        pylab.ylabel(r'$\bm v$', fontsize=fontsize, rotation='horizontal', labelpad=2, alpha=1)
        from matplotlib.patches import FancyArrowPatch
        disx = (xlim[1] - xlim[0]) * 1.4
        disy = (ylim[1] - ylim[0]) * 1.5
        la.add_patch(FancyArrowPatch((xlim[0] - disx * 0.025, ylim[0]),(xlim[0] + disx, ylim[0]),arrowstyle='-|>',mutation_scale=20, lw=2, clip_on=False, facecolor='black'))
        la.add_patch(FancyArrowPatch((xlim[0], ylim[0] - ylim[1] / 20),(xlim[0], ylim[0] + disy),arrowstyle='-|>',mutation_scale=20, lw=2, clip_on=False, facecolor='black'))

    pylab.plot(X, depositionRate, 'k', lw=2, alpha=alpha)
    pylab.setp(la, ylim=ylim, xlim=xlim, xticks=[], yticks=[], alpha=alpha)
    if align:
        pylab.plot((-parameters.featureDepth / 2,), (1e-8 / 2,), 'ks', lw=3) 
    if align or showPosition:
        pylab.axes(mainaxes)
        pylab.loglog((bulkSuppressor,), (-appliedPotential,), 'kx', lw=3, ms=8) 


def plot(filesuffix='.png'):

    from  potentialVsuppressorStudy import appliedPotentials, bulkSuppressors
    figureOfMerits = numpy.zeros((len(appliedPotentials), len(bulkSuppressors), 4), 'd')
    #print bulkSuppressors
    #print appliedPotentials
    #raw_input()
    for i, appliedPotential in enumerate(appliedPotentials):
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

    pylab.xlabel(r'$C_{\text{Supp}}$ $\left(\mole\per\power{\metre}{3}\right)$', fontsize=fontsize)
    pylab.ylabel(r'$-E_{\text{Applied}}$ $\left(\volt\right)$', fontsize=fontsize)

    pylab.text(0.1, 0.04, r'I', fontsize=fontsize)
    pylab.text(0.003, 0.04, r'II', fontsize=fontsize)
    pylab.text(0.0004, 0.04, r'III', fontsize=fontsize)
    pylab.text(0.0002, 0.5, r'IV', fontsize=fontsize)

    for fP, appliedPotential in ((0.885, -2.477e-1), (1.1, -0.02535364)):
        for fM, bulkSuppressor in ((1.8, 2.10490414e-04), (1.04, 1.96304065e-03), (0.925, 2.009e-02), (0.88, 2.05651231e-01)):
            makeBackGroudPlot(appliedPotential, bulkSuppressor, fP, fM, mainaxes=a.ax)

    makeBackGroudPlot(-0.097701, 6.57933225e-03, 0.927, 0.952, axislabel=True, mainaxes=a.ax)
    makeBackGroudPlot(-0.39442061, 6.73415066e-02, 0.875, 0.895, mainaxes=a.ax)
    makeBackGroudPlot(-8.697e-01, 1.024e-03, 0.844, 1.1, mainaxes=a.ax)

    ## experimental data points

    pylab.axes(a.ax)
    pylab.loglog((.08,), (0.275,), 'ko', ms=8, mfc='none', mew=2) 
    #pylab.arrow(.08, 0.275, -0.05, 0, length_includes_head=True, head_width=.01, fill=True)
    pylab.loglog((.08,), (0.25,), 'ko', ms=8, mfc='none', mew=2) 
    #pylab.arrow(.08, 0.25, -0.05, 0, length_includes_head=True)

    from matplotlib.patches import FancyArrowPatch
    a.ax.add_patch(FancyArrowPatch((.08, 0.275),(.03, 0.275),arrowstyle='->',mutation_scale=30, lw=2))
    a.ax.add_patch(FancyArrowPatch((.08, 0.25),(.03, 0.25),arrowstyle='->',mutation_scale=30, lw=2))

    pylab.loglog((.02,), (0.25,), 'ko', ms=8, mew=2) 
    pylab.loglog((.01,), (0.2,), 'ko', ms=8, mew=2) 

    for fs in filesuffix:
        pylab.savefig('appliedPotentialVBulkSuppressor' + fs)

if __name__ == '__main__':
    plot()
    pylab.show()
