## Need to import PyTables before importing fipy for some reason.
import tables

import pylab
from fipy import Grid1D
import numpy
import matplotlib

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

from plot import DepositionViewer
from generate import generateDataSet

class KPlusViewer(DepositionViewer):

    def __init__(self):
        parameter = 'kPlus'
        super(KPlusViewer, self).__init__(parameter, (1e-2, 5e0, 2.5e1, 5e1, 1e2, 1e3), r'$k^+=%4.2f$ $\power{\metre}{3}\per\mole\cdot\second$')
        self.inlayDataSet = generateDataSet(parameter=parameter, values=['%1.2e' % kPlus for kPlus in 10**numpy.linspace(0, 3, 100)])

    def plot(self):
        super(KPlusViewer, self).plot(legend=3, lfs=8, subplot=True, filesuffix='.png', replaceString='.00$', inlayDataSet=self.inlayDataSet) 

    # plotDeposition((1e7, 1.5e7, 2e7, 2.5e7, 3e7),
    #                'tmp/base-kMinus-',
    #                r'$k^-=%1.1f\times 10^{%i}$ $1\per\metre$',
    #                'kMinus', filesuffix=filesuffix)

    # plotDeposition((0.001, 0.005, 0.01, 0.02, 0.03, 0.04),
    #                'tmp/base-deltaRef-',
    #                r'$L=%1.3f$ $\metre$',
    #                'deltaRef',
    #                legend=2, 
    #                loc='upper right', filesuffix=filesuffix)

    # plotDeposition((0.005, 0.01, 0.02, 0.04),
    #                'tmp/base-bulkSuppressor-',
    #                r'$C_{\text{Supp}}^{\infty}=%1.3f$ $\mole\per\power{\metre}{3}$',
    #                'bulkSuppressor',
    #                maxSuppressor=0.04, filesuffix=filesuffix)

    # plotDeposition((-0.200, -0.250, -0.300),
    #                'tmp/base-appliedPotential-',
    #                r'$E_{\text{Applied}}=%1.2f$ $\volt$',
    #                'appliedPotential', filesuffix=filesuffix)

    # plotDeposition(numpy.array((15e-6, 25e-6, 35e-6, 45e-6, 55e-6, 65e-6, 75e-6, 85e-6))[::-1],
    #                'tmp/base-featureDepth-',
    #                r'$h=%1.0f$ $\micro\metre$',
    #                'featureDepth',
    #                mulFactor=1000000,
    #                legend=3,
    #                filesuffix=filesuffix,
    #                xticks=(-80, -60, -40, -20, 0),
    #                colors = numpy.array(['b', 'g', 'r', 'c', 'm', 'y', 'k', '#663300'])[::-1])

if __name__ == '__main__':
    KPlusViewer().plot()
