import pylab
import numpy

from extremefill.contourViewer import ContourViewer
import matplotlib
fontsize=16
matplotlib.rcParams['legend.fontsize'] = 11
matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize

class AppliedPotentialVbulkSuppressorViewer(ContourViewer):
    fontsize = fontsize

    def V(self, data):
        return data['omega'] * data['i0'] * self.E(data['appliedPotential'], data) / data['charge'] / data['faradaysConstant']

    def figureOfMerit(self, data):
        depositionRate, theta, potential, cupric, suppressor, X, dx = self.getDepositionRate(data)
        derivative = (depositionRate[2:] - depositionRate[:-2]) / dx / 2
        depositionRate = depositionRate[1:-1]
        mask = (depositionRate > (self.V(data) / 10000)) & (depositionRate > (max(depositionRate) / 100))

        if mask.any():
            if mask.all():
                localized = 0
            else:
                localized = 1
            return max(derivative[mask]), min(derivative[mask]), max(depositionRate), localized
        else:
            return -1e-10, -1e-10, max(depositionRate), 0

    def _makeBackGroundPlot(self, appliedPotential, bulkSuppressor, fP, fM, axeslabel=False):
        xpos = (numpy.log10(bulkSuppressor) + 4.) / 4. *  fM
        ypos = (numpy.log10(-appliedPotential) + 2) / 2. * fP   
        self.makeBackGroundPlot({'appliedPotential' : appliedPotential, 'bulkSuppressor' : bulkSuppressor}, xpos, ypos, axeslabel=axeslabel)
        
    def plot(self, filesuffix='.png'):

        appliedPotentials = -10**numpy.linspace(-2, 0, 100)
        bulkSuppressors = 10**numpy.linspace(-4, 0, 100)
        figureOfMerits = numpy.zeros((len(appliedPotentials), len(bulkSuppressors), 4), 'd')

        for i, appliedPotential in enumerate(appliedPotentials):
            for j, bulkSuppressor in enumerate(bulkSuppressors):
                figureOfMerits[i, j, :] = self.figureOfMerit(self.generateData({'appliedPotential' : appliedPotential, 'bulkSuppressor' : bulkSuppressor}))
        data = self.generateData({'appliedPotential' : appliedPotentials[0], 'bulkSuppressor' : bulkSuppressors[0]})

        a = self._contourf(bulkSuppressors, -appliedPotentials, figureOfMerits, data)

        pylab.xlabel(r'$C_{\text{Supp}}$ $\left(\mole\per\power{\metre}{3}\right)$', fontsize=fontsize)
        pylab.ylabel(r'$-E_{\text{Applied}}$ $\left(\volt\right)$', fontsize=fontsize)

        pylab.text(0.1, 0.04, r'I', fontsize=fontsize)
        pylab.text(0.003, 0.04, r'II', fontsize=fontsize)
        pylab.text(0.0004, 0.04, r'III', fontsize=fontsize)
        pylab.text(0.0002, 0.5, r'IV', fontsize=fontsize)

        for fP, appliedPotential in ((0.885, -2.477e-1), (1.1, -0.02535364)):
            for fM, bulkSuppressor in ((1.8, 2.10490414e-04), (1.04, 1.96304065e-03), (0.925, 2.009e-02), (0.88, 2.05651231e-01)):
                self._makeBackGroundPlot(appliedPotential, bulkSuppressor, fP, fM)

        self._makeBackGroundPlot(-0.097701, 6.57933225e-03, 0.927, 0.952, axeslabel=True)
        self._makeBackGroundPlot(-0.39442061, 6.73415066e-02, 0.875, 0.895)
        self._makeBackGroundPlot(-8.697e-01, 1.024e-03, 0.844, 1.1)

        ## experimental data points

        pylab.axes(a.ax)
        pylab.loglog((.08,), (0.275,), 'ko', ms=8, mfc='none', mew=2) 
        pylab.loglog((.08,), (0.25,), 'ko', ms=8, mfc='none', mew=2) 

        from matplotlib.patches import FancyArrowPatch
        a.ax.add_patch(FancyArrowPatch((.08, 0.275),(.03, 0.275),arrowstyle='->',mutation_scale=30, lw=2))
        a.ax.add_patch(FancyArrowPatch((.08, 0.25),(.03, 0.25),arrowstyle='->',mutation_scale=30, lw=2))

        pylab.loglog((.02,), (0.25,), 'ko', ms=8, mew=2) 
        pylab.loglog((.01,), (0.2,), 'ko', ms=8, mew=2) 

        for fs in filesuffix:
            pylab.savefig('appliedPotentialVbulkSuppressor' + fs)

if __name__ == '__main__':
    AppliedPotentialVbulkSuppressorViewer().plot()

