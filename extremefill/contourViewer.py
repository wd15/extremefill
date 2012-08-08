import pylab
from extremefill.viewer import Viewer

class ContourViewer(Viewer):
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
    
    def makeBackGroundPlot(self, parameters, xpos, ypos, axeslabel=False):
        alpha = 0.5
        a = pylab.axes((xpos, ypos, 0.08, 0.08), frame_on=True, axisbg='w', alpha=alpha)
        a.patch.set_alpha(alpha)
        data = self.generateData(parameters)
        depositionRate, theta, potential, cupric, suppressor, X, dx = self.getDepositionRate(data)

        xlim = (-data['featureDepth'], 0)
        ylim = (0,1e-8)

        if axeslabel is True:
            pylab.xlabel(r'$\bm z$', fontsize=self.fontsize, labelpad=3.5, alpha=1)
            pylab.ylabel(r'$\bm v$', fontsize=self.fontsize, rotation='horizontal', labelpad=2, alpha=1)

            from matplotlib.patches import FancyArrowPatch
            disx = (xlim[1] - xlim[0]) * 1.4
            disy = (ylim[1] - ylim[0]) * 1.5
            a.add_patch(FancyArrowPatch((xlim[0] - disx * 0.025, ylim[0]),(xlim[0] + disx, ylim[0]),arrowstyle='-|>',mutation_scale=20, lw=2, clip_on=False, facecolor='black'))
            a.add_patch(FancyArrowPatch((xlim[0], ylim[0] - ylim[1] / 20),(xlim[0], ylim[0] + disy),arrowstyle='-|>',mutation_scale=20, lw=2, clip_on=False, facecolor='black'))

        pylab.plot(X, depositionRate, 'k', lw=2, alpha=alpha)
        pylab.setp(a, ylim=ylim, xlim=xlim, xticks=[], yticks=[], alpha=alpha)

        a.patch.set_alpha(alpha)

    def _contourf(self, var1, var2, figureOfMerits, data):
        a = pylab.contourf(var1, var2, figureOfMerits[...,0], (-1e+10, 1e+10), colors=('red',))
        a = pylab.contourf(var1, var2, figureOfMerits[...,1], (-1e+10, 0), colors=('orange',))
        a = pylab.contourf(var1, var2, figureOfMerits[...,0], (-1e+10, 0), colors=('green',))
        a = pylab.contourf(var1, var2, figureOfMerits[...,2], (-1e+10, self.V(data) / 10000), colors=('blue',))
        a = pylab.contour(var1, var2, figureOfMerits[...,3], (-1e+10, 0.5, 1e+10), linestyles='solid', linewidths=4, colors='k')
        a.ax.set_xscale('log')
        a.ax.set_yscale('log')
        return a
