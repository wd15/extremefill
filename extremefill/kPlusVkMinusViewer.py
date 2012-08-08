import pylab
import numpy


from contourViewer import ContourViewer

class KPlusVKMinusViewer(ContourViewer):
    fontsize=16

    def plot(self, filesuffix=('.png',)):

        pylab.figure()

        kPluses = 10**numpy.linspace(0, 3, 100)
        kMinuses = 10**numpy.linspace(6, 9, 100)
        figureOfMerits = numpy.zeros((len(kPluses), len(kMinuses), 4), 'd')
        for i, kPlus in enumerate(kPluses):
            for j, kMinus in enumerate(kMinuses):
                figureOfMerits[i, j, :] = self.figureOfMerit(self.generateData({'kPlus' : kPlus, 'kMinus' : kMinus}))
        data = self.generateData({'kPlus' : kPluses[0], 'kMinus' : kMinuses[0]})

        self._contourf(kMinuses, kPluses, figureOfMerits, data)

        pylab.xticks((10**6, 10**7, 10**8, 10**9), fontsize=self.fontsize)
        pylab.yticks((10**0, 10**1, 10**2, 10**3), fontsize=self.fontsize)
        pylab.xlabel(r'$k^-$ $\left(1\per\metre\right)$', fontsize=self.fontsize)
        pylab.ylabel(r'$k^+$ $\left(\power{\metre}{3}\per\mole\cdot\second\right)$', fontsize=self.fontsize)

        pylab.text(2 * 10**6, 7 * 10**2, r'I', fontsize=self.fontsize)
        pylab.text(3 * 10**7, 7 * 10**2, r'II', fontsize=self.fontsize)
        pylab.text(6 * 10**8, 7 * 10**2, r'III', fontsize=self.fontsize)
        pylab.text(6 * 10**8, 7 * 10**1, r'IV', fontsize=self.fontsize)

        for fP, kPlus, paxeslabel in ((1.143,3.51e+00, False), (0.975, 9.33e+00, False), (0.916, 3.51e+01, False), (0.89, 9.33e+01, False), (0.87, 3.51e+02, True)):
            for fM, kMinus, maxeslabel in ((1.4, 2.48e+06, False), (1.07, 7.05e+6, False), (0.96, 2.48e+07, False), (0.91, 7.05e+7, False), (0.88, 2.48e+08, True)):
                xpos = (numpy.log10(kMinus) - 6.) / 3. *  fM
                ypos = numpy.log10(kPlus) / 3. * fP        
                self.makeBackGroundPlot({'kPlus' : kPlus, 'kMinus' : kMinus}, xpos, ypos, axeslabel=paxeslabel and maxeslabel)

        for fs in filesuffix:
            pylab.savefig('kPlusVkMinus' + fs)

if __name__ == '__main__':
    KPlusVKMinusViewer().plot()

