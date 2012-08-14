## Need to import PyTables before importing fipy for some reason.
# import matplotlib
# from matplotlib import rc
# matplotlib.use('Agg')
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
import pylab
pylab.ioff()
pylab.hold(True)
import numpy

from extremefill.depositionViewer import DepositionViewer

class FeatureDepthViewer(DepositionViewer):

    def __init__(self, datafile='data.h5'):
        super(FeatureDepthViewer, self).__init__('featureDepth', numpy.array((15e-6, 25e-6, 35e-6, 45e-6, 55e-6, 65e-6, 75e-6, 85e-6))[::-1], r'$h=%1.0f$ $\micro\metre$', datafile=datafile)

    def _legend(self, ax):
        if ax.colNum == 0 and ax.rowNum == 1:
            self.legend =  pylab.legend(loc='upper left')

    def scale(self, variable):
        return variable * 1000000

    def _xticks(self):
        return (-80, -60, -40, -20, 0)

    def _colors(self):
        return numpy.array(['b', 'g', 'r', 'c', 'm', 'y', 'k', '#663300'])[::-1]

if __name__ == '__main__':
    FeatureDepthViewer().plot()
