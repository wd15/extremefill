## Need to import PyTables before importing fipy for some reason.
import pylab
import numpy

from depositionViewer import DepositionViewer

class AppliedPotentialViewer(DepositionViewer):

    def __init__(self):
        super(AppliedPotentialViewer, self).__init__('appliedPotential', (-0.200, -0.250, -0.300), r'$E_{\text{Applied}}=%1.2f$ $\volt$')

if __name__ == '__main__':
    AppliedPotentialViewer().plot()
