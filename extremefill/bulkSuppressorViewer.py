## Need to import PyTables before importing fipy for some reason.
from depositionViewer import DepositionViewer
import pylab

class BulkSuppressorViewer(DepositionViewer):

    def __init__(self):
        super(BulkSuppressorViewer, self).__init__('bulkSuppressor', (0.005, 0.01, 0.02, 0.04), r'$C_{\text{Supp}}^{\infty}=%1.3f$ $\mole\per\power{\metre}{3}$')

    def _maxSuppressor(self):
        return 0.04

if __name__ == '__main__':
    BulkSuppressorViewer().plot()
