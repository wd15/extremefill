## Need to import PyTables before importing fipy for some reason.
import pylab

from depositionViewer import DepositionViewer

class DeltaRefViewer(DepositionViewer):
    def __init__(self):
        super(DeltaRefViewer, self).__init__('deltaRef', (0.001, 0.005, 0.01, 0.02, 0.03, 0.04),  r'$L=%1.3f$ $\metre$')

    def _legend(self, ax):
        if ax.colNum == 1 and ax.rowNum == 0:
            self.legend = pylab.legend(loc='upper right')

if __name__ == '__main__':
    DeltaRefViewer().plot()
