from extremefill.depositionViewer import DepositionViewer

class BulkSuppressorViewer(DepositionViewer):

    def __init__(self, datafile='data.h5'):
        super(BulkSuppressorViewer, self).__init__('bulkSuppressor', (0.005, 0.01, 0.02, 0.04), r'$C_{\text{Supp}}^{\infty}=%1.3f$ $\mole\per\power{\metre}{3}$', datafile=datafile)

    def _maxSuppressor(self):
        return 0.04

if __name__ == '__main__':
    BulkSuppressorViewer().plot()
