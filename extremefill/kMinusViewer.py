from extremefill.depositionViewer import DepositionViewer

class KMinusViewer(DepositionViewer):
    def __init__(self, datafile='data.h5'):
        super(KMinusViewer, self).__init__('kMinus', (1e7, 1.5e7, 2e7, 2.5e7, 3e7), r'$k^-=%1.1f\times 10^{%i}$ $1\per\metre$', datafile=datafile)

if __name__ == '__main__':
    KMinusViewer().plot()
