from extremefill.depositionViewer import DepositionViewer

class AppliedPotentialViewer(DepositionViewer):

    def __init__(self, datafile='data.h5'):
        super(AppliedPotentialViewer, self).__init__('appliedPotential', (-0.200, -0.250, -0.300), r'$E_{\text{Applied}}=%1.2f$ $\volt$', datafile=datafile)

if __name__ == '__main__':
    AppliedPotentialViewer().plot()
