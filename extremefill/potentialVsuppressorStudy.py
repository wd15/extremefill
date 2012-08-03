from feature1D import feature
import numpy

appliedPotentials = -10**numpy.linspace(-2, 0, 100)
bulkSuppressors = 10**numpy.linspace(-4, 0, 100)

if __name__ == '__main__':
    for appliedPotential in appliedPotentials:
        for bulkSuppressor in bulkSuppressors:
            filename='tmp/base-appliedPotential-%1.3e-bulkSuppressor-%1.3e' % (appliedPotential, bulkSuppressor)
            print filename
            feature(appliedPotential=appliedPotential, bulkSuppressor=bulkSuppressor, filename=filename, totalSteps=2, sweeps=100, dt=1, tol=1e-4)

