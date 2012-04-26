from feature1D import feature
import numpy

for kPlus in 10**numpy.linspace(0, 3, 100):
    for kMinus in 10**numpy.linspace(6, 9, 100):
        filename='tmp/base-kPlus-%1.2e-kMinus-%1.2e' % (kPlus, kMinus)
        print filename
        feature(kPlus=kPlus, kMinus=kMinus, filename=filename, totalSteps=1, sweeps=100, dt=1e20, tol=1e-4)

