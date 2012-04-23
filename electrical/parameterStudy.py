from feature1D import feature

for kPlus in (0.01, 0.1, 1., 5., 10., 25., 50., 75., 100., 125., 150.):
    feature(kPlus=kPlus, filename='tmp/base-kPlus-' + str(kPlus) + '.gz', totalSteps=400)

for kMinus in (0.1e7, 0.5e7, 1e7, 1.5e7, 2e7, 2.5e7, 3e7)
    feature(kMinus=kMinus, filename='tmp/base-kMinus-' + str(kMinus) + '.gz', totalSteps=400)


