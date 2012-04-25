from feature1D import feature

#for kPlus in (300., 600., 1000.):
#    feature(kPlus=kPlus, filename='tmp/base-kPlus-' + str(kPlus) + '.gz', totalSteps=400)
#for kPlus in (0.01, 0.1, 1., 5., 10., 25., 50., 75., 100., 125., 150.):
#    feature(kPlus=kPlus, filename='tmp/base-kPlus-' + str(kPlus) + '.gz', totalSteps=400)

##for kMinus in (6e7, 1e8, 5e8):
##    feature(kMinus=kMinus, filename='tmp/base-kMinus-' + str(kMinus) + '.gz', totalSteps=400)

# for deltaRef in (0.001, 0.01, 0.02, 0.03, 0.04):
#     feature(deltaRef=deltaRef, filename='tmp/base-deltaRef-' + str(deltaRef) + '.gz', totalSteps=400)

# for bulkSuppressor in (0.005, 0.01, 0.02, 0.04, 0.08):
#     feature(bulkSuppressor=bulkSuppressor, filename='tmp/base-bulkSuppressor-' + str(bulkSuppressor) + '.gz', totalSteps=400)

for appliedPotential in (-0.200, -0.250, -0.300):
    feature(appliedPotential=appliedPotential, filename='tmp/base-appliedPotential-' + str(appliedPotential) + '.gz', totalSteps=30, dt=1.e+20)

# for featureDepth in (0e-6, 5e-6, 15e-6, 25e-6, 35e-6, 45e-6, 55e-6, 65e-6, 75e-6, 85e-6):
#     feature(featureDepth=featureDepth, filename='tmp/base-featureDepth-' + str(featureDepth) + '.gz', totalSteps=400) 



