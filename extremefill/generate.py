import tables
from dicttable import DictTable
from main import run

def generateDataSet(parameter=None, values=None, datafile='data.h5'):

    h5data = DictTable(datafile)
    dataset = []
    for value in values:
        key = parameter + value.replace('-', 'm').replace('+', 'p') ## replace due to PyTables natural naming scheme
        if h5data.haskey(key):
            dataset.append(h5data[key])
        else:
            print 'generating data for ' + parameter + ' = ' + value
            print
            data = run(totalSteps=10, **{parameter : float(value)})
            h5data[key] = data
            dataset.append(data)

    return dataset
        
##for kMinus in (6e7, 1e8, 5e8):
##    feature(kMinus=kMinus, filename='tmp/base-kMinus-' + str(kMinus) + '.gz', totalSteps=400)

# for deltaRef in (0.001, 0.01, 0.02, 0.03, 0.04):
#     feature(deltaRef=deltaRef, filename='tmp/base-deltaRef-' + str(deltaRef) + '.gz', totalSteps=400)

# for bulkSuppressor in (0.005, 0.01, 0.02, 0.04, 0.08):
#     feature(bulkSuppressor=bulkSuppressor, filename='tmp/base-bulkSuppressor-' + str(bulkSuppressor) + '.gz', totalSteps=400)

# for appliedPotential in (-0.200, -0.250, -0.300):
#     feature(appliedPotential=appliedPotential, filename='tmp/base-appliedPotential-' + str(appliedPotential) + '.gz', totalSteps=30, dt=1.e+20)

# for featureDepth in (0e-6, 5e-6, 15e-6, 25e-6, 35e-6, 45e-6, 55e-6, 65e-6, 75e-6, 85e-6):
#     feature(featureDepth=featureDepth, filename='tmp/base-featureDepth-' + str(featureDepth) + '.gz', totalSteps=400) 



