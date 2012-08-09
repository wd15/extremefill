## Need to import PyTables before importing fipy for some reason.
import tables

from fipy import Grid1D
import numpy

from extremefill.dicttable import DictTable
from extremefill.simulation import Simulation

class Viewer(object):
    def __init__(self, datafile='data.h5'):
        self.datafile = datafile


    def generateData(self, parameters):
        h5data = DictTable(self.datafile)
        h5key = ''
        for key in parameters.keys():
            h5key += key + ('%1.2e' % parameters[key]).replace('-', 'm').replace('+', 'p') ## replace due to PyTables natural naming scheme
        if h5data.haskey(h5key):
            data = h5data[h5key]
        else:
            print 'generating data for ' + h5key
            parameters['totalSteps'] = 10
            simulation = Simulation(**parameters)
            simulation.run()
            data = simulation.parameters
            h5data[h5key] = data
        return data

    def getX(self, data):
        L = data['delta'] + data['featureDepth']
        N = 1000
        dx = L / N 
        mesh = Grid1D(nx=N, dx=dx) - [[data['featureDepth']]]
        x = 0
        ID = numpy.argmin(abs(mesh.x - x))
        if mesh.x[ID] < x:
            ID = ID + 1
        X = mesh.x[:ID + 1].value
        return X, ID, dx

    def E(self, potential, data):
        return numpy.exp(-data['alpha'] * data['Fbar'] * potential) \
            - numpy.exp((2 - data['alpha']) * data['Fbar'] * potential)

    def getDepositionRate(self, data):
        X, ID, dx = self.getX(data)
        theta = data['theta'][:ID + 1]
        potential = data['potential'][:ID + 1]
        cupric = data['cupric'][:ID + 1]
        suppressor = data['suppressor'][:ID + 1]
        I0 = data['i0'] + data['i1'] * theta
        cbar = cupric / data['bulkCupric']
        current = cbar * I0 * self.E(-potential, data)
        depositionRate = data['omega'] * current / data['charge'] / data['faradaysConstant']
        return depositionRate, theta, potential, cupric, suppressor, X, dx
