## Need to import PyTables before importing fipy for some reason.
import tables

import pylab
from fipy import Grid1D
import numpy
import matplotlib

matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['legend.labelspacing'] = 0.1
matplotlib.rcParams['figure.subplot.wspace'] = 0.3
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
from dicttable import DictTable
from main import run

class DepositionViewer(object):
    def __init__(self, parameter, values, label, lfs=10):
        self.parameter = parameter
        self.dataset = self.generateDataSet(parameter=self.parameter, values=['%1.2e' % kPlus for kPlus in values])
        self.label = label
        self.lfs = 8

    def generateDataSet(self, parameter=None, values=None, datafile='data.h5'):

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
        return X, ID

    def E(self, potential, data):
        return numpy.exp(-data['alpha'] * data['Fbar'] * potential) \
            - numpy.exp((2 - data['alpha']) * data['Fbar'] * potential)

    def _legend(self, ax):
        if ax.colNum == 0 and ax.rowNum == 0:
            self.legend =  pylab.legend(loc='upper right')
        
    def plot(self, mulFactor=1, filesuffix='.png', xticks=(-50, -40, -30, -20, -10, 0), colors=None):        
        
        pylab.figure()
        figDeposition = pylab.subplot(221)
        figTheta = pylab.subplot(222)
        figSuppressor = pylab.subplot(223)
        figCupric = pylab.subplot(224)
        maxFeatureDepth = 0
        print '**************'
        print 'varying ' + self.parameter

        xlabel = r'$z$ $\left(\micro\metre\right)$'
        scaleFactor = 1000000

        if colors is None:
            colors = [None] * len(self.dataset)

        maxSuppressor = self._maxSuppressor()


        for color, data in zip(colors, self.dataset):
            variable = data[self.parameter]
            featureDepth = data['featureDepth']
            X, ID = self.getX(data)
            maxFeatureDepth = max(featureDepth, maxFeatureDepth)
            theta = data['theta'][:ID + 1]
            potential = data['potential'][:ID + 1]
            cupric = data['cupric'][:ID + 1]
            suppressor = data['suppressor'][:ID + 1]
            I0 = data['i0'] + data['i1'] * theta
            cbar = cupric / data['bulkCupric']
            current = cbar * I0 * self.E(-potential, data)
            depositionRate = data['omega'] * current / data['charge'] / data['faradaysConstant']

            if 'times' in self.label:
                def me(n):
                    s = '%1.1e' % n
                    m, e = s.split('e')
                    return float(m), int(e)
                Label = self.label % me(variable * mulFactor)
            else:
                Label = self.label % (variable * mulFactor)
                Label = self.replaceString(Label)

            kwargs = {'label' : Label}
            if color is not None:
                kwargs['color'] = color

            pylab.axes(figDeposition)        
            pylab.plot(X * scaleFactor, depositionRate, **kwargs)

            pylab.axes(figTheta)
            pylab.plot(X * scaleFactor, theta, **kwargs)

            pylab.axes(figSuppressor)
            pylab.plot(X * scaleFactor, suppressor, **kwargs)

            pylab.axes(figCupric)
            pylab.plot(X * scaleFactor, cupric, **kwargs)

            print '------------'
            print self.parameter + ' value: ' + str(data[self.parameter])
            print 'voltage drop',-data['appliedPotential'] - potential[ID]

        ax = pylab.axes(figDeposition)
        self._legend(ax)
        pylab.ylabel(r'$v$ $\left(\metre\per\second\right)$', rotation='vertical', fontsize=14)
        pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
        pylab.xlim(xmax=0)
        pylab.ylim(ymax=8e-9)
        pylab.ylim(ymin=0.0)
        pylab.xticks(xticks)
        pylab.title('(a)')

        ax = pylab.axes(figTheta)
        self._legend(ax)
        pylab.ylabel(r'$\theta$', rotation='vertical')
        pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
        pylab.xlim(xmax=0)
        pylab.ylim(ymax=1* 1.05)
        pylab.ylim(ymin=0)
        pylab.xticks(xticks)
        pylab.title('(b)')

        ax = pylab.axes(figSuppressor)
        self._legend(ax)
        pylab.xlabel(xlabel, fontsize=14)
        pylab.ylabel(r'$C_{\text{Supp}}$ $\left(\mole\per\power{\metre}{3}\right)$', rotation='vertical', fontsize=14)
        pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
        pylab.xlim(xmax=0)
        pylab.ylim(ymin=0)
        pylab.ylim(ymax=maxSuppressor* 1.05)
        pylab.xticks(xticks)
        pylab.title('(c)')

        ax = pylab.axes(figCupric)
        self._legend(ax)
        pylab.xlabel(xlabel, fontsize=14)
        pylab.ylabel(r'$C_{\text{Cu}}$ $\left(\mole\per\power{\metre}{3}\right)$', rotation='vertical', fontsize=14, labelpad=-3)
        pylab.ylim(ymax=data['bulkCupric'] * 1.05)
        pylab.ylim(ymin=0)
        pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
        pylab.xlim(xmax=0)
        pylab.xticks(xticks)
        pylab.title('(d)')

        self.legend.labelspacing = 0
        self.legend.columnspacing = 0.1
        for t in self.legend.texts:
            t.set_size(self.lfs)

        self.subplot(figDeposition)


        if type(filesuffix) is str:
            filesuffix = (filesuffix,)

        for fs in filesuffix:
            print self.parameter  + fs
            pylab.savefig(self.parameter + fs)

    def subplot(self, fig):
        pass

    def replaceString(self, Label):
        return Label

    def _maxSuppressor(self):
        return self.dataset[0]['bulkSuppressor']        
