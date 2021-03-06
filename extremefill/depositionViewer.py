## Need to import PyTables before importing fipy for some reason.
import matplotlib
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['legend.labelspacing'] = 0.1
matplotlib.rcParams['figure.subplot.wspace'] = 0.3
matplotlib.rcParams['figure.subplot.hspace'] = 0.3

# from matplotlib import rc
# matplotlib.use('Agg')
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
import pylab
# pylab.ioff()
# pylab.hold(True)
from extremefill.viewer import Viewer

class DepositionViewer(Viewer):
    def __init__(self, parameter, values, label, lfs=10, datafile='data.h5'):
        super(DepositionViewer, self).__init__(datafile=datafile)
        self.parameter = parameter
        self.dataset = []
        for value in values:
            self.dataset.append(self.generateData({parameter : value}))
        self.label = label
        self.lfs = lfs

    def _legend(self, ax):
        if ax.colNum == 0 and ax.rowNum == 0:
            self.legend =  pylab.legend(loc='upper right')
        
    def plot(self, mulFactor=1, filesuffix='.png', colors=None, deposition_only=False):        
        
        if deposition_only:
            IDs = (221, 221, 221, 221)
        else:
            IDs = (221, 222, 223, 224)

        figDeposition = pylab.subplot(IDs[0])
        figTheta = pylab.subplot(IDs[1])
        figSuppressor = pylab.subplot(IDs[2])
        figCupric = pylab.subplot(IDs[3])


        maxFeatureDepth = 0

        xlabel = r'$z$ $\left(\micro\metre\right)$'
        scaleFactor = 1000000

        maxSuppressor = self._maxSuppressor()

        for color, data in zip(self._colors(), self.dataset):
            variable = data[self.parameter]
            featureDepth = data['featureDepth']
            maxFeatureDepth = max(featureDepth, maxFeatureDepth)
            
            depositionRate, theta, potential, cupric, suppressor, X, dx = self.getDepositionRate(data)
            
            if 'times' in self.label:
                def me(n):
                    s = '%1.1e' % n
                    m, e = s.split('e')
                    return float(m), int(e)
                Label = self.label % me(self.scale(variable))
            else:
                Label = self.label % (self.scale(variable))
                Label = self.replaceString(Label)

            kwargs = {'label' : Label}
            if color is not None:
                kwargs['color'] = color

            pylab.axes(figDeposition)        
            pylab.plot(X * scaleFactor, depositionRate, **kwargs)

            if not deposition_only:
                pylab.axes(figTheta)

                pylab.plot(X * scaleFactor, theta, **kwargs)

                pylab.axes(figSuppressor)
                pylab.plot(X * scaleFactor, suppressor, **kwargs)

                pylab.axes(figCupric)
                pylab.plot(X * scaleFactor, cupric, **kwargs)
 
            # print '------------'
            # print self.parameter + ' value: ' + str(data[self.parameter])

        xticks = self._xticks()

        
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

        ax = pylab.axes(figDeposition)
        self._legend(ax)
        pylab.ylabel(r'$v$ $\left(\metre\per\second\right)$', rotation='vertical', fontsize=14)
        pylab.xlim(xmin=-maxFeatureDepth * scaleFactor)
        pylab.xlim(xmax=0)
        pylab.ylim(ymax=8e-9)
        pylab.ylim(ymin=0.0)
        pylab.xticks(xticks)
        pylab.title('(a)')
        if deposition_only:
            pylab.title('')
            self.legend =  pylab.legend(loc='upper right')

        if hasattr(self, 'legend'):
            self.legend.labelspacing = 0
            self.legend.columnspacing = 0.1
            for t in self.legend.texts:
                t.set_size(self.lfs)

        if not deposition_only:
            self.subplot(figDeposition)

        if type(filesuffix) is str:
            filesuffix = (filesuffix,)

        for fs in filesuffix:
            print self.parameter  + fs
            pylab.savefig(self.parameter + fs)
            pylab.close('all')

    def subplot(self, fig):
        pass

    def replaceString(self, Label):
        return Label

    def _maxSuppressor(self):
        return self.dataset[0]['bulkSuppressor']        

    def scale(self, variable):
        return variable

    def _xticks(self):
        return (-50, -40, -30, -20, -10, 0)

    def _colors(self):
        return [None] * len(self.dataset)
