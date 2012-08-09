#!/usr/bin/env python

r"""

"""
__docformat__ = 'restructuredtext'

import fipy
from fipy import numerix

class Simulation(object):

    def __init__(self,  **kwargs):  
        self.baseParameters = {'dt' : .5e-7,
                               'dtMax' : 1e+20,
                               'dtMin' : .5e-7,
                               'totalSteps' : 400,
                               'view' : False,
                               'PRINT' : False,
                               'sweeps' : 5,
                               'tol' : 1e-10,
                               'delta' : 150e-6,
                               'deltaRef' : 0.03,
                               'featureDepth' : 56e-6,
                               'i1' : -40.,
                               'i0' : 40.,
                               'diffusionCupric' : 2.65e-10,
                               'appliedPotential' : -0.25,
                               'faradaysConstant' : 9.6485e4,
                               'gasConstant' : 8.314,
                               'temperature' : 298.,
                               'alpha' : 0.4,
                               'charge' : 2,
                               'bulkCupric' : 1000.,
                               'bulkSuppressor' : .02,
                               'diffusionSuppressor' : 9.2e-11,
                               'kappa' : 15.26,
                               'kPlus' : 150.,
                               'kMinus' : 2.45e7,
                               'omega' : 7.1e-6,
                               'gamma' : 2.5e-7,
                               'perimeterRatio' : 1. / 2.8e-6 * 0.093,
                               'areaRatio' : 0.093,
                               'capacitance' : 0.3}
        self.parameters = self.baseParameters.copy()
        for k, v in kwargs.items():
            self.parameters[k] = v
        
    def run(self):
        dt = self.parameters['dt']
        dtMax = self.parameters['dtMax']
        dtMin = self.parameters['dtMin']
        totalSteps = self.parameters['totalSteps']
        view = self.parameters['view']
        PRINT = self.parameters['PRINT']
        sweeps = self.parameters['sweeps']
        tol = self.parameters['tol']
        delta = self.parameters['delta']
        deltaRef = self.parameters['deltaRef']
        featureDepth = self.parameters['featureDepth']
        i1 = self.parameters['i1']
        i0 = self.parameters['i0']
        diffusionCupric = self.parameters['diffusionCupric']
        appliedPotential = self.parameters['appliedPotential']
        faradaysConstant = self.parameters['faradaysConstant']
        gasConstant = self.parameters['gasConstant']
        temperature = self.parameters['temperature']
        alpha = self.parameters['alpha']
        charge = self.parameters['charge']
        bulkCupric = self.parameters['bulkCupric']
        bulkSuppressor = self.parameters['bulkSuppressor']
        diffusionSuppressor = self.parameters['diffusionSuppressor']
        kappa = self.parameters['kappa']
        kPlus = self.parameters['kPlus']
        kMinus = self.parameters['kMinus']
        omega = self.parameters['omega']
        gamma = self.parameters['gamma']
        perimeterRatio = self.parameters['perimeterRatio']
        areaRatio = self.parameters['areaRatio']
        capacitance = self.parameters['capacitance']

        Fbar = faradaysConstant / gasConstant / temperature
        self.parameters['Fbar'] = Fbar
        epsilon = 1e-30 

        L = delta + featureDepth
        N = 1000
        dx = L / N 
        mesh = fipy.Grid1D(nx=N, dx=dx) - [[featureDepth]]

        potential = fipy.CellVariable(mesh=mesh, hasOld=True, name=r'$\psi$')
        potential[:] = -appliedPotential

        cupric = fipy.CellVariable(mesh=mesh, hasOld=True, name=r'$c_{cu}$')
        cupric[:] = bulkCupric
        cupric.constrain(bulkCupric, mesh.facesRight)

        suppressor = fipy.CellVariable(mesh=mesh, hasOld=True, name=r'$c_{\theta}$')
        suppressor[:] = bulkSuppressor
        suppressor.constrain(bulkSuppressor, mesh.facesRight)

        theta = fipy.CellVariable(mesh=mesh, hasOld=True, name=r'$\theta$')

        I0 = (i0 + i1 * theta)
        baseCurrent = I0 * (numerix.exp(alpha * Fbar * potential) \
                                - numerix.exp(-(2 - alpha) * Fbar * potential))
        cbar =  cupric / bulkCupric
        current = cbar * baseCurrent
        currentDerivative = cbar * I0 * (alpha * Fbar *  numerix.exp(alpha * Fbar * potential) \
                                             + (2 - alpha) * Fbar * numerix.exp(-(2 - alpha) * Fbar * potential))

        def dirac(x):
            value = numerix.zeros(mesh.numberOfCells, 'd')
            ID = numerix.argmin(abs(mesh.x - x))
            if mesh.x[ID] < x:
                ID = ID + 1
            value[ID] = 1. / dx
            return value

        THETA = (mesh.x < 0) * perimeterRatio + dirac(0) * (1 - areaRatio) + dirac(-featureDepth) * areaRatio
        AREA = (mesh.x < 0) * (areaRatio - 1) + 1 
        THETA_UPPER = fipy.CellVariable(mesh=mesh)
        THETA_UPPER[-1] = kappa / dx / (deltaRef - delta)

        potentialEq = fipy.TransientTerm(capacitance * THETA) == fipy.DiffusionTerm(kappa * AREA) \
            - THETA * (current - potential * currentDerivative) \
            - fipy.ImplicitSourceTerm(THETA * currentDerivative) \
            - THETA_UPPER * appliedPotential - fipy.ImplicitSourceTerm(THETA_UPPER) 

        cupricEq = fipy.TransientTerm(AREA) == fipy.DiffusionTerm(diffusionCupric * AREA) \
            - fipy.ImplicitSourceTerm(baseCurrent * THETA / (bulkCupric * charge * faradaysConstant))

        suppressorEq = fipy.TransientTerm(AREA) == fipy.DiffusionTerm(diffusionSuppressor * AREA) \
            - fipy.ImplicitSourceTerm(gamma * kPlus * (1 - theta) * THETA)

        thetaEq =  fipy.TransientTerm(THETA + epsilon) == kPlus * suppressor * THETA \
            - fipy.ImplicitSourceTerm(THETA * (kPlus * suppressor + kMinus * current * (omega / charge / faradaysConstant)))

        t = 0.

        if view:
            potentialBar = -potential / appliedPotential
            potentialBar.name = r'$\bar{\eta}$'
            cbar.name = r'$\bar{c_{cu}}$'
            suppressorBar = suppressor / bulkSuppressor
            suppressorBar.name = r'$\bar{c_{\theta}}$'

            viewer = fipy.Viewer((theta, suppressorBar, cbar, potentialBar), datamax=1, datamin=0.0)

        for step in range(totalSteps):
            if view:
                viewer.axes.set_title(r'$t=%1.2e$' % t)
                viewer.plot()

            potential.updateOld()
            cupric.updateOld()
            suppressor.updateOld()
            theta.updateOld()

            for sweep in range(sweeps):
                potentialRes = potentialEq.sweep(potential, dt=dt)
                cupricRes = cupricEq.sweep(cupric, dt=dt)
                suppressorRes = suppressorEq.sweep(suppressor, dt=dt)
                thetaRes = thetaEq.sweep(theta, dt=dt)
                import numpy
                res = numpy.array((potentialRes, cupricRes, suppressorRes, thetaRes))
                if sweep == 0:
                    res0 = res
                else:
                    if ((res / res0) < tol).all():
                        break

                if PRINT:
                    print res / res0

            if sweep == sweeps - 1 and PRINT:
                print 'Did not reach sufficient tolerance'
                print 'kPlus',kPlus
                print 'kMinus',kMinus
                print 'res',res

            if PRINT:
                print 'theta',theta[0]
                print 'cBar_supp',suppressor[0] / bulkSuppressor
                print 'cBar_cu',cupric[0] / bulkCupric
                print 'potentialBar',-potential[0] / appliedPotential
                print 'dt',dt
                print 'step',step

            t += dt

            dt = dt * 1e+10
            dt = min((dt, dtMax))
            dt = max((dt, dtMin))

        if view:
            viewer.plot()

        self.parameters['potential'] = numpy.array(potential)
        self.parameters['cupric'] = numpy.array(cupric)
        self.parameters['suppressor'] = numpy.array(suppressor)
        self.parameters['theta'] = numpy.array(theta)

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())


