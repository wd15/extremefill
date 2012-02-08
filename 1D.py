#!/usr/bin/env python

"""

Here we just integrate the 1D ODE that solves the electrical problem.
See http://matforge.org/wd15/blog/MoffatElectrical for details.

"""

import numpy

i0 = 40. ## A / m**2
alpha = 0.4
F = 9.6485e4 ## C / mol = J / V / mol
V_CELL = -0.275 ## V
R = 8.314 ## J / K / mol
T = 298.0 ## K
kappa = 15.26 ## S / m = A / V / m
delta = 100.e-6 ## m
C_DL = 0.3 ## F / m**2 = A s / V / m**2  
alphabar = alpha * F / R / T

def iF(V_IR):
    return i0 * numpy.exp(-alphabar * (V_CELL - V_IR))

def RHS(t, y):
    V_IR = y[0]
    return ((iF(V_IR) - kappa * V_IR / delta) / C_DL,)

def jacobian(t, y):
    V_IR = y[0]
    return ((alphabar * iF(V_IR) - kappa / delta) / C_DL,)


from scipy.integrate import ode

integrator = ode(RHS, jacobian).set_integrator('vode', method='bdf', with_jacobian=True)
integrator.set_initial_value(0., 0.)

times = []
V_IRs = []
dt = 1e-6
while integrator.successful() and integrator.t < 100.0e-6:
    integrator.integrate(integrator.t + dt)
    print integrator.t, integrator.y
    print 'accuracy',RHS(integrator.t, integrator.y)
    times.append(integrator.t)
    V_IRs.append(integrator.y[0])

import pylab
pylab.figure()
pylab.plot(times, V_IRs)
pylab.ylabel(r'$V_{IR}$')
pylab.xlabel(r'$t$')
pylab.savefig('voltageVersusTime.png')

pylab.figure()
pylab.plot(times, iF(numpy.array(V_IRs)))
pylab.ylabel(r'$i_F$')
pylab.xlabel(r'$t$')
pylab.savefig('currentVersusTime.png')

