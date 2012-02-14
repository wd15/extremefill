#!/usr/bin/env python

"""

Here we just integrate the 1D ODE that solves the electrical problem.
See http://matforge.org/wd15/blog/MoffatElectrical for details.

"""

import numpy

## numbers
i0 = 40. ## A / m**2
alpha = 0.4
F = 9.6485e4 ## C / mol = J / V / mol
V_CELL = -0.275 ## V
R = 8.314 ## J / K / mol
T = 298.0 ## K
kappa = 15.26 ## S / m = A / V / m
delta = 100.e-6 ## m
C_DL = 0.3 ## F / m**2 = A s / V / m**2  
Fbar = F / R / T


def iT(VW, VC):
    return -kappa / delta * (V_CELL - VW + VC)

def iF(V):
    return i0 * (numpy.exp(-alpha * Fbar * V) - numpy.exp((2 - alpha) * Fbar * V))

def iF_V(V):
    return i0 * (-alpha * Fbar * numpy.exp(-alpha * Fbar * V) \
                      - (2 - alpha) * Fbar *  numpy.exp((2 - alpha) * Fbar * V))

def RHS(t, y):
    VW, VC = y
    return numpy.array([(iF(VW) - iT(VW, VC)) / C_DL,
                        (iF(VC) + iT(VW, VC)) / C_DL])

def jacobian(t, y):
    VW, VC = y
    return numpy.array([[(iF_V(VW) - kappa / delta) / C_DL,
                         kappa / delta / C_DL],
                        [kappa / delta / C_DL,
                         (iF_V(VC) - kappa / delta) / C_DL]])

from scipy.integrate import ode

integrator = ode(RHS, jacobian).set_integrator('vode', method='bdf', with_jacobian=True)
integrator.set_initial_value((V_CELL, 0.), 0.)

times = []
V_DLs = []
dt = 1e-6
initialStep = True
while integrator.successful() and integrator.t < 100.0e-6:
    integrator.integrate(integrator.t + dt)
    if initialStep is True:
        residuals0 = RHS(integrator.t, integrator.y)
        initialStep = False    
    print 'residual',RHS(integrator.t, integrator.y) / residuals0
    times.append(integrator.t)
    V_DLs.append(integrator.y)

print 'V_DL',V_DLs[-1]
print 'iT',iT(V_DLs[-1][0], V_DLs[-1][1])

import pylab
pylab.figure()
pylab.plot(times, numpy.array(V_DLs)[:,0], times, numpy.array(V_DLs)[:,1])
pylab.ylabel(r'$V_{DL}$')
pylab.xlabel(r'$t$')
pylab.savefig('voltageVersusTimeCounter.png')

# pylab.figure()
# pylab.plot(times, iF(numpy.array(V_IRs)))
# pylab.ylabel(r'$i_F$')
# pylab.xlabel(r'$t$')
# pylab.savefig('currentVersusTime.png')


