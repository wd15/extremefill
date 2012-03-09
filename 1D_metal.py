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
V_APPLIED = -0.275 ## V
R = 8.314 ## J / K / mol
T = 298.0 ## K
kappa = 15.26 ## S / m = A / V / m
delta = 100e-8
delta_W = delta / 2 ## m
C_DL = 0.3 ## F / m**2 = A s / V / m**2  
D_c = 5e-10 ## m**2 / s
c_inf = 250. ## mol / m**3

Fbar = F / R / T
k = numpy.pi / 2 / delta

def iF0(V):
    return i0 * (numpy.exp(-alpha * Fbar * V) - numpy.exp((2 - alpha) * Fbar * V))

def iT(VW):
    return -kappa / delta_W * (V_APPLIED - VW)

def coeff(t):
    brac = numpy.exp(-k * D_c * t) / k - delta
    return brac / 2 / F / D_c / c_inf

def iFW(t, V):
    return iF0(V) / (1 - coeff(t) * iF0(V))

def cu(t, x, V):
    c_x = iFW(t, V) / 2 / F / D_c
    print 't',t
    print numpy.exp(-k * D_c * t)
    return c_x / k * numpy.cos(k * x) * numpy.exp(-k * D_c * t) + c_x * (x - delta) + c_inf

def iF0_V(V):
    return i0 * (-alpha * Fbar * numpy.exp(-alpha * Fbar * V) \
                      - (2 - alpha) * Fbar * numpy.exp((2 - alpha) * Fbar * V))

def iFW_V(t, V):
    return iF0_V(V) / (1 - coeff(t) * iF0(V))**2

def RHS(t, y):
    VW, VC = y
    return numpy.array([(iFW(t, VW) - iT(VW)) / C_DL,
                        (iF0(VC) + iT(VW)) / C_DL])

def jacobian(t, y):
    VW, VC = y
    return numpy.array([[(iFW_V(t, VW) - kappa / delta_W) / C_DL,
                         0],
                        [kappa / delta_W / C_DL,
                         iF0_V(VC) / C_DL]])

from scipy.integrate import ode

integrator = ode(RHS, jacobian).set_integrator('vode', method='bdf', with_jacobian=True)
integrator.set_initial_value((0., 0.), 0.)

times = []
V_DLs = []
dt = 1e-11
initialStep = True
    
while integrator.successful() and integrator.t < 1e6:
    integrator.integrate(integrator.t + dt)
    if initialStep is True:
        residuals0 = RHS(integrator.t, integrator.y)
        initialStep = False    
    print 'residual',RHS(integrator.t, integrator.y) / residuals0
    print 'integrator.t',integrator.t
    print 'V_DLs',integrator.y
    print
    times.append(integrator.t)
    V_DLs.append(integrator.y)
    dt *= 1.1

VW_DL = V_DLs[-1][0]
VC_DL = V_DLs[-1][1]
VW_IR = V_APPLIED - VW_DL
VC_IR = -(delta - delta_W) * VW_IR / delta_W
print 'VW_DL',VW_DL
print 'VC_DL',VC_DL 
print 'VW_IR',VW_IR
print 'VC_IR',VC_IR
print 'V_CELL',VW_DL + VW_IR - VC_IR - VC_DL
print 'i_T',iT(VW_DL)
print 'dimensionless parameter << 1:',delta * iF0(V_APPLIED) / c_inf / 2 / F / D_c

import pylab
pylab.figure()
pylab.semilogx(times, numpy.array(V_DLs)[:,0], times, numpy.array(V_DLs)[:,1])
pylab.ylabel(r'$V_{DL}$')
pylab.xlabel(r'$t$')
pylab.savefig('voltageVersusTimeCounter.png')

N = 1000
x = numpy.arange(N) / float(N - 1) * delta
fig1 = pylab.figure()

for i in range(len(V_DLs)):
    if i % 10 == 0:
        VW_DL = V_DLs[i][0]
        VC_DL = V_DLs[i][1]
        t = times[i]
        pylab.plot(x, cu(t, x, VW_DL))

pylab.ylabel(r'$c$')
pylab.xlabel(r'$x$')
pylab.savefig('copper.png')



# pylab.figure()
# pylab.plot(times, iF(numpy.array(V_IRs)))
# pylab.ylabel(r'$i_F$')
# pylab.xlabel(r'$t$')
# pylab.savefig('currentVersusTime.png')


