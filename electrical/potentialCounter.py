#!/usr/bin/env python

r"""


Including the effects of the counter electrode.

At the working electrode :math:`W`:

.. math::

    c_{DL} \frac{\partial \psi}{\partial t} + i_{F} \left( \psi \right)  = 
    -\kappa \vec{n} \cdot \nabla \psi \\
    i_F = i_0 \left[\exp{\left(-\frac{\alpha F \left(V_{\text{APPLIED}} - \psi\right) }{R T} \right)} -  \exp{\left(\frac{\left(2 -\alpha\right) F \left(V_{\text{APPLIED}} - \psi\right)}{R T} \right)}  \right]

At the counter electrode :math:`C`:

.. math::

    i_T^C =& -\kappa \vec{n}^C \cdot \nabla V_{IR} \\
    i_T^C =& -c_{DL} \frac{\partial V_{DL}^C}{\partial t} + i_{F}^C \\
    i_F^C =& i_0^C \left[ \exp{\left(-\frac{\alpha F V_{DL}^C}{R T} \right)} -  \exp{\left(\frac{\left(2 -\alpha\right) F V_{DL}^C}{R T} \right)}  \right] \\


The normals point out of the electrolyte.

In the bulk:

$ \nabla^2 V_{IR} = 0 $

and

$V_{CELL} = V_{DL}^W + V_{IR}^W - V_{DL}^C - V_{IR}^C$

$V_{APPLIED} = V_{DL}^W + V_{IR}^W$

$\delta_W V_{IR}^C + \delta_C V_{IR}^W = 0 $

$\delta_C + \delta_W = \delta$

where $\delta_C$ and $\delta_W$ are the distances from the counter and
working electrodes to the reference electode, respectively. $\delta$
is the total distance between the electrodes.

If we do this the equations can be reduced to two coupled ODEs for
$V_{DL}^C$ and $V_{DL}^W$:

$ -\frac{\kappa}{\delta_W} \left(V_{APPLIED} - V_{DL}^W\right) = -c_{DL} \frac{\partial V_{DL}^W}{\partial t} +  i_0^W \left[ \exp{\left(-\frac{\alpha F V_{DL}^W}{R T} \right)} -  \exp{\left(\frac{\left(2 -\alpha\right) F V_{DL}^W}{R T} \right)}  \right]$

and 

$ \frac{\kappa}{\delta_W} \left(V_{APPLIED} - V_{DL}^W\right) = -c_{DL} \frac{\partial V_{DL}^C}{\partial t} +  i_0^C \left[ \exp{\left(-\frac{\alpha F V_{DL}^C}{R T} \right)} -  \exp{\left(\frac{\left(2 -\alpha\right) F V_{DL}^C}{R T} \right)}  \right]$

The material parameters are:

$i_0^W = i_0^C = 40\;\text{A / m}^2$

$\alpha = 0.4$

$F = 9.6485 \times 10^{-4}\;\text{J / V / mol}$

$V_{APPLIED} = -0.275\;\text{V}$

$R = 8.314\;\text{J / K / mol}$

$T = 298.0\;\text{K}$

$\kappa = 15.26\;\text{A / V / m}$

$\delta = 100 \times10^{-6}\;\text{m}$

$\delta_W = 50 \times10^{-6}\;\text{m}$

$C_{DL} = 0.3\;\text{A s / V / m}^2$

Integrating this gives:

$V_{DL}^W = -0.2666$

$V_{DL}^C = 0.0667$

$i_T^C = 2546$

$V_{IR}^W =  -0.008345$

$V_{IR}^C =  0.008345$

$V_{CELL} = -0.350$

The code for this is source:trunk/moffat/electrical/misc/1D_counter.py.

The evolution of the $V_{DL}$:

[[Image(source:trunk/moffat/electrical/misc/voltageVersusTimeCounter.png, 600px)]]



This example is a simple test for solving the potential equation in
1D with the counter electrode. We solve the potential equation,

.. math::

    \nabla^2  \psi = 0

where :math:`\psi` is the electric potential. The boundary
condition at :math:`\phi=0` (the working electrode) is given by,

.. math::

    c_{DL} \frac{\partial \psi}{\partial t} + i_{F} \left( \psi \right)  = 
    -\kappa \vec{n} \cdot \nabla \psi

where :math:`c_{DL}` is the capacitance, :math:`i_F` is the
current density, :math:`\kappa` is the conductance and the normal
points out of the electrolyte. The boundary condition at the counter electode is :math:`\psi=0`. The current is given by

.. math::

    i_F = i_0 \left[\exp{\left(-\frac{\alpha F \left(V_{\text{APPLIED}} - \psi\right) }{R T} \right)} -  \exp{\left(\frac{\left(2 -\alpha\right) F \left(V_{\text{APPLIED}} - \psi\right)}{R T} \right)}  \right]

The paramaters are as follows.

>>> import numpy
>>> delta = 50e-6 ## m
>>> faradaysConstant = 9.6485e4 ## C / mol = J / V / mol
>>> gasConstant = 8.314 ## J / K / mol
>>> temperature = 298. ## K
>>> Fbar = faradaysConstant / gasConstant / temperature ## 1 / V
>>> alpha = 0.4
>>> appliedVoltage = -0.275  ## V
>>> i0 = 40. ## A / m**2 
>>> capacitance = 0.3 ## F / m**2 = A s / V / m**2  
>>> kappa = 15.26 ## S / m = A / V / m

In 1D, we can write an ODE for :math:`\psi` at the interface.

.. math::

    c_{DL} \frac{\partial \psi}{\partial t} + i_{F} \left( \psi \right)  = 
    -\frac{\kappa \psi}{\delta}

where :math:`\psi` is now only defined at the interface. Using scipy we can solve the ODE.

>>> def iF(psi):
...     V = appliedVoltage - psi
...     return i0 * (numpy.exp(-alpha * Fbar * V) - numpy.exp((2 - alpha) * Fbar * V))

>>> def iFderivative(psi):
...     V = appliedVoltage - psi
...     return i0 * (alpha * Fbar * numpy.exp(-alpha * Fbar * V) \
...            + (2 - alpha) * Fbar * numpy.exp((2 - alpha) * Fbar * V))

>>> def RHS(t, y):
...     psi = y[0]
...     return numpy.array((-iF(psi) / capacitance  - kappa * psi / delta / capacitance,))

>>> def jacobian(t, y):
...     psi = y[0]
...     return numpy.array((-iFderivative(psi) / capacitance  - kappa / delta / capacitance,))

>>> from scipy.integrate import ode
>>> integrator = ode(RHS, jacobian).set_integrator('vode', method='bdf', with_jacobian=True)
>>> s = integrator.set_initial_value((appliedVoltage,), 0.)

>>> totalSteps = 1000
>>> dt = 1e-8
>>> times = numpy.zeros(totalSteps)
>>> potentialSciPy = numpy.zeros(totalSteps)
>>> step = 0

We iterate for ``totalSteps`` and save the results for comparison with
the FiPy solution.

>>> while integrator.successful() and step < totalSteps:
...     null = integrator.integrate(integrator.t + dt)
...     times[step] = integrator.t
...     potentialSciPy[step] = integrator.y[0]
...     step += 1

We solve this problem with FiPy using the ``PotentialEquation``.

>>> from fipy import Grid1D, DistanceVariable, CellVariable
>>> from potentialEquation import PotentialEquation

>>> L = delta * 1.1
>>> N = 110
>>> interfaceID = 10

>>> mesh = Grid1D(Lx=L, nx=N) - [[L - delta]]

>>> distanceVar = DistanceVariable(mesh=mesh)
>>> distanceVar[:] = -1.
>>> distanceVar.setValue(1., where=mesh.x > 0)
>>> distanceVar.calcDistanceFunction()

>>> potentialVar = CellVariable(mesh=mesh, hasOld=True)
>>> potentialVar[:] = appliedVoltage
>>> potentialVar.constrain(0, mesh.facesRight)

>>> potentialEqn = PotentialEquation(var=potentialVar,
...                                  distanceVar=distanceVar,
...                                  currentDensity=iF(potentialVar),
...                                  conductivity=kappa,
...                                  capacitance=capacitance)

>>> potentialFiPy = numpy.zeros(totalSteps)

>>> t = 0.
>>> for step in range(totalSteps):
...     potentialVar.updateOld()
...     for sweep in range(5):
...         residual = potentialEqn.sweep(potentialVar, dt=dt)
...     potentialFiPy[step] = potentialVar.value[interfaceID]

>>> print numpy.allclose(potentialFiPy, potentialSciPy, atol=1e-4)
True

>>> if __name__ == '__main__':
...     import pylab
...     f = pylab.figure()
...     pylab.plot(times, potentialSciPy, times, potentialFiPy)
...     pylab.savefig('potential.png')
...     pylab.show()


The results are in close agreement.

.. image:: potential.*
   :width: 90%
   :align: center
   :alt: comparison of FiPy and SciPy for the potential equation



"""

__docformat__ = 'restructuredtext'

__all__ = []
if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
