#!/usr/bin/env python

r"""

This example is a simple test for solving the potential equation in
1D. We solve the potential equation,

.. math::

    \nabla^2  \psi = 0

where :math:`\psi` is the electric potential. The boundary
condition at :math:`\phi=0` (the working electrode) is given by,

.. math::

    c_{DL} \frac{\partial \psi}{\partial t} + i_{F} \left( \psi - \psi_0 \right)  = 
    -\kappa \vec{n} \cdot \nabla \psi

where :math:`c_{DL}` is the capacitance, :math:`i_F` is the
current density, :math:`\kappa` is the conductance and the normal
points out of the electrolyte. The quantity :math:`\psi_0` is the value
of :math:`\psi` at the reference electrode. 


The boundary condition at the counter electode is :math:`\psi=0`. The
current is given by

.. math::

    i_F \left( \xi \right) = i_0 \left[\exp{\left(-\frac{\alpha F \left(V_{\text{APPLIED}} - \xi\right) }{R T} \right)} -  \exp{\left(\frac{\left(2 -\alpha\right) F \left(V_{\text{APPLIED}} - \xi \right)}{R T} \right)}  \right]

The paramaters are as follows.

>>> import numpy
>>> delta = 100e-6 ## m
>>> deltaRef = 50e-6 ## m
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

    c_{DL} \frac{\partial \psi}{\partial t} + i_{F} \left( \psi - \psi_0 \right)  = 
    -\frac{\kappa \psi}{\delta}

where :math:`\psi` is now only defined at the interface. Using scipy we can solve the ODE.

>>> def iF(xi):
...     V = appliedVoltage - xi
...     return i0 * (numpy.exp(-alpha * Fbar * V) - numpy.exp((2 - alpha) * Fbar * V))

>>> def iFderivative(xi):
...     V = appliedVoltage - xi
...     return i0 * (alpha * Fbar * numpy.exp(-alpha * Fbar * V) \
...            + (2 - alpha) * Fbar * numpy.exp((2 - alpha) * Fbar * V))

>>> def RHS(t, y):
...     psi = y[0]
...     psi0 = (1 - deltaRef / delta) * psi
...     return numpy.array((-iF(psi - psi0) / capacitance  - kappa * psi / delta / capacitance,))

>>> def jacobian(t, y):
...     psi = y[0]
...     psi0 = (1 - deltaRef / delta) * psi
...     return numpy.array((-iFderivative(psi - psi0) / capacitance  - kappa / delta / capacitance,))

>>> from scipy.integrate import ode
>>> integrator = ode(RHS, jacobian).set_integrator('vode', method='bdf', with_jacobian=True)
>>> initialValue = appliedVoltage + (delta - deltaRef) / deltaRef * appliedVoltage
>>> s = integrator.set_initial_value((initialValue,), 0.)

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
>>> potentialVar[:] = initialValue
>>> potentialVar.constrain(0, mesh.facesRight)

>>> potentialRef = CellVariable(mesh=mesh, hasOld=True)

>>> xi = potentialVar - potentialRef

>>> potentialEqn = PotentialEquation(var=potentialVar,
...                                  distanceVar=distanceVar,
...                                  currentDensity=iF(xi),
...                                  conductivity=kappa,
...                                  capacitance=capacitance)

>>> potentialFiPy = numpy.zeros(totalSteps)
>>> potentialRefs = numpy.zeros(totalSteps)

>>> t = 0.
>>> for step in range(totalSteps):
...     potentialVar.updateOld()
...     for sweep in range(5):
...         potentialRef[:] = potentialVar([[deltaRef]])
...         residual = potentialEqn.sweep(potentialVar, dt=dt)
...     potentialFiPy[step] = potentialVar.value[interfaceID]
...     potentialRefs[step] = potentialVar([[deltaRef]])

>>> print numpy.allclose(potentialFiPy, potentialSciPy, atol=1e-3)
True

>>> if __name__ == '__main__':
...     import pylab
...     f = pylab.figure()
...     pylab.plot(times, potentialSciPy, times, potentialFiPy)
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
