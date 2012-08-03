#!/usr/bin/env python

r"""

This example is the 1D solution, but with a via or trench, where we
assume that there is no lateral variation in any of the fields and the
deposition rate is negligible. Let's start with the equations in order
of potential, cupric concentration and suppressor concentration and
surfactant supressor integrated over space.

.. math::

    \int_{-h}^{\delta}  dz \left[ c_{DL} \frac{\partial \psi}{\partial t} \Theta \left( z \right) =
    \kappa \frac{\partial^2 \psi}{\partial z^2} A \left(z \right) -
    i_F \left(z\right) \Theta \left(z \right)  \right]

..
.. math::

    \int_{-h}^{\delta}  dz  \left[\frac{\partial c_{\text{cu}}}{\partial t} A \left( z \right)  =
    D_{cu} \frac{\partial^2 c_{\text{cu}}}{\partial z^2} A \left(z \right)  -
    \frac{i_F\left(z\right)}{n F} \Theta \left(z \right) \right]

..
.. math::

    \int_{-h}^{\delta}  dz \left[\frac{\partial c_{\theta}}{\partial t} A \left( z \right) =
    D_{\theta} \frac{\partial^2 c_{\theta}}{\partial z^2} A \left(z \right)  -
    \Gamma k^+ c_{\theta} \left(1 - \theta \right)  \Theta \left(z \right)  \right]

..
.. math::

    \int_{-h}^{\delta}  dz \left[\frac{\partial \theta}{\partial t} \Theta \left( z \right) =
    k^+ c_{\theta} \left(1 - \theta \right)  \Theta \left(z \right)  -
    k^- \theta \frac{i_F \Omega}{n F} \Theta \left(z \right) \right]

where :math:`\Theta\left( z \right) = H\left(- z \right) \frac{l
\left( z \right)}{A_F} + \delta \left( z \right)\left(1 -
\frac{A_T}{A_F} \right) + \delta(z + h) \frac{A_T}{A_F}`. The
cross-sectional area ratio is given by,

.. math::

   A \left( z \right) = \begin{cases}
   1 & \text{when $z > 0$,} \\
   \frac{A_T}{A_F} & \text{when $z \le 0$,}
   \end{cases}

where :math:`A_F` is the cross-sectional area above of the modeling
domain and :math:`A_S` is the cross-sectional area in the
trench/via. The length of the perimeter is given by
:math:`l\left(z\right)` and is a step-function through 0. Also,
:math:`\delta\left(z\right)` is the Dirac delta function and
:math:`H\left(z\right)` is the Heaviside step function with
:math:`z=-h` at the bottom of the trench.  The current density is
given by,

.. math::

    i_F \left( \xi, \theta \right) = \left(i_0 + i_{\theta} \theta\right) \left[\exp{\left(-\frac{\alpha F \left(\eta - \xi\right) }{R T} \right)} -  \exp{\left(\frac{\left(2 -\alpha\right) F \left(\eta - \xi \right)}{R T} \right)}  \right]

where :math:`\xi = \psi - \psi_{\text{ref}}`. The boundary conditions on the working electrode are
included in the volume integrals. Additionally,

.. math::
 
     \psi = 0 \;\; & \text{at $z = \delta$} \\
     \psi = \frac{\delta}{\delta_{\text{ref}}} \eta \;\; & \text{at $t = 0$ forcing $\xi=\eta$ initially} \\
     c_{\text{cu}} = c_{\text{cu}}^{\infty} \;\; & \text{at $z = \delta$} \\
     c_{\text{cu}} = c_{\text{cu}}^{\infty} \;\; & \text{at $t = 0$} \\
     c_{\theta} = c_{\theta}^{\infty} \;\; & \text{at $z = \delta$} \\
     c_{\theta} = c_{\theta}^{\infty} \;\; & \text{at $t = 0$}\\
     \theta=0 \;\; & \text{at $t = 0$}

and

.. math::

     \psi_{\text{ref}} = \psi\left(0\right) \left(1 - \frac{\delta_{\text{ref}}}{\delta}\right)

The following code compares the full 1D feature model (but with no
feature) with the simple 1D ODE for solving the electrical equation
with no suppressor and no cupric depeletion.

>>> import pylab
>>> times, potentials = feature1D(delta=100e-6,
...                               deltaRef=50e-6,
...                               featureDepth=0.0,
...                               i0=40.,
...                               i1=0.0,
...                               diffusionCupric=1e+10)

>>> timesScipy, potentialsScipy = noFeatureODE()

>>> print max(abs(potentialsScipy - potentials))
>>> numerix.allclose(potentials, potentialsScipy)
True

>>> pylab.figure()
>>> pylab.plot(times, potentials, timesScipy, potentialsScipy)
>>> pylab.ylabel(r'$\psi\left(0\right)$ (V)')
>>> pylab.xlabel(r'$t$ (s)')
>>> pylab.savefig('FiPyVScipy.png')

Agreement is good for :math:`\psi`.

.. image:: FiPyVScipy.*
   :width: 90%
   :align: center
   :alt: comparison of FiPy and SciPy for the potential equation

"""
__docformat__ = 'restructuredtext'

from fipy import Grid1D, CellVariable, Variable, numerix, TransientTerm, DiffusionTerm, ImplicitSourceTerm, Viewer

def feature1D(delta=150e-6,
              deltaRef=0.03,
              featureDepth=56e-6,
              i1=-40.,
              i0=40.,
              diffusionCupric=5.6e-10,
              dt=.5e-7,
              dtMax=.5e-7,
              totalSteps=200,
              eta=-0.275,
              view=False):

    F = 9.6485e4 ## C / mol = J / V / mol
    R = 8.314 ## J / K / mol
    T = 298. ## K
    Fbar = F / R / T ## 1 / V
    alpha = 0.4
    capicatance = 0.3 ## F / m**2 = A s / V / m**2  
    kappa = 15.26 ## S / m = A / V / m
    kPlus = 125. ## m**3 / mol / s
    kMinus = 2.45e7 ## 1 / m
    charge = 2
    areaRatio = 0.093
    perimeterRatio = 1. / 2.8e-6  
    cupric_inf = 1000.
    suppressor_inf = .02
    diffusion_suppressor = 1e-9
    gamma = 2.5e-7
    epsilon = 1e-30
    omega = 7.1e-6

    L = delta + featureDepth
    N = 400
    dx = L / N 
    mesh = Grid1D(nx=N, dx=dx) - [[featureDepth]]

    potential = CellVariable(mesh=mesh, hasOld=True, name=r'$\psi$')
    potential[:] = 0 ##delta * eta / deltaRef 
    potential.constrain(0, mesh.facesRight)
    potentialRef = Variable()

    cupric = CellVariable(mesh=mesh, hasOld=True, name=r'$c_{cu}$')
    cupric[:] = cupric_inf
    cupric.constrain(cupric_inf, mesh.facesRight)

    suppressor = CellVariable(mesh=mesh, hasOld=True, name=r'$c_{\theta}$')
    suppressor[:] = suppressor_inf
    suppressor.constrain(suppressor_inf, mesh.facesRight)

    theta = CellVariable(mesh=mesh, hasOld=True, name=r'$\theta$')

    V = eta - potential + potentialRef
    baseCurrent = (i0 + i1 * theta) * (numerix.exp(-alpha * Fbar * V) - numerix.exp((2 - alpha) * Fbar * V))
    current = cupric / cupric_inf * baseCurrent

    def dirac(x):
        value = numerix.zeros(mesh.numberOfCells, 'd')
        ID = numerix.argmin(abs(mesh.x - x))
        if mesh.x[ID] < x:
            ID = ID + 1
        value[ID] = 1. / dx
        return value

    THETA = (mesh.x < 0) * perimeterRatio + dirac(0) * (1 - areaRatio) + dirac(-featureDepth) * areaRatio 
    AREA = (mesh.x < 0) * (areaRatio - 1) + 1 

    potentialEq = TransientTerm(capicatance * THETA) == DiffusionTerm(kappa * AREA) - \
        current * THETA

    cupricEq = TransientTerm(AREA) == DiffusionTerm(diffusionCupric * AREA) \
        - ImplicitSourceTerm(baseCurrent * THETA / cupric_inf)

    suppressorEq = TransientTerm(AREA) == DiffusionTerm(diffusion_suppressor * AREA) \
        - ImplicitSourceTerm(gamma * kPlus * (1 - theta) * THETA)

    thetaEq =  TransientTerm(THETA + epsilon) == kPlus * suppressor * THETA \
        - ImplicitSourceTerm(THETA * (kPlus * suppressor + kMinus * current * (omega / charge / F)))

    t = 0.

    if view:
        viewers = (Viewer(V, datamax=0, datamin=eta), Viewer(cupric), Viewer(suppressor), Viewer(theta))

    potentials = []
    times = []

    for step in range(totalSteps):
        potential.updateOld()
        cupric.updateOld()
        suppressor.updateOld()
        theta.updateOld()
##        print 'step',step
        for sweep in range(5):
            potentialRef.setValue(potential([[0]]) * (1 - deltaRef / delta))
            potentialRes = potentialEq.sweep(potential, dt=dt)
            cupricRes = cupricEq.sweep(cupric, dt=dt)
            suppressorRes = suppressorEq.sweep(suppressor, dt=dt)
            thetaRes = thetaEq.sweep(theta, dt=dt)
            
        if view:
            for viewer in viewers:
                viewer.plot()
                  
##        print 'potential([[0]]) - potentialRef',potential([[0]]) - potentialRef
##        print 'cupric([[0]]) / cupric_inf',cupric([[0]]) / cupric_inf
        t += dt
        times += [t]
        print potential([[0]])
        raw_input('stopped')
        potentials += [potential([[0]])]
        dt = dt * 1.1
        dt = min((dt, dtMax))
##        print 'time',t

    if view:
        for viewer in viewers:
            viewer.plot()
        
    return times, potentials 

def noFeatureODE():
    import numpy
    delta = 100e-6 ## m
    deltaRef = 50e-6 ## m
    faradaysConstant = 9.6485e4 ## C / mol = J / V / mol
    gasConstant = 8.314 ## J / K / mol
    temperature = 298. ## K
    Fbar = faradaysConstant / gasConstant / temperature ## 1 / V
    alpha = 0.4
    appliedVoltage = -0.275  ## V
    i0 = 40. ## A / m**2 
    capacitance = 0.3 ## F / m**2 = A s / V / m**2  
    kappa = 15.26 ## S / m = A / V / m

    def iF(xi):
        V = appliedVoltage - xi
        return i0 * (numpy.exp(-alpha * Fbar * V) - numpy.exp((2 - alpha) * Fbar * V))

    def iFderivative(xi):
        V = appliedVoltage - xi
        return i0 * (alpha * Fbar * numpy.exp(-alpha * Fbar * V) \
               + (2 - alpha) * Fbar * numpy.exp((2 - alpha) * Fbar * V))

    def RHS(t, y):
        psi = y[0]
        psi0 = (1 - deltaRef / delta) * psi
        return numpy.array((-iF(psi - psi0) / capacitance  - kappa * psi / delta / capacitance,))

    def jacobian(t, y):
        psi = y[0]
        psi0 = (1 - deltaRef / delta) * psi
        return numpy.array((-iFderivative(psi - psi0) / capacitance  - kappa / delta / capacitance,))

    from scipy.integrate import ode
    integrator = ode(RHS, jacobian).set_integrator('vode', method='bdf', with_jacobian=True)
    initialValue = 0 ##appliedVoltage + (delta - deltaRef) / deltaRef * appliedVoltage
    s = integrator.set_initial_value((initialValue,), 0.)

    totalSteps = 200
    dt = .5e-7
    times = numpy.zeros(totalSteps)
    potentialSciPy = numpy.zeros(totalSteps)
    step = 0

    while integrator.successful() and step < totalSteps:
        null = integrator.integrate(integrator.t + dt)
        times[step] = integrator.t
        potentialSciPy[step] = integrator.y[0]
        step += 1

    return times, potentialSciPy

if __name__ == '__main__':
##    times, potentials = feature1D(featureDepth=50e-6, dt=1e-7, dtMax=1., totalSteps=1000, eta=-0.275)
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())


