#!/usr/bin/env python

r"""

This example solves the full steady state system. The steady state solution
is given by,

.. math::

    \theta =& \frac{1}{\bar{B}} \frac{\bar{c}_{\theta}}{E \bar{c}_{\text{cu}}}\;\; \text{or} \;\; 1
    \\
    \bar{c}_{\theta} =& \frac{1}{1 + \bar{K} \left(1 - \theta \right)}
    \\
    \bar{c}_{\text{cu}} =& \frac{1}{1 + \bar{C} E \left(1 - \theta \right)}
    \\
    \bar{\xi} =& -\bar{G} \left(1 - \theta \right) 
    E \bar{c}_{\text{cu}}
    
where

.. math::

    \bar{B} = \frac{k^-}{k^+} \frac{\Omega}{n F} \frac{i_0}{c_{\theta}^{\infty}} 
    \text{,}\;\;
    \bar{K} = \frac{\Gamma k^+ \delta}{D_{\theta}}
    \text{,}\;\;
    \bar{C} = \frac{\delta i_0}{D_{\text{cu}} n F c_{\text{cu}}^{\infty}}
    \text{,}\;\;
    \bar{G} = \frac{\delta_{\text{ref}} i_0}{\kappa \eta}
    \text{,}\;\;
    \bar{F} = \frac{F \eta }{R T}

and 

.. math::

    E\left(\bar{\xi}\right) = \left[\exp{\left(-\alpha \bar{F} \left(1 - \bar{\xi}\right) \right)} - 
    \exp{\left(2 - \alpha\right) \bar{F} \left(1 - \bar{\xi}\right)}  \right]

..
.. math::

    \bar{c}_{\text{cu}} = \frac{c_{\text{cu}}}{c_{\text{cu}}^{\infty}}
    \text{,}\;\;
    \bar{c}_{\theta} = \frac{c_{\theta}}{c_{\theta}^{\infty}}
    \text{,}\;\;
    \bar{\xi} = \frac{\xi}{\eta} 

The whole system is controlled by only six parameters including
:math:`\alpha`. All over bars indicate dimensionless quantities whether
variables or parameters.

If :math:`B = \frac{\bar{B} E
\bar{c}_{\text{cu}}}{\bar{c}_{\theta}}`, for :math:`B<1`, there is
only one solution at :math:`\theta=1`, this is an attractor (for all
:math:`0 \le \theta \le 1`). For :math:`B>1`, the solution is
:math:`\frac{1}{B}`. This is also a stable solution as
:math:`\dot{\theta}` is positive for :math:`\theta<\frac{1}{B}` and
negative for :math:`\theta>\frac{1}{B}`.

When we solve the trench we need :math:`B < 1` on the top surface and
:math:`B > 1`, at the bottom of the trench. That is essentially the
whole problem.

For :math:`\theta=1`, the whole system shuts down
:math:`\bar{c}_{\theta}=1`, :math:`\bar{c}_{\text{cu}}=1` and
:math:`\bar{\psi}=0`. If :math:`B>1` this is not a stable solution
though so the system should migrate to some value of :math:`\theta<1`
with any perturbation.

I've set up a function ``solve`` that returns the :math:`\theta`,
:math:`\bar{c}_{\theta}`, :math:`\bar{c}_{\text{cu}}` and
:math:`\bar{\psi}` in that order.  As a very trivial test case, lets
find the roots for the most trivial system

>>> B = 10.
>>> F = 9.6485e4
>>> R = 8.314
>>> T = 298.
>>> eta = -0.275
>>> alpha = 0.4
>>> Fbar = F * eta / R / T
>>> Bbar = B / E(0, Fbar, alpha)
>>> print solve(Bbar=Bbar, Kbar=0, Cbar=0, Gbar=0, Fbar=Fbar, alpha=alpha)
[ 0.1  1.   1.   0. ]

The above case is for just the suppressor surfactant
concentration. Everything else is switched off. Let's try another
case. Non-physical root is the attractor.

>>> print solve(Bbar=0.1 / E(0, Fbar, alpha), Kbar=0, Cbar=0, Gbar=0, Fbar=Fbar, alpha=alpha)
[ 10.   1.   1.   0.]

Let's include the effects of the potential and vary the position of
the reference electrode (:math:`G`).

>>> import numpy
>>> Garray = -10**numpy.linspace(-3, 3, 1000)
>>> thetas = [solve(Bbar=Bbar, Kbar=0, Cbar=0, Gbar=Gbar, Fbar=Fbar, alpha=alpha)[0] for Gbar in Garray]

>>> import pylab
>>> pylab.figure()
>>> pylab.semilogx(abs(Garray), thetas)
>>> pylab.xlabel(r'$\bar{G}=\frac{\delta_{ref} i_0 E_A}{\kappa \eta}$')
>>> pylab.ylabel(r'$\theta$')
>>> pylab.savefig('Gvtheta.png') 

.. image:: Gvtheta.*
   :width: 90%
   :align: center
   :alt: Stuff

Let's try a more physical system and vary the potential.

>>> A = 7.73e-8
>>> kMinus = 2.45e7
>>> kPlus = 125.
>>> omega = 7.1e-6
>>> n = 2
>>> sinf = .02
>>> i0 = A * n * F / omega
>>> deltaRef = 0.03
>>> kappa = 15.26

>>> Fbar = eta * F / R / T
>>> Bbar = kMinus / kPlus * omega / n / F * i0 / sinf
>>> Gbar = deltaRef * i0 / kappa / eta

The physical solution for a flat interface is close to 1.

>>> print solve(Bbar=Bbar, Kbar=0, Cbar=0, Gbar=Gbar, Fbar=Fbar, alpha=alpha)
[ 0.95669784  1.          1.          0.89737542]

If we allow some suppressor diffusion.

>>> gamma = 2.5e-7
>>> delta = 150e-6
>>> Ds = 9.2e-11  
>>> Kbar = gamma * kPlus * delta / Ds

>>> print solve(Bbar=Bbar, Kbar=Kbar, Cbar=0, Gbar=Gbar, Fbar=Fbar, alpha=alpha)
[ 0.37914907  0.0306438   1.          0.99485398]

and some cupric diffusion

>>> Dc = 2.65e-10
>>> cinf = 1000.

>>> Cbar = delta * i0 / Dc / n / F / cinf / 100

>>> print solve(Bbar=Bbar, Kbar=Kbar, Cbar=Cbar, Gbar=Gbar, Fbar=Fbar, alpha=alpha)
[ 0.37914907  0.0306438   1.          0.99485398]
##.. image:: etavtheta.*
##   :width: 90%
##   :align: center
##   :alt: Stuff


"""
__docformat__ = 'restructuredtext'


def E(xi, Fbar, alpha):
     import numpy
     return numpy.exp(-alpha * Fbar * (1 - xi)) - numpy.exp((2 - alpha) * Fbar * (1 - xi))

def func(x, Bbar, Kbar, Cbar, Gbar, Fbar, alpha):
    theta, suppressor, cupric, potential = x
    EE = E(potential, Fbar, alpha)
    return (theta - suppressor / Bbar / EE / cupric,
            suppressor - 1 / (1 + Kbar * (1 - theta)),
            cupric - 1 / (1 + Cbar * EE * (1 - theta)),
            potential + Gbar * (1 - theta) * EE * cupric)
    
def solve(theta=0, suppressor=1, cupric=1, potential=0,
          Bbar=None, Kbar=None, Cbar=None, Gbar=None, Fbar=None, alpha=None):
    from scipy.optimize import fsolve
    return fsolve(func,
                  (theta, suppressor, cupric, potential),
                  (Bbar, Kbar, Cbar, Gbar, Fbar, alpha))


if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())



