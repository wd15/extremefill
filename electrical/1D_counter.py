#!/usr/bin/env python

def integrate(i0=None,
              alpha=None,
              F=None,
              V_APPLIED=None,
              R=None,
              T=None,
              kappa=None,
              delta=None,
              delta_W=None,
              C_DL=None):

    """
    Here we just integrate the 1D ODE that solves the electrical problem.
    See http://matforge.org/wd15/blog/MoffatElectrical for details.
    """

    import numpy

    ## numbers
    Fbar = F / R / T

    def iT(VW, VC):
        return -kappa / delta_W * (V_APPLIED - VW)

    def iF(V):
        return i0 * (numpy.exp(-alpha * Fbar * V) - numpy.exp((2 - alpha) * Fbar * V))

    def iF_V(V):
        return i0 * (-alpha * Fbar * numpy.exp(-alpha * Fbar * V) \
                          - (2 - alpha) * Fbar * numpy.exp((2 - alpha) * Fbar * V))

    def RHS(t, y):
        VW, VC = y
        return numpy.array([(iF(VW) - iT(VW, VC)) / C_DL,
                            (iF(VC) + iT(VW, VC)) / C_DL])

    def jacobian(t, y):
        VW, VC = y
        return numpy.array([[(iF_V(VW) - kappa / delta_W) / C_DL,
                             0],
                            [kappa / delta_W / C_DL,
                             iF_V(VC) / C_DL]])

    from scipy.integrate import ode

    integrator = ode(RHS, jacobian).set_integrator('vode', method='bdf', with_jacobian=True)
    integrator.set_initial_value((0., 0.), 0.)

    times = []
    V_DLs = []
    dt = 1e-8
    while integrator.successful() and integrator.t < 200.0e-7:
        integrator.integrate(integrator.t + dt)
        times.append(integrator.t)
        V_DLs.append(integrator.y)

    VW_DL = V_DLs[-1][0]
    VW_IR = V_APPLIED - VW_DL
    VC_IR = -(delta - delta_W) * VW_IR / delta_W

    return numpy.array((VW_IR, VC_IR)) 


