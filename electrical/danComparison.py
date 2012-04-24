from feature1D import feature
from flatSteadyState1D import solve, func, thetaFunc

times, potentials = feature(view=True,
                            dt=1e-8,
                            dtMax=1e+20,
                            totalSteps=300,
                            PRINT=True,
                            kMinus=1e+7
                            )


