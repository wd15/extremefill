from feature1D import feature
from flatSteadyState1D import solve, func, thetaFunc

A = 7.73e-8
n = 2
F = 9.6485e4
omega = 7.1e-6
i0 = A * n * F / omega
kMinus = 2.45e7
kPlus = 125.
bulkSuppressor = 0.02
eta = -0.275
R = 8.314
T = 298.
deltaRef = 0.03
kappa = 15.26
gamma = 2.5e-7
delta = 150e-6
alpha = 0.4
diffusionCupric = 2.65e-10
diffusionSuppressor = 9.2e-11
bulkCupric=1000.

Fbar = F * eta / R / T
Bbar = kMinus * omega * i0 / kPlus / n / F / bulkSuppressor
Gbar = deltaRef * i0 / kappa / eta
Kbar = gamma * kPlus * delta / diffusionSuppressor
Cbar = delta * i0 / diffusionCupric / n / F / bulkCupric

##print solve(Bbar=Bbar, Kbar=0, Cbar=0, Gbar=0, Fbar=Fbar, alpha=alpha)
##print solve(Bbar=Bbar, Kbar=0, Cbar=0, Gbar=Gbar, Fbar=Fbar, alpha=alpha)
##print solve(Bbar=Bbar, Kbar=Kbar, Cbar=0, Gbar=Gbar, Fbar=Fbar, alpha=alpha)
x = solve(Bbar=Bbar, Kbar=Kbar, Cbar=Cbar, Gbar=Gbar, Fbar=Fbar, alpha=alpha)
print func(x[1:], Bbar, Kbar, Cbar, Gbar, Fbar, alpha)
print x
print func((0.220681744779, 0.593413596545, 1 - 0.00786830954277), Bbar, Kbar, Cbar, Gbar, Fbar, alpha)
print thetaFunc(0.220681744779,  Bbar,  1 - 0.00786830954277, Fbar, alpha, 0.593413596545)

raw_input('stopped')

times, potentials = feature(featureDepth=0.,
                            view=True,
                            dt=1e-8,
                            dtMax=10.,
                            totalSteps=10000000,
                            PRINT=True,
                            relaxation=1.0,
                            bulkSuppressor=bulkSuppressor,
                            diffusionCupric=diffusionCupric,
                            diffusionSuppressor=diffusionSuppressor,
                            kappa=kappa,
                            i0=i0,
                            i1=-i0,
                            faradaysConstant=F,
                            charge=n,
                            omega=omega,
                            kPlus=kPlus,
                            kMinus=kMinus,
                            alpha=alpha,
                            appliedPotential=eta,
                            gasConstant=R,
                            temperature=T,
                            deltaRef=deltaRef,
                            gamma=gamma,
                            delta=delta)


