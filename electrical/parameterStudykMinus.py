from feature1D import feature

A = 7.73e-8
n = 2
F = 9.6485e4
omega = 7.1e-6
i0 = A * n * F / omega
i0 = 40.
kMinus = 2.45e7
kPlus = 150.
bulkSuppressor = 0.02 
eta = -0.25
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

for kPlus in (0.01, 0.1, 1., 5., 10., 20., 40., 100., 200.):
    print 'kPlus',kPlus
    feature(featureDepth=56e-6,
            view=False,
            dt=1e-8,
            dtMax=1e+20,
            totalSteps=400,
            PRINT=False,
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
            delta=delta,
            filename='tmp/base-kplus-' + str(kPlus) + '.gz')


