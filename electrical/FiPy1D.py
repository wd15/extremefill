from fipy import Grid1D, DistanceVariable, CellVariable
from fipy import buildMetalIonDiffusionEquation, numerix
from potentialEquation import PotentialEquation

delta = 100e-9
L = delta * 1.1
N = 110
metalDiffusion = 5e-10
bulkMetalConcentration = 250.
molarVolume = 7.1e-6
faradaysConstant = 9.6485e4
charge = 2
gasConstant = 8.314
temperature = 298.
Fbar = faradaysConstant / gasConstant / temperature
alpha = 0.4
appliedVoltage = -0.275

mesh = Grid1D(Lx=L, nx=N) + ((L - delta,),)

distanceVar = DistanceVariable(mesh=mesh)
distanceVar[:] = -1.
distanceVar[mesh.x > 0] = 1.

metalVar = CellVariable(mesh=mesh, hasOld=True)
metalVar[:] = bulkMetalConcentration
metalVar.constrain(bulkMetalConcentration, where=mesh.facesLeft)

potentialVar = CellVariable(mesh=mesh)
potentialVar[:] = appliedVoltage

doubleLayerVoltage = CellVariable(mesh=mesh, hasOld=True)

def exchangeCurrentDensity(VDL):
    return numerix.exp(-alpha * Fbar * VDL) - numerix.exp((2 - alpha) * Fbar * VDL)

workingCurrent = metalVar / bulkMetalConcentration * exchangeCurrentDensity(doubleLayerVoltage)

depositionRate = workingCurrent * molarVolume / charge / faradaysConstant

metalEquation = buildMetalIonDiffusionEquation(ionVar=metalVar,
                                               distanceVar=distanceVar,
                                               depositionRate=depositionRate,
                                               diffusionCoeff=metalDiffusion,
                                               molarVolume=molarVolume)

potentialEqn = PotentialEquation(var=potentialVar,
                                 distanceVar=distanceVar,
                                 diffusionCoeff=kappa,
                                 transientCoeff=0.)

 = diffusionVar._cellInterfaceFlag
potentialEqn += ImplicitSourceTerm(diffusionMask * largeValue) - diffusionMask * (appliedVoltage - 
