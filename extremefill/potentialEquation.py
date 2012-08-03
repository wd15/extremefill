#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "potentialEquation.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##



__docformat__ = 'restructuredtext'


#from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
from fipy import CellVariable
from fipy.models.levelSet.distanceFunction.levelSetDiffusionEquation import _buildLevelSetDiffusionEquation
#from fipy.models.levelSet.electroChem.metalIonSourceVariable import _MetalIonSourceVariable

__all__ = ["PotentialEquation"]

class _InterfaceAreas(CellVariable):
    def __init__(self, distanceVar):
        super(_InterfaceAreas, self).__init__(distanceVar.mesh, hasOld=0)
        self.distanceVar = self._requires(distanceVar)

    def _calcValue(self):
        return self.distanceVar.cellInterfaceAreas.value

def PotentialEquation(var=None,
                      distanceVar=None,
                      currentDensity=None,
                      conductivity=None,
                      capacitance=None):

    r"""
    The `PotentialEquation` function is a factory function for
    generating a FiPy equation that models the electric potential in
    an electrolyte with a specific type of boundary condition. The
    zero level set marks the working electrode / electrolyte
    interface. The governing equation is,

    .. math::

       \nabla^2  \psi = 0

    where :math:`\psi` is the electric potential. The boundary
    condition at :math:`\phi=0` is given by,

    .. math::

       c_{DL} \frac{\partial \psi}{\partial t} + i_{F} \left( \psi \right)  = 
       -\kappa \vec{n} \cdot \nabla \psi

    where :math:`c_{DL}` is the capacitance, :math:`i_F` is the
    current density, :math:`\kappa` is the conductance and the normal
    points out of the electrolyte. In the level set formulation the
    boundary condition can be included as a source term in the bulk
    equation.

    .. math::
       
       \int \left[ |\nabla \phi| \delta \left( \phi \right) \left(
       c_{DL} \frac{\partial \psi }{\partial t} + i_F  \left( \psi
       \right)\right)\right] dV = \int \nabla \cdot \left(\kappa_{\phi}
       \nabla \psi \right)dV

    where,

    .. math::

       \kappa_{\phi} = \begin{cases}
           \kappa & \text{when $\phi > 0$} \\
           0  & \text{when $\phi \le 0$}
       \end{cases}

    Using the approach for modeling surfactants in FiPy, the
    electrical equation can be written,

    .. math::

       \int \left[
       \frac{A_{\phi=0} c_{DL}}{V}
       \frac{\partial \psi }{\partial t} =
       \nabla \cdot \left(\kappa_{\phi}
       \nabla\right)
       - \frac{A_{\phi=0}}{V} i_F
       \right] dV

    where :math:`A_{\phi=0}` is the area of electrolyte / electrode
    interface in the control volume.

    :Parameters:
      - `var`: The potential variable.
      - `distanceVar`: A `DistanceVariable` object.
      - `currentDensity`: A float or a `CellVariable` representing the current Density.
      - `conductivity`: Float
      - `capacitance`: Float

    """
    A = _InterfaceAreas(distanceVar)
    V = var.mesh.cellVolumes

    eq = _buildLevelSetDiffusionEquation(ionVar=var,
                                         distanceVar=distanceVar,
                                         transientCoeff=capacitance * A / V + (distanceVar < 0.),
                                         diffusionCoeff=conductivity)
    
    return eq + currentDensity * A / V

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
