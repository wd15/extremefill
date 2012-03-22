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


from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
from fipy.models.levelSet.distanceFunction.levelSetDiffusionEquation import _buildLevelSetDiffusionEquation
from fipy.models.levelSet.electroChem.metalIonSourceVariable import _MetalIonSourceVariable

__all__ = ["potentialEquation"]

def potentialEquation(var = None,
                      distanceVar = None,
                      currentDensity = None,
                      conductivity=None,
                      capacitance=None)

    r"""

    The `potentialEquation` function generates an equation that models the
    potential distribution and double layer voltage build up at the zero level
    set. The governing equation is given by,

    .. math::

       \nabla \cdot \nabla  \psi = 0

    The boundary condition at $\phi=0$ is given by,

    .. math::

       c_{DL} \frac{\partial \psi}{\partial t} + i_{F} \left( \psi \right)  = 
       -\kappa \vec{n} \cdot \nabla \psi

    where the normal points out of the electrolyte.

    Using the a level set implementation the boundary condition is included in the main equation
       
    \int \left[ |nabla \phi| \delta \left( \phi \right) \frac{\partial \psi }{\partial t} \right] dV = 
    \int \nabla \cdot  \left(\kappa_{\phi} \nabla\right)  \psi dV - \int i_F \left( \psi \right) dV

    where,

    .. math::

       \kappa_{\phi} = \begin{cases}
           \kappa & \text{when $\phi > 0$} \\
           0  & \text{when $\phi \le 0$}
       \end{cases}

    :Parameters:
      - `var`: The potential variable.
      - `distanceVar`: A `DistanceVariable` object.
      - `currentDensity`: A float or a `CellVariable` representing the current Density.
      - `conductivity`: Float
      - `capictance`: Float

    """

    # eq = _buildLevelSetDiffusionEquation(ionVar = ionVar,
    #                                      distanceVar = distanceVar,
    #                                      transientCoeff = transientCoeff,
    #                                      diffusionCoeff = diffusionCoeff)
    
    # coeff = _MetalIonSourceVariable(ionVar = ionVar,
    #                                 distanceVar = distanceVar,
    #                                 depositionRate = depositionRate,
    #                                 metalIonMolarVolume = metalIonMolarVolume)

    # return eq + ImplicitSourceTerm(coeff)

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
