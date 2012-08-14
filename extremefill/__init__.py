from extremefill.kPlusViewer import KPlusViewer
from extremefill.kMinusViewer import KMinusViewer
from extremefill.appliedPotentialViewer import AppliedPotentialViewer
from extremefill.deltaRefViewer import DeltaRefViewer
from extremefill.featureDepthViewer import FeatureDepthViewer
from extremefill.kPlusVkMinusViewer import KPlusVkMinusViewer
from extremefill.appliedPotentialVbulkSuppressorViewer import AppliedPotentialVbulkSuppressorViewer
from extremefill.bulkSuppressorViewer import BulkSuppressorViewer
from extremefill.schematicViewer import SchematicViewer
import extremefill.simulation
from extremefill.simulation import Simulation

## figure ordering from the paper
viewers = (None,
           SchematicViewer,
           KPlusViewer,
           KMinusViewer,
           KPlusVkMinusViewer,
           FeatureDepthViewer,
           DeltaRefViewer,
           AppliedPotentialViewer,
           BulkSuppressorViewer,
           AppliedPotentialVbulkSuppressorViewer)

def generateFigures(filesuffix=('.png',), datafile='data.h5', fignumbers=(2, 3, 4, 5, 6, 7, 8, 9, 10)):
    r"""
    Generate all the figures in the paper. By default PNG images are
    generated.

    To generate Figure 2 from the paper:

    >>> import extremefill
    >>> extremefill.generateFigures(fignumbers=3)

    .. image:: kPlus.*
       :width: 50%
       :align: center
       :alt: Figure 3 from the paper

    Figure 3 as it appears in the paper.

    :Parameters:
      - `filesuffix`: tuple of the file suffixes of the generated images.
      - `datafile`: path to the cached HDF5 data file to either read from or write to.
      - `fignumbers` : tuple of figure numbers to generate. 
    """
    if fignumbers is int:
        fignumbers = (fignumbers,)

    for number in fignumbers:
        Viewer = viewers[number - 1]
        if Viewer is not None:
            Viewer(datafile=datafile).plot(filesuffix=filesuffix)

def test():
    r"""
    Run all the doctests available.
    """
    import doctest
    doctest.testmod(extremefill.simulation)
       
def run(view=True, **parameters):
    r"""
    Run an individual simulation using
    ``Simulation().run(**parameters)``.

    Run the default simulation for 10 time steps

    >>> import extremefill
    >>> extremefill.run(totalSteps=10)

    .. image:: kPlus10.*
       :width: 50%
       :align: center
       :alt: Default simulation after 10 time steps

    The negative x values represents the trench domain while the
    positive values represent the electrolyte domain. The red, blue,
    green and cyan curves represent normalized values for the cupric
    concentration, adsorbed suppressor, suppressor concentration and
    potential, respectively.

    :Parameters:
      - `view` : view the simulation as it is running
      - `parameters` : any of the parameters used in
        ``Simulation.run``
    """
    Simulation().run(view=view, **parameters)
    raw_input('press key to continue')

def _getVersion():
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        version = get_distribution(__name__).version
    except DistributionNotFound:
        version = "unknown, try running `python setup.py egg_info`"
        
    return version
    
__version__ = _getVersion()           

