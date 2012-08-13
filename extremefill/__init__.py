from extremefill.kPlusViewer import KPlusViewer
from extremefill.kMinusViewer import KMinusViewer
from extremefill.appliedPotentialViewer import AppliedPotentialViewer
from extremefill.deltaRefViewer import DeltaRefViewer
from extremefill.featureDepthViewer import FeatureDepthViewer
from extremefill.kPlusVkMinusViewer import KPlusVkMinusViewer
from extremefill.appliedPotentialVbulkSuppressorViewer import AppliedPotentialVbulkSuppressorViewer
from extremefill.bulkSuppressorViewer import BulkSuppressorViewer
from extremefill.schematicViewer import SchematicViewer
from extremefill import simulation

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

def generateFigures(filesuffix=('.png',), datafile='data.h5'):
    for Viewer in viewers:
        if Viewer is not None:
            Viewer(datafile=datafile).plot(filesuffix=filesuffix)

def generateFigure(number, filesuffix=('.png',), datafile='data.h5'):
    Viewer = viewers[number - 1]
    if Viewer is not None:
        Viewer(datafile=datafile).plot(filesuffix=filesuffix)

def test():
    print 'hello'
           
def _getVersion():
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        version = get_distribution(__name__).version
    except DistributionNotFound:
        version = "unknown, try running `python setup.py egg_info`"
        
    return version
    
__version__ = _getVersion()           

