from extremefill.kPlusViewer import KPlusViewer
from extremefill.kMinusViewer import KMinusViewer
from extremefill.appliedPotentialViewer import AppliedPotentialViewer
from extremefill.deltaRefViewer import DeltaRefViewer
from extremefill.featureDepthViewer import FeatureDepthViewer
from extremefill.kPlusVkMinusViewer import KPlusVkMinusViewer
from extremefill.appliedPotentialVbulkSuppressorViewer import AppliedPotentialVbulkSuppressorViewer
from extremefill.bulkSuppressorViewer import BulkSuppressorViewer

## figure ordering from the paper
viewers = (None,
           None,
           KPlusViewer,
           KMinusViewer,
           KPlusVkMinusViewer,
           FeatureDepthViewer,
           DeltaRefViewer,
           AppliedPotentialViewer,
           BulkSuppressorViewer,
           AppliedPotentialVbulkSuppressorViewer)

def generateFigures(filesuffix=('.png',)):
    for Viewer in viewers:
        if Viewer is not None:
            Viewer().plot(filesuffix=filesuffix)

def generateFigure(number, filesuffix=('.png',)):
    Viewer = viewers[number]
    if Viewer is not None:
        Viewer().plot(filesuffix=filesuffix)
           
           

