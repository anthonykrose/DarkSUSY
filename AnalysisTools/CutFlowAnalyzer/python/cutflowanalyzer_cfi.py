import FWCore.ParameterSet.Config as cms


# From AR:
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from Configuration.Geometry.GeometryIdeal_cff import *
from Configuration.StandardSequences.MagneticField_cff import *
from TrackingTools.GeomPropagators.SmartPropagator_cff import *

from AnalysisAlgos.MuJetProducer.MuJetProducer_cfi import *

demo = cms.EDAnalyzer('CutFlowAnalyzer'
)
