from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from Configuration.Geometry.GeometryIdeal_cff import *
from Configuration.StandardSequences.MagneticField_cff import *
from TrackingTools.GeomPropagators.SmartPropagator_cff import *
from TrackingTools.MaterialEffects.MaterialPropagator_cfi import *
from TrackingTools.MaterialEffects.OppositeMaterialPropagator_cfi import *

import PhysicsTools.PatAlgos.patSequences_cff

muonMatch = PhysicsTools.PatAlgos.patSequences_cff.muonMatch.clone(src = cms.InputTag("muons"),
                                                                   resolveByMatchQuality = cms.bool(True))
patMuons = PhysicsTools.PatAlgos.patSequences_cff.patMuons.clone(muonSource = cms.InputTag("muons"),
                                                                 genParticleMatch = cms.InputTag("muonMatch"),
                                                                 addTeVRefits = cms.bool(False),
                                                                 embedTrack = cms.bool(True),
                                                                 embedCombinedMuon = cms.bool(True),
                                                                 embedStandAloneMuon = cms.bool(True),
                                                                 embedPickyMuon = cms.bool(False),
                                                                 embedTpfmsMuon = cms.bool(False),
                                                                 embedHighLevelSelection = cms.bool(True),
                                                                 usePV = cms.bool(True),
                                                                 beamLineSrc = cms.InputTag("offlineBeamSpot"),
                                                                 pvSrc = cms.InputTag("offlinePrimaryVertices"),
                                                                 isolation = cms.PSet(),
                                                                 isoDeposits = cms.PSet(),
                                                                 embedCaloMETMuonCorrs = cms.bool(False),
                                                                 embedTcMETMuonCorrs = cms.bool(False),
                                                                 )
# Tracker Muons Part
selectedPatTrackerMuons = PhysicsTools.PatAlgos.patSequences_cff.selectedPatMuons.clone(src = cms.InputTag("patMuons"),
                cut = cms.string("pt > 3.5 && isTrackerMuon() && numberOfMatches() > 1 && -2.4 < eta() && eta() < 2.4"))
cleanPatTrackerMuons = PhysicsTools.PatAlgos.patSequences_cff.cleanPatMuons.clone(src = cms.InputTag("selectedPatTrackerMuons"))
countPatTrackerMuons = PhysicsTools.PatAlgos.patSequences_cff.countPatMuons.clone(src = cms.InputTag("cleanPatTrackerMuons"))

# PF Muons Part
selectedPatPFMuons = PhysicsTools.PatAlgos.patSequences_cff.selectedPatMuons.clone(src = cms.InputTag("patMuons"),
    #"Loose Muon" requirement on PF muons as recommended by Muon POG:
    #https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
    cut = cms.string("pt > 3.5 && isPFMuon() && ( isTrackerMuon() || isGlobalMuon() ) && -2.4 < eta() && eta() < 2.4"))
cleanPatPFMuons = PhysicsTools.PatAlgos.patSequences_cff.cleanPatMuons.clone(src = cms.InputTag("selectedPatPFMuons"))
countPatPFMuons = PhysicsTools.PatAlgos.patSequences_cff.countPatMuons.clone(src = cms.InputTag("cleanPatPFMuons"))

from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi import *

# Trigger match
#    First matcher from PhysicsTools/PatAlgos/python/triggerLayer1/triggerMatcher_cfi.py
#    is cleanMuonTriggerMatchHLTMu20 . Clone it!
#    Note in 2012 wildcard HLT_Mu* includes ONLY muon trigger. No more HLT_MultiVertex6 and such!
# This is trigger match for Tracker muons

cleanTrackerMuonTriggerMatchHLTMu = cleanMuonTriggerMatchHLTMu20.clone(src = cms.InputTag( "cleanPatTrackerMuons" ),
                                                               matchedCuts = cms.string('path("HLT_Mu*")'))
cleanTrackerMuonTriggerMatchHLTIsoMu = cleanMuonTriggerMatchHLTMu20.clone(src = cms.InputTag( "cleanPatTrackerMuons" ),
                                                               matchedCuts = cms.string('path("HLT_IsoMu*")'))
cleanTrackerMuonTriggerMatchHLTDoubleMu = cleanMuonTriggerMatchHLTMu20.clone(src = cms.InputTag( "cleanPatTrackerMuons" ),
                                                               matchedCuts = cms.string('path("HLT_DoubleMu*_v*")'))
cleanTrackerMuonTriggerMatchHLTMuOnia = cleanMuonTriggerMatchHLTMu20.clone(src = cms.InputTag( "cleanPatTrackerMuons" ),
                                                               matchedCuts = cms.string('path("HLT_Dimuon0_Jpsi_Muon_v*")'))
cleanPatTrackerMuonsTriggerMatch = cms.EDProducer("PATTriggerMatchMuonEmbedder",
                                                  src = cms.InputTag("cleanPatTrackerMuons"),
                                                  matches = cms.VInputTag("cleanTrackerMuonTriggerMatchHLTMu",
                                                                          "cleanTrackerMuonTriggerMatchHLTIsoMu",
                                                                          "cleanTrackerMuonTriggerMatchHLTDoubleMu",
                                                                          "cleanTrackerMuonTriggerMatchHLTMuOnia"))


#cleanPatTrackerMuonsTriggerMatch = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
#                                                  src = cms.InputTag("cleanPatTrackerMuons"),
#                                                  matched = cms.InputTag( "patTrigger" ),
#                                                  matchedCuts = cms.string( 'path( "HLT_Dimuon0_Jpsi_Muon*" )' ),
#                                                  maxDPtRel = cms.double( 0.5 ),
#                                                  maxDeltaR = cms.double( 0.5 ),
#                                                  resolveAmbiguities    = cms.bool( True ),  # only one match per trigger object
#                                                  resolveByMatchQuality = cms.bool( True ),
#                                                  )


# This is trigger match for PF muons

cleanPFMuonTriggerMatchHLTMu = cleanMuonTriggerMatchHLTMu20.clone(src = cms.InputTag( "cleanPatPFMuons" ),
                                                           matchedCuts = cms.string('path("HLT_Mu*")')) 
cleanPFMuonTriggerMatchHLTIsoMu = cleanMuonTriggerMatchHLTMu20.clone(src = cms.InputTag( "cleanPatPFMuons" ),
                                                           matchedCuts = cms.string('path("HLT_IsoMu*")'))
cleanPFMuonTriggerMatchHLTDoubleMu = cleanMuonTriggerMatchHLTMu20.clone(src = cms.InputTag( "cleanPatPFMuons" ),
                                                          matchedCuts = cms.string('path("HLT_DoubleMu*_v*")'))
cleanPFMuonTriggerMatchHLTMuOnia = cleanMuonTriggerMatchHLTMu20.clone(src = cms.InputTag( "cleanPatPFMuons" ),
                                                          matchedCuts = cms.string('path("HLT_Dimuon0_Jpsi_Muon_v*")'))
cleanPatPFMuonsTriggerMatch = cms.EDProducer("PATTriggerMatchMuonEmbedder",
                                             src = cms.InputTag("cleanPatPFMuons"),
                                             matches = cms.VInputTag("cleanPFMuonTriggerMatchHLTMu",
                                                                     "cleanPFMuonTriggerMatchHLTIsoMu",
                                                                     "cleanPFMuonTriggerMatchHLTDoubleMu",
                                                                     "cleanPFMuonTriggerMatchHLTMuOnia"))


#cleanPatPFMuonsTriggerMatch = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
#                                             src = cms.InputTag("cleanPatPFMuons"),
#                                             matched = cms.InputTag( "patTrigger" ),
#                                             matchedCuts = cms.string( 'path( "HLT_Dimuon0_Jpsi_Muon*" )' ),
#                                             maxDPtRel = cms.double( 0.5 ),
#                                             maxDeltaR = cms.double( 0.5 ),
#                                             resolveAmbiguities    = cms.bool( True ),  # only one match per trigger object
#                                             resolveByMatchQuality = cms.bool( True ),
#                                             )


                                                                   
patifyTrackerMuon = cms.Sequence(selectedPatTrackerMuons * cleanPatTrackerMuons * countPatTrackerMuons * cleanTrackerMuonTriggerMatchHLTMu * cleanTrackerMuonTriggerMatchHLTIsoMu * cleanTrackerMuonTriggerMatchHLTDoubleMu * cleanTrackerMuonTriggerMatchHLTMuOnia * cleanPatTrackerMuonsTriggerMatch)
patifyPFMuon = cms.Sequence(selectedPatPFMuons * cleanPatPFMuons * countPatPFMuons * cleanPFMuonTriggerMatchHLTMu * cleanPFMuonTriggerMatchHLTIsoMu * cleanPFMuonTriggerMatchHLTDoubleMu * cleanPFMuonTriggerMatchHLTMuOnia * cleanPatPFMuonsTriggerMatch)

#patifyTrackerMuon = cms.Sequence(selectedPatTrackerMuons * cleanPatTrackerMuons * countPatTrackerMuons * cleanPatTrackerMuonsTriggerMatch)
#patifyPFMuon = cms.Sequence(selectedPatPFMuons * cleanPatPFMuons * countPatPFMuons * cleanPatPFMuonsTriggerMatch)


patifyMC = cms.Sequence(muonMatch * patMuons * patTrigger * patTriggerEvent * patifyTrackerMuon * patifyPFMuon)
patifyData = cms.Sequence(patMuons * patTrigger * patTriggerEvent * patifyTrackerMuon * patifyPFMuon)
