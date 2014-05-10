import FWCore.ParameterSet.Config as cms

MuJetPlusFakeMuon = cms.EDProducer("MuJetPlusFakeMuon",
                                   numberOfFakes = cms.uint32(1),
                                   muJets = cms.InputTag("MuJetProducer"),
                                   muons = cms.InputTag("cleanPatMuons"),
                                   tracks = cms.InputTag("generalTracks"),
                                   caloTowers = cms.InputTag("towerMaker"),
                                   
                                   minPt = cms.double(3.),
                                   minPmag = cms.double(0.),
                                   maxAbsEta = cms.double(2.4),
                                   minTrackerHits = cms.int32(8),
                                   maxTrackerNormChi2 = cms.double(4.),
                                   maxTrackerDxy = cms.double(-1.),
                                   maxTrackerDz = cms.double(-1.),
                                   maxQoverpError = cms.double(-1.),
                                   maxPhiError = cms.double(-1.),
                                   maxEtaError = cms.double(-1.),
                                   maxDxyError = cms.double(-1.),
                                   maxDzError = cms.double(-1.),
                                   
                                   calculateVertex = cms.bool(True),
                                   calculateIsolation = cms.bool(True),
                                   
                                   groupingMode = cms.string("GroupByMassAndVertexProbOrDeltaR"),
                                   maxDeltaR = cms.double(0.01),
                                   maxMass = cms.double(5.),
                                   minVertexProb = cms.double(0.01),
                                   groupByCharge = cms.string("opposite"),
                                   )
