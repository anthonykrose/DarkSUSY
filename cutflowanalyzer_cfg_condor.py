import FWCore.ParameterSet.Config as cms
import os

PROCESS = str(os.getenv("PROCESS"))
if ( PROCESS == None or PROCESS == "" ): PROCESS = int(0)
else : PROCESS = int(PROCESS)
JOBS = 100

process = cms.Process("Analysis")

#process.load("FileLists.PAT.DoubleJPsiDPSto4mu_pTJPsi3GeV_8TeV-pythia8_537p4_PATv3_cff")
#process.load("DoubleJPsi_8TeV_4muon_1hpt_Filelist2_cff")

#process.load("AnalysisDataFormats.MuJetAnalysis.filelists.SPS_cff")


process.load("AnalysisDataFormats.MuJetAnalysis.filelists.mGamma_0250_ctau02_cff");



lenFileNames = len(process.source.fileNames)

process.source.fileNames = process.source.fileNames[lenFileNames*PROCESS/JOBS:lenFileNames*(PROCESS+1)/JOBS]

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger", destinations = cms.untracked.vstring("cout"), cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR")))

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "FT_53_V6_AN3::All"

process.load("AnalysisAlgos.MuJetProducer.MuJetProducer_cff")
process.MuJetProducer.muons = cms.InputTag("cleanPatPFMuonsTriggerMatch")

process.MuJetProducer05 = process.MuJetProducer.clone(maxMass = cms.double(5.))
process.MuJetProducer05.groupingMode = cms.string("GroupByMassAndVertexProbOrDeltaR")

# For PF muons
process.MuJetProducer05.selectTrackerMuons = cms.bool(True)
process.MuJetProducer05.selectGlobalMuons = cms.bool(False)
process.MuJetProducer05.minTrackerHits = cms.int32(-1)
process.MuJetProducer05.maxTrackerNormChi2 = cms.double(-1)
process.MuJetProducer05.minSegmentMatches = cms.int32(-1)
process.MuJetProducer05.maxDeltaR = cms.double(0.01)
process.MuJetProducer05.minPt = cms.double(8.)
process.MuJetProducer05.maxAbsEta = cms.double(2.4)

process.Analysis = cms.EDAnalyzer('Analysis',
  genParticles = cms.InputTag("genParticles"),
#  muons = cms.InputTag("cleanPatTrackerMuonsTriggerMatch"),
  muons = cms.InputTag("cleanPatPFMuonsTriggerMatch"),
  
  muJets = cms.InputTag("MuJetProducer05"),
#  muJetOrphans = cms.InputTag("MuJetProducer05", "Orphans"),
#  muJetPlusTracks = cms.InputTag("MuJetPlusTracks15"),

)

process.p = cms.Path(process.MuJetProducer05 * process.Analysis)
#process.p = cms.Path(process.Analysis)
#process.p = cms.Path(process.MuJetProducer05)

process.TFileService = cms.Service("TFileService", fileName = cms.string("subJob_mGamma_0250_ctau02_%03d.root" % PROCESS) )

