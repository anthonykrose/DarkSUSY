import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)


readFiles.extend( [

"/store/data/Run2012C/MuOnia/AOD/24Aug2012-v1/00000/007177CD-97FB-E111-B68A-003048FFCC2C.root",
"/store/data/Run2012C/MuOnia/AOD/24Aug2012-v1/00000/02599017-50FB-E111-973E-001A9281172C.root"


]);


readFiles.extend( [

])


secFiles.extend( [
    ] )
