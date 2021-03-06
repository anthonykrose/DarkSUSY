import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [

"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_000.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_001.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_002.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_003.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_004.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_005.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_006.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_007.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_008.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_009.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_010.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_011.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_012.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_013.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_014.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_015.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_016.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_017.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_018.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_019.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_020.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_021.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_022.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_023.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_024.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_025.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_026.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_027.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_028.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_029.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_030.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_031.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_032.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_033.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_034.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_035.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_036.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_037.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_038.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_039.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_040.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_041.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_042.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_043.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_044.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_045.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_046.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_047.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_048.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_049.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_050.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_051.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_052.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_053.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_054.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_055.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_056.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_057.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_058.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_059.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_060.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_061.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_062.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_063.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_064.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_065.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_066.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_067.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_068.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_069.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_070.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_071.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_072.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_073.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_074.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_075.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_076.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_077.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_078.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_079.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_080.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_081.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_082.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_083.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_084.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_085.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_086.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_087.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_088.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_089.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_090.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_091.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_092.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_093.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_094.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_095.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_096.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_097.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_098.root",
"file:/fdata/hepx/store/user/castaned/DarkSUSY_mH_125_mGammaD_0700_ctauExp_5_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/Ntup_099.root",


   ]);

secFiles.extend( [
    ] )

