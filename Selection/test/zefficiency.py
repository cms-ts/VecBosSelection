import FWCore.ParameterSet.Config as cms

import os

process = cms.Process("ZSelection")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V20::All'

process.load("MagneticField.Engine.uniformMagneticField_cfi") 

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) )

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True),
                                     makeTriggerResults=cms.untracked.bool(True),
                                     )

###########
# HLT Summary
#########

#process.MessageLogger.destinations = ['HLTreport_Mu_All.txt']
#from HLTrigger.HLTanalyzers.hlTrigReport_cfi import hlTrigReport
#process.hltReport = hlTrigReport.clone(
#    HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
#    )

#process.endpath = cms.EndPath(process.hltReport) 
#process.MessageLogger.categories.append("HLTrigReport")



readFiles = cms.untracked.vstring()
readFiles.extend([
    "/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-May10ReReco-v1/0000/0234F556-657C-E011-9556-002618943948.root",
    "/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-May10ReReco-v1/0000/FE8C3F99-D97B-E011-BEA4-0018F3D096EE.root",
    ])

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
     reportEvery = cms.untracked.int32(500),
 )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            #fileNames = readFiles,
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames =cms.untracked.vstring('file:/tmp/marone/0EA02EE5-617C-E011-BC49-00304867926C.root'),
                            )

trigger2011v1  = cms.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3","HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3","HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3","HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3","HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2","HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v3","HLT_Ele45_CaloIdVT_TrkIdT_v3","HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v4","HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4","HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v4")

trigger2010    = cms.vstring("HLT_Ele17_CaloIdl_Ele8_CaloIsoIdL_CaloIsoVL_v3","HLT_Ele15_SW_L1R","HLT_Ele15_SW_CaloEleId_L1R","HLT_Ele17_SW_CaloEleId_L1R","HLT_Ele17_SW_TightEleId_L1R","HLT_Ele17_SW_TightEleId_L1R_v2","HLT_Ele17_SW_TightEleId_L1R_v3","HLT_Photon10_L1R","HLT_Photon15_L1R","HTL_Photon15_Cleaned_L1R")

alltriggers    = cms.vstring() # In this way, the HLT string is empty and it will trigger every event

process.load("JetCollections_cfi")


process.Selection = cms.EDFilter('EfficiencyFilter',
                                 electronCollection = cms.InputTag("gsfElectrons"),
                                 triggerCollectionTag = cms.untracked.InputTag("TriggerResults","","HLT"),
                                 filename=cms.untracked.string("ZAnalysisFilter.root"),
                                 UseCombinedPrescales = cms.bool(True),
                                 TriggerNames = alltriggers,
                                 removePU=  cms.bool(False),
                                 )

process.p = cms.Path(process.PFJetPath
                     *process.Selection
                     )
