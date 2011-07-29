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

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
     reportEvery = cms.untracked.int32(500),
 )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            #fileNames = readFiles,
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames =cms.untracked.vstring('file:/tmp/marone/0EA02EE5-617C-E011-BC49-00304867926C.root'),
                            )

process.Selection = cms.EDFilter('ZanalyzerFilter',
                                 electronCollection = cms.InputTag("gsfElectrons"),
                                 triggerCollectionTag = cms.untracked.InputTag("TriggerResults","","HLT"),
                                 filename=cms.untracked.string("ZAnalysisFilter.root"),

                                 )


process.demo = cms.EDAnalyzer('HistoAnalyzer'
                              #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
                              , tracks = cms.InputTag('generalTracks')
                              # immagino di dover mettere qualcosa tipo:
                              , electronCollection = cms.InputTag('gsfElectrons')
                              , triggerCollection = cms.InputTag("TriggerResults","","HLT"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo.root')
)



process.p = cms.Path(process.Selection
#                    *process.demo
                     )

process.source.fileNames = ['/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-May10ReReco-v1/0000/402BCC31-C47C-E011-8DBA-003048679048.root']
