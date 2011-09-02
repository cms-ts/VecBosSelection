import FWCore.ParameterSet.Config as cms

##########################################################################
##########################################################################
##########################################################################
##################### Clean PFJets
##-------------------- Import the Jet RECO modules -----------------------
##-------------------- Turn-on the FastJet density calculation -----------------------
from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJets = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True)
kt6PFJets.Rho_EtaMax = cms.double(4.5)

# to compute FastJet rho to correct isolation (note: EtaMax restricted to 2.5)
kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True)
kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

##--------- Turn-on the FastJet jet area calculation for your favorite algorithm ----
#from RecoJets.JetProducers.ak5PFJets_cfi import *
#ak5PFJets.Rho_EtaMax = cms.double(4.5)
#ak5PFJets.doAreaFastjet = True
##########################################################################
PFJetPath = cms.Sequence(
    kt6PFJets +
    kt6PFJetsForIsolation
# ak5PFJets +
# ak5PFJetsClean +
# ak5PFJetsLooseId 
    )
##########################################
