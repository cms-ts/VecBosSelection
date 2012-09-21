// -*- C++ -*-

// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/analyzerMuCuts.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include <string>
#include <cmath>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TString.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// member functions
//

// ------------ method called for each event  ------------
void
analyzerMuCuts::analyze (const edm::Event & iEvent, const edm::EventSetup & iSetup)
{
 
  Handle < pat::MuonCollection > muonCollection;
  iEvent.getByLabel (theMuCollectionLabel, muonCollection);

  if (!muonCollection.isValid ()) {
     
     if (muonCollection->size()<=1) {}
     
     for (pat::MuonCollection::const_iterator recoMu = muonCollection->begin (); recoMu != muonCollection->end (); recoMu++) {
	
	//  INSERIRE LE SELEZIONI
	
     }

  }
}


// ------------ method called once each job just before starting event loop  ------------
void
analyzerMuCuts::beginJob (){


}


// ------------ method called once each job just after ending the event loop  ------------
void
analyzerMuCuts::endJob ()
{

  //endJob

}

DEFINE_FWK_MODULE (analyzerMuCuts);
