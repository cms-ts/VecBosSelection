#include "FWCore/Framework/interface/EDFilter.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"
#include <string>
#include <cmath>
#include "TH1.h"
#include <iostream>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"

class ZanalyzerFilter : public edm::EDFilter {
   public:
      explicit ZanalyzerFilter(const edm::ParameterSet &);
      ~ZanalyzerFilter();

       virtual void beginJob();

   private:
      virtual bool filter(edm::Event&, edm::EventSetup const&);
      virtual void endJob() ;

      // ----------member data ---------------------------

edm::InputTag theElectronCollectionLabel;
  edm::InputTag triggerCollectionTag_; 
};

  std::string outputFile_;
TFile *fOFile;
TH1D* eventMultip;
TH1D* eventAccept;
TH1F* gsfelEt;
TH1D* Conversion;
TH1D* Isolation;
TH1D* Identification;
TH1D* Selected;
TH1F* h_invMass;
TH1F* h_invMassEE;
TH1F* h_invMassEB;
TH1F* h_invMassBB;
