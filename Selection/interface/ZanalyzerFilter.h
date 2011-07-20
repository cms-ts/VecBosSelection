#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TFile.h"
#include "TH1.h"
#include <iostream>

class ZanalyzerFilter : public edm::EDFilter {
   public:
      explicit ZanalyzerFilter(const edm::ParameterSet &);
      ~ZanalyzerFilter();

       virtual void beginJob();

   private:
       virtual bool filter(edm::Event&, edm::EventSetup const&);
       virtual void endJob() ;
       std::string outputFile_;
      // ----------member data ---------------------------

edm::InputTag theElectronCollectionLabel;

};

