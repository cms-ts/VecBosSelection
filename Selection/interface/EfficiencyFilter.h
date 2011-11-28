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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

class EfficiencyFilter : public edm::EDFilter {
   public:
      explicit EfficiencyFilter(const edm::ParameterSet &);
      ~EfficiencyFilter();

       virtual void beginJob();
       virtual bool beginRun(edm::Run &, edm::EventSetup const&);

   private:
      virtual bool filter(edm::Event&, edm::EventSetup const&);
      virtual void endJob() ;
      virtual bool DoWP80(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent);
      virtual bool DoHLTMatch(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent);

      // ----------member data ---------------------------

      edm::InputTag theElectronCollectionLabel;
      edm::InputTag triggerCollection_; 
      edm::InputTag electronIsolatedProducer_;
      edm::InputTag candTag_ ;
      edm::InputTag theJetCollectionLabel_;
      bool useCombinedPrescales_; // switch between HLT only and L1*HLT prescales
      bool useAllTriggers_; // if no trigger names are provided, use all triggers to find event weight
      HLTConfigProvider hltConfig_;        // to get configuration for L1s/Pre
      std::vector<std::string> triggerNames_; // name of the algorithms selected by our analysis
      std::vector<unsigned int> triggerIndices_; // index of the algorithms selected by our analysis
      bool removePU_;

      TH1D *probepass;
      TH1D *probefail;
      TH1D *probeall;
      TH1D *probepass1jet;
      TH1D *probepass2jet;
      TH1D *probepass3jet;
      TH1D *probepass4jet;
      TH1D *probepass0jet;

};

//
// constructors and destructor
//
EfficiencyFilter::EfficiencyFilter (const edm::ParameterSet & parameters)
{
  theElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("electronCollection");
  std::string outputfile_D = parameters.getUntrackedParameter<std::string>("filename");
  triggerCollection_=parameters.getUntrackedParameter<edm::InputTag>("triggerCollectionTag");
  useCombinedPrescales_ = parameters.getParameter<bool>("UseCombinedPrescales");
  triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");
  useAllTriggers_       = (triggerNames_.size()==0);
  removePU_             = parameters.getParameter<bool>("removePU");
  electronIsolatedProducer_ = parameters.getParameter< edm::InputTag > ("electronIsolatedProducer");
  candTag_ = parameters.getParameter< edm::InputTag > ("candTag");
  theJetCollectionLabel_       = parameters.getParameter<edm::InputTag>("JetCollectionLabel");
  triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");



  //Initializations...
  edm::Service<TFileService> fs;
  probefail = fs->make<TH1D>("Probefail","Invariant mass when probe fails", 60, 60.0, 120.0);
  probepass = fs->make<TH1D>("Probepass","Invariant mass when probe passes", 60, 60.0, 120.0);
  probeall = fs->make<TH1D>("Probeall","Invariant mass when probe fails or passes", 60, 60.0, 120.0);
  probepass0jet = fs->make<TH1D>("Probepass0Jet","Invariant mass when probe passes + no Jet", 60, 60.0, 120.0);
  probepass1jet = fs->make<TH1D>("Probepass1Jet","Invariant mass when probe passes + 1 Jet", 60, 60.0, 120.0);
  probepass2jet = fs->make<TH1D>("Probepass2Jet","Invariant mass when probe passes + 2 Jets", 60, 60.0, 120.0);
  probepass3jet = fs->make<TH1D>("Probepass3Jet","Invariant mass when probe passes + 3 Jets", 60, 60.0, 120.0);
  probepass4jet = fs->make<TH1D>("Probepass4Jet","Invariant mass when probe passes + 4 Jets", 60, 60.0, 120.0);


}



EfficiencyFilter::~EfficiencyFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}



