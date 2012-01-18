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
#include "TTree.h"
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
#include "../interface/SelectionUtils.h"

class EfficiencyFilter : public edm::EDFilter, public SelectionUtils {
   public:
      explicit EfficiencyFilter(const edm::ParameterSet &);
      ~EfficiencyFilter();

       virtual void beginJob();
       virtual bool beginRun(edm::Run &, edm::EventSetup const&);

   private:
      virtual bool filter(edm::Event&, edm::EventSetup const&);
      virtual void endJob() ;

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

      TH1D *probeall;

      TH1D *WP80_probepass;
      TH1D *WP80_probefail;
      TH1D *HLT_probepass;
      TH1D *HLT_probefail;

      TH1D *WP80_probepass1jet;
      TH1D *WP80_probepass2jet;
      TH1D *WP80_probepass3jet;
      TH1D *WP80_probepass4jet;
      TH1D *WP80_probepass0jet;
      TH1D *WP80_probefail0jet ;
      TH1D *WP80_probefail1jet ;
      TH1D *WP80_probefail2jet ;
      TH1D *WP80_probefail3jet ;
      TH1D *WP80_probefail4jet; 

      TH1D *HLT_probepass1jet;
      TH1D *HLT_probepass2jet;
      TH1D *HLT_probepass3jet;
      TH1D *HLT_probepass4jet;
      TH1D *HLT_probepass0jet;
      TH1D *HLT_probefail0jet ;
      TH1D *HLT_probefail1jet ;
      TH1D *HLT_probefail2jet ;
      TH1D *HLT_probefail3jet ;
      TH1D *HLT_probefail4jet; 

      TH1D *tagall;

      TH1D *WP80_tagpass;
      TH1D *WP80_tagfail;
      TH1D *HLT_tagpass;
      TH1D *HLT_tagfail;

      TH1D *WP80_tagpass1jet;
      TH1D *WP80_tagpass2jet;
      TH1D *WP80_tagpass3jet;
      TH1D *WP80_tagpass4jet;
      TH1D *WP80_tagpass0jet;
      TH1D *WP80_tagfail0jet ;
      TH1D *WP80_tagfail1jet ;
      TH1D *WP80_tagfail2jet ;
      TH1D *WP80_tagfail3jet ;
      TH1D *WP80_tagfail4jet; 

      TH1D *HLT_tagpass1jet;
      TH1D *HLT_tagpass2jet;
      TH1D *HLT_tagpass3jet;
      TH1D *HLT_tagpass4jet;
      TH1D *HLT_tagpass0jet;
      TH1D *HLT_tagfail0jet ;
      TH1D *HLT_tagfail1jet ;
      TH1D *HLT_tagfail2jet ;
      TH1D *HLT_tagfail3jet ;
      TH1D *HLT_tagfail4jet; 

      TH1D *probept;
      TH1D *probept_passWP80;
      TH1D *probept_failWP80;
      TH1D *tagpt;
      TH1D *tagpt_passWP80;
      TH1D *tagpt_failWP80;

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

  probeall  = fs->make<TH1D>("probeall","Invariant mass when probe fails or passes", 60, 60.0, 120.0);

  WP80_probefail = fs->make<TH1D>("WP80_probefail","Invariant mass when probe fails", 60, 60.0, 120.0);
  WP80_probepass = fs->make<TH1D>("WP80_probepass","Invariant mass when probe passes", 60, 60.0, 120.0);
  HLT_probefail = fs->make<TH1D>("HLT_probefail","Invariant mass when probe fails", 60, 60.0, 120.0);
  HLT_probepass = fs->make<TH1D>("HLT_probepass","Invariant mass when probe passes", 60, 60.0, 120.0);

  WP80_probepass0jet = fs->make<TH1D>("WP80_probepass0Jet","Invariant mass when probe passes + no Jet", 60, 60.0, 120.0);
  WP80_probepass1jet = fs->make<TH1D>("WP80_probepass1Jet","Invariant mass when probe passes + 1 Jet", 60, 60.0, 120.0);
  WP80_probepass2jet = fs->make<TH1D>("WP80_probepass2Jet","Invariant mass when probe passes + 2 Jets", 60, 60.0, 120.0);
  WP80_probepass3jet = fs->make<TH1D>("WP80_probepass3Jet","Invariant mass when probe passes + 3 Jets", 60, 60.0, 120.0);
  WP80_probepass4jet = fs->make<TH1D>("WP80_probepass4Jet","Invariant mass when probe passes + 4 Jets", 60, 60.0, 120.0);

  WP80_probefail0jet = fs->make<TH1D>("WP80_probefail0Jet","Invariant mass when probe fail + no Jet", 60, 60.0, 120.0);
  WP80_probefail1jet = fs->make<TH1D>("WP80_probefail1Jet","Invariant mass when probe fail + 1 Jet", 60, 60.0, 120.0);
  WP80_probefail2jet = fs->make<TH1D>("WP80_probefail2Jet","Invariant mass when probe fail + 2 Jets", 60, 60.0, 120.0);
  WP80_probefail3jet = fs->make<TH1D>("WP80_probefail3Jet","Invariant mass when probe fail + 3 Jets", 60, 60.0, 120.0);
  WP80_probefail4jet = fs->make<TH1D>("WP80_probefail4Jet","Invariant mass when probe fail + 4 Jets", 60, 60.0, 120.0);

  HLT_probepass0jet = fs->make<TH1D>("HLT_probepass0Jet","Invariant mass when probe passes + no Jet", 60, 60.0, 120.0);
  HLT_probepass1jet = fs->make<TH1D>("HLT_probepass1Jet","Invariant mass when probe passes + 1 Jet", 60, 60.0, 120.0);
  HLT_probepass2jet = fs->make<TH1D>("HLT_probepass2Jet","Invariant mass when probe passes + 2 Jets", 60, 60.0, 120.0);
  HLT_probepass3jet = fs->make<TH1D>("HLT_probepass3Jet","Invariant mass when probe passes + 3 Jets", 60, 60.0, 120.0);
  HLT_probepass4jet = fs->make<TH1D>("HLT_probepass4Jet","Invariant mass when probe passes + 4 Jets", 60, 60.0, 120.0);

  HLT_probefail0jet = fs->make<TH1D>("HLT_probefail0Jet","Invariant mass when probe fail + no Jet", 60, 60.0, 120.0);
  HLT_probefail1jet = fs->make<TH1D>("HLT_probefail1Jet","Invariant mass when probe fail + 1 Jet", 60, 60.0, 120.0);
  HLT_probefail2jet = fs->make<TH1D>("HLT_probefail2Jet","Invariant mass when probe fail + 2 Jets", 60, 60.0, 120.0);
  HLT_probefail3jet = fs->make<TH1D>("HLT_probefail3Jet","Invariant mass when probe fail + 3 Jets", 60, 60.0, 120.0);
  HLT_probefail4jet = fs->make<TH1D>("HLT_probefail4Jet","Invariant mass when probe fail + 4 Jets", 60, 60.0, 120.0);

  tagall  = fs->make<TH1D>("tagall","Invariant mass when tag fails or passes", 60, 60.0, 120.0);

  WP80_tagfail = fs->make<TH1D>("WP80_tagfail","Invariant mass when tag fails WP80", 60, 60.0, 120.0);
  WP80_tagpass = fs->make<TH1D>("WP80_tagpass","Invariant mass when tag passes WP80", 60, 60.0, 120.0);
  HLT_tagfail = fs->make<TH1D>("HLT_tagfail","Invariant mass when tag fails HLT", 60, 60.0, 120.0);
  HLT_tagpass = fs->make<TH1D>("HLT_tagpass","Invariant mass when tag passes HLT", 60, 60.0, 120.0);

  WP80_tagpass0jet = fs->make<TH1D>("WP80_tagpass0Jet","Invariant mass when tag passes + no Jet", 60, 60.0, 120.0);
  WP80_tagpass1jet = fs->make<TH1D>("WP80_tagpass1Jet","Invariant mass when tag passes + 1 Jet", 60, 60.0, 120.0);
  WP80_tagpass2jet = fs->make<TH1D>("WP80_tagpass2Jet","Invariant mass when tag passes + 2 Jets", 60, 60.0, 120.0);
  WP80_tagpass3jet = fs->make<TH1D>("WP80_tagpass3Jet","Invariant mass when tag passes + 3 Jets", 60, 60.0, 120.0);
  WP80_tagpass4jet = fs->make<TH1D>("WP80_tagpass4Jet","Invariant mass when tag passes + 4 Jets", 60, 60.0, 120.0);

  WP80_tagfail0jet = fs->make<TH1D>("WP80_tagfail0Jet","Invariant mass when tag fail + no Jet", 60, 60.0, 120.0);
  WP80_tagfail1jet = fs->make<TH1D>("WP80_tagfail1Jet","Invariant mass when tag fail + 1 Jet", 60, 60.0, 120.0);
  WP80_tagfail2jet = fs->make<TH1D>("WP80_tagfail2Jet","Invariant mass when tag fail + 2 Jets", 60, 60.0, 120.0);
  WP80_tagfail3jet = fs->make<TH1D>("WP80_tagfail3Jet","Invariant mass when tag fail + 3 Jets", 60, 60.0, 120.0);
  WP80_tagfail4jet = fs->make<TH1D>("WP80_tagfail4Jet","Invariant mass when tag fail + 4 Jets", 60, 60.0, 120.0);

  HLT_tagpass0jet = fs->make<TH1D>("HLT_tagpass0Jet","Invariant mass when tag passes + no Jet", 60, 60.0, 120.0);
  HLT_tagpass1jet = fs->make<TH1D>("HLT_tagpass1Jet","Invariant mass when tag passes + 1 Jet", 60, 60.0, 120.0);
  HLT_tagpass2jet = fs->make<TH1D>("HLT_tagpass2Jet","Invariant mass when tag passes + 2 Jets", 60, 60.0, 120.0);
  HLT_tagpass3jet = fs->make<TH1D>("HLT_tagpass3Jet","Invariant mass when tag passes + 3 Jets", 60, 60.0, 120.0);
  HLT_tagpass4jet = fs->make<TH1D>("HLT_tagpass4Jet","Invariant mass when tag passes + 4 Jets", 60, 60.0, 120.0);

  HLT_tagfail0jet = fs->make<TH1D>("HLT_tagfail0Jet","Invariant mass when tag fail + no Jet", 60, 60.0, 120.0);
  HLT_tagfail1jet = fs->make<TH1D>("HLT_tagfail1Jet","Invariant mass when tag fail + 1 Jet", 60, 60.0, 120.0);
  HLT_tagfail2jet = fs->make<TH1D>("HLT_tagfail2Jet","Invariant mass when tag fail + 2 Jets", 60, 60.0, 120.0);
  HLT_tagfail3jet = fs->make<TH1D>("HLT_tagfail3Jet","Invariant mass when tag fail + 3 Jets", 60, 60.0, 120.0);
  HLT_tagfail4jet = fs->make<TH1D>("HLT_tagfail4Jet","Invariant mass when tag fail + 4 Jets", 60, 60.0, 120.0);

  probept= fs->make<TH1D>("probept","Pt of the electron probe", 120, 0, 120.0);
  probept_passWP80= fs->make<TH1D>("probept_passWP80","Pt of the electron probe when passing WP80", 120, 0, 120.0);
  probept_failWP80= fs->make<TH1D>("probept_failWP80","Pt of the electron probe when failing WP80", 120, 0, 120.0);
  tagpt= fs->make<TH1D>("tagpt","Pt of the electron tag", 120, 0, 120.0);
  tagpt_passWP80= fs->make<TH1D>("tagpt_passWP80","Pt of the electron tag when passing WP80", 120, 0, 120.0);
  tagpt_failWP80= fs->make<TH1D>("tagpt_failWP80","Pt of the electron tag when failing WP80", 120, 0, 120.0);

  HLTnumberOfMatches_PASS    = fs->make<TH1D>("HLTnumberOfMatches_PASS","Total # of electron matching the HLT",100,0,100);
  HLTnumberOfMatches_FAIL    = fs->make<TH1D>("HLTnumberOfMatches_FAIL","Total # of electron NOT matching the HLT",100,0,100);
  HLTnumberOfMatches_TOTALELE= fs->make<TH1D>("HLTnumberOfMatches_TOTALELE","Total # of electrons per event",100,0,100);
}



EfficiencyFilter::~EfficiencyFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}



