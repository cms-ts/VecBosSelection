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
      edm::InputTag superClusterCollection_EB_;
      edm::InputTag superClusterCollection_EE_;
      edm::InputTag triggerCollection_; 
      edm::InputTag electronIsolatedProducer_;
      edm::InputTag candTag_ ;
      edm::InputTag theJetCollectionLabel_;
      edm::InputTag VertexCollectionTag_;
      bool useCombinedPrescales_; // switch between HLT only and L1*HLT prescales
      bool useAllTriggers_; // if no trigger names are provided, use all triggers to find event weight
      HLTConfigProvider hltConfig_;        // to get configuration for L1s/Pre
      std::vector<std::string> triggerNames_; // name of the algorithms selected by our analysis
      std::vector<unsigned int> triggerIndices_; // index of the algorithms selected by our analysis
      bool removePU_;

      bool WP80_efficiency_;
      bool HLTele8_efficiency_;
      bool HLTele17_efficiency_;
      bool RECO_efficiency_;
      bool New_HE_;

      TH1D *probeall_mee;
      TH1D *probepass_mee;
      TH1D *probefail_mee;

      TH1D *probeall_leadjetpt;

      TH1D *probeall_pt;
      TH1D *probepass_pt;
      TH1D *probefail_pt;

      TH1D *probeall_eta;
      TH1D *probepass_eta;
      TH1D *probefail_eta;

      TH1D *probepass0jet;
      TH1D *probepass1jet;
      TH1D *probepass2jet;
      TH1D *probepass3jet;
      TH1D *probepass4jet;
      TH1D *probepass5jet;
      TH1D *probefail0jet ;
      TH1D *probefail1jet ;
      TH1D *probefail2jet ;
      TH1D *probefail3jet ;
      TH1D *probefail4jet; 
      TH1D *probefail5jet; 

      TH1D *probepass0pu;
      TH1D *probepass1pu;
      TH1D *probepass2pu;
      TH1D *probepass3pu;
      TH1D *probepass4pu;
      TH1D *probefail0pu ;
      TH1D *probefail1pu ;
      TH1D *probefail2pu ;
      TH1D *probefail3pu ;
      TH1D *probefail4pu; 

      TH1D *probepass0eta;
      TH1D *probepass1eta;
      TH1D *probepass2eta;
      TH1D *probepass3eta;
      TH1D *probepass4eta;
      TH1D *probefail0eta ;
      TH1D *probefail1eta ;
      TH1D *probefail2eta ;
      TH1D *probefail3eta ;
      TH1D *probefail4eta; 

      TH1D *probepass0leadjetpt;
      TH1D *probepass1leadjetpt;
      TH1D *probepass2leadjetpt;
      TH1D *probepass3leadjetpt;
      TH1D *probepass4leadjetpt;
      TH1D *probepass5leadjetpt;
      TH1D *probepass6leadjetpt;
      TH1D *probepass7leadjetpt;
      TH1D *probepass8leadjetpt;
      TH1D *probefail0leadjetpt ;
      TH1D *probefail1leadjetpt ;
      TH1D *probefail2leadjetpt ;
      TH1D *probefail3leadjetpt ;
      TH1D *probefail4leadjetpt; 
      TH1D *probefail5leadjetpt; 
      TH1D *probefail6leadjetpt; 
      TH1D *probefail7leadjetpt; 
      TH1D *probefail8leadjetpt; 

      TH1D *probepass0leadjetpt0nJet;
      TH1D *probepass1leadjetpt0nJet;
      TH1D *probepass2leadjetpt0nJet;
      TH1D *probepass3leadjetpt0nJet;
      TH1D *probepass4leadjetpt0nJet;
      TH1D *probepass5leadjetpt0nJet;
      TH1D *probepass6leadjetpt0nJet;
      TH1D *probepass7leadjetpt0nJet;
      TH1D *probepass8leadjetpt0nJet;
      TH1D *probefail0leadjetpt0nJet;
      TH1D *probefail1leadjetpt0nJet;
      TH1D *probefail2leadjetpt0nJet;
      TH1D *probefail3leadjetpt0nJet;
      TH1D *probefail4leadjetpt0nJet; 
      TH1D *probefail5leadjetpt0nJet; 
      TH1D *probefail6leadjetpt0nJet; 
      TH1D *probefail7leadjetpt0nJet; 
      TH1D *probefail8leadjetpt0nJet; 

      TH1D *probepass0leadjetpt1nJet;
      TH1D *probepass1leadjetpt1nJet;
      TH1D *probepass2leadjetpt1nJet;
      TH1D *probepass3leadjetpt1nJet;
      TH1D *probepass4leadjetpt1nJet;
      TH1D *probepass5leadjetpt1nJet;
      TH1D *probepass6leadjetpt1nJet;
      TH1D *probepass7leadjetpt1nJet;
      TH1D *probepass8leadjetpt1nJet;
      TH1D *probefail0leadjetpt1nJet;
      TH1D *probefail1leadjetpt1nJet;
      TH1D *probefail2leadjetpt1nJet;
      TH1D *probefail3leadjetpt1nJet;
      TH1D *probefail4leadjetpt1nJet; 
      TH1D *probefail5leadjetpt1nJet; 
      TH1D *probefail6leadjetpt1nJet; 
      TH1D *probefail7leadjetpt1nJet; 
      TH1D *probefail8leadjetpt1nJet; 

      TH1D *probepass0leadjetpt2nJet;
      TH1D *probepass1leadjetpt2nJet;
      TH1D *probepass2leadjetpt2nJet;
      TH1D *probepass3leadjetpt2nJet;
      TH1D *probepass4leadjetpt2nJet;
      TH1D *probepass5leadjetpt2nJet;
      TH1D *probepass6leadjetpt2nJet;
      TH1D *probepass7leadjetpt2nJet;
      TH1D *probepass8leadjetpt2nJet;
      TH1D *probefail0leadjetpt2nJet;
      TH1D *probefail1leadjetpt2nJet;
      TH1D *probefail2leadjetpt2nJet;
      TH1D *probefail3leadjetpt2nJet;
      TH1D *probefail4leadjetpt2nJet; 
      TH1D *probefail5leadjetpt2nJet; 
      TH1D *probefail6leadjetpt2nJet; 
      TH1D *probefail7leadjetpt2nJet; 
      TH1D *probefail8leadjetpt2nJet; 

      TH1D *probepass0leadjetpt3nJet;
      TH1D *probepass1leadjetpt3nJet;
      TH1D *probepass2leadjetpt3nJet;
      TH1D *probepass3leadjetpt3nJet;
      TH1D *probepass4leadjetpt3nJet;
      TH1D *probepass5leadjetpt3nJet;
      TH1D *probepass6leadjetpt3nJet;
      TH1D *probepass7leadjetpt3nJet;
      TH1D *probepass8leadjetpt3nJet;
      TH1D *probefail0leadjetpt3nJet;
      TH1D *probefail1leadjetpt3nJet;
      TH1D *probefail2leadjetpt3nJet;
      TH1D *probefail3leadjetpt3nJet;
      TH1D *probefail4leadjetpt3nJet; 
      TH1D *probefail5leadjetpt3nJet; 
      TH1D *probefail6leadjetpt3nJet; 
      TH1D *probefail7leadjetpt3nJet; 
      TH1D *probefail8leadjetpt3nJet; 

      TH1D *probepass0leadjetpt4nJet;
      TH1D *probepass1leadjetpt4nJet;
      TH1D *probepass2leadjetpt4nJet;
      TH1D *probepass3leadjetpt4nJet;
      TH1D *probepass4leadjetpt4nJet;
      TH1D *probepass5leadjetpt4nJet;
      TH1D *probepass6leadjetpt4nJet;
      TH1D *probepass7leadjetpt4nJet;
      TH1D *probepass8leadjetpt4nJet;
      TH1D *probefail0leadjetpt4nJet;
      TH1D *probefail1leadjetpt4nJet;
      TH1D *probefail2leadjetpt4nJet;
      TH1D *probefail3leadjetpt4nJet;
      TH1D *probefail4leadjetpt4nJet; 
      TH1D *probefail5leadjetpt4nJet; 
      TH1D *probefail6leadjetpt4nJet; 
      TH1D *probefail7leadjetpt4nJet; 
      TH1D *probefail8leadjetpt4nJet; 

      TH1D *probepass0leadjetpt5nJet;
      TH1D *probepass1leadjetpt5nJet;
      TH1D *probepass2leadjetpt5nJet;
      TH1D *probepass3leadjetpt5nJet;
      TH1D *probepass4leadjetpt5nJet;
      TH1D *probepass5leadjetpt5nJet;
      TH1D *probepass6leadjetpt5nJet;
      TH1D *probepass7leadjetpt5nJet;
      TH1D *probepass8leadjetpt5nJet;
      TH1D *probefail0leadjetpt5nJet;
      TH1D *probefail1leadjetpt5nJet;
      TH1D *probefail2leadjetpt5nJet;
      TH1D *probefail3leadjetpt5nJet;
      TH1D *probefail4leadjetpt5nJet; 
      TH1D *probefail5leadjetpt5nJet; 
      TH1D *probefail6leadjetpt5nJet; 
      TH1D *probefail7leadjetpt5nJet; 
      TH1D *probefail8leadjetpt5nJet; 

      TH1D *tagall_mee;
      TH1D *tagpass_mee;
      TH1D *tagfail_mee;

      TH1D *tagall_leadjetpt;

      TH1D *tagall_pt;
      TH1D *tagpass_pt;
      TH1D *tagfail_pt;

      TH1D *tagall_eta;
      TH1D *tagpass_eta;
      TH1D *tagfail_eta;

      TH1D *tagpass0jet;
      TH1D *tagpass1jet;
      TH1D *tagpass2jet;
      TH1D *tagpass3jet;
      TH1D *tagpass4jet;
      TH1D *tagpass5jet;
      TH1D *tagfail0jet ;
      TH1D *tagfail1jet ;
      TH1D *tagfail2jet ;
      TH1D *tagfail3jet ;
      TH1D *tagfail4jet; 
      TH1D *tagfail5jet; 

      TH1D *tagpass0pu;
      TH1D *tagpass1pu;
      TH1D *tagpass2pu;
      TH1D *tagpass3pu;
      TH1D *tagpass4pu;
      TH1D *tagfail0pu ;
      TH1D *tagfail1pu ;
      TH1D *tagfail2pu ;
      TH1D *tagfail3pu ;
      TH1D *tagfail4pu; 

      TH1D *tagpass0eta;
      TH1D *tagpass1eta;
      TH1D *tagpass2eta;
      TH1D *tagpass3eta;
      TH1D *tagpass4eta;
      TH1D *tagfail0eta ;
      TH1D *tagfail1eta ;
      TH1D *tagfail2eta ;
      TH1D *tagfail3eta ;
      TH1D *tagfail4eta; 

      TH1D *tagpass0leadjetpt;
      TH1D *tagpass1leadjetpt;
      TH1D *tagpass2leadjetpt;
      TH1D *tagpass3leadjetpt;
      TH1D *tagpass4leadjetpt;
      TH1D *tagpass5leadjetpt;
      TH1D *tagpass6leadjetpt;
      TH1D *tagpass7leadjetpt;
      TH1D *tagpass8leadjetpt;
      TH1D *tagfail0leadjetpt;
      TH1D *tagfail1leadjetpt;
      TH1D *tagfail2leadjetpt;
      TH1D *tagfail3leadjetpt;
      TH1D *tagfail4leadjetpt; 
      TH1D *tagfail5leadjetpt; 
      TH1D *tagfail6leadjetpt; 
      TH1D *tagfail7leadjetpt; 
      TH1D *tagfail8leadjetpt; 

      TH1D *tagpass0leadjetpt0nJet;
      TH1D *tagpass1leadjetpt0nJet;
      TH1D *tagpass2leadjetpt0nJet;
      TH1D *tagpass3leadjetpt0nJet;
      TH1D *tagpass4leadjetpt0nJet;
      TH1D *tagpass5leadjetpt0nJet;
      TH1D *tagpass6leadjetpt0nJet;
      TH1D *tagpass7leadjetpt0nJet;
      TH1D *tagpass8leadjetpt0nJet;
      TH1D *tagfail0leadjetpt0nJet;
      TH1D *tagfail1leadjetpt0nJet;
      TH1D *tagfail2leadjetpt0nJet;
      TH1D *tagfail3leadjetpt0nJet;
      TH1D *tagfail4leadjetpt0nJet; 
      TH1D *tagfail5leadjetpt0nJet; 
      TH1D *tagfail6leadjetpt0nJet; 
      TH1D *tagfail7leadjetpt0nJet; 
      TH1D *tagfail8leadjetpt0nJet; 

      TH1D *tagpass0leadjetpt1nJet;
      TH1D *tagpass1leadjetpt1nJet;
      TH1D *tagpass2leadjetpt1nJet;
      TH1D *tagpass3leadjetpt1nJet;
      TH1D *tagpass4leadjetpt1nJet;
      TH1D *tagpass5leadjetpt1nJet;
      TH1D *tagpass6leadjetpt1nJet;
      TH1D *tagpass7leadjetpt1nJet;
      TH1D *tagpass8leadjetpt1nJet;
      TH1D *tagfail0leadjetpt1nJet;
      TH1D *tagfail1leadjetpt1nJet;
      TH1D *tagfail2leadjetpt1nJet;
      TH1D *tagfail3leadjetpt1nJet;
      TH1D *tagfail4leadjetpt1nJet; 
      TH1D *tagfail5leadjetpt1nJet; 
      TH1D *tagfail6leadjetpt1nJet; 
      TH1D *tagfail7leadjetpt1nJet; 
      TH1D *tagfail8leadjetpt1nJet; 

      TH1D *tagpass0leadjetpt2nJet;
      TH1D *tagpass1leadjetpt2nJet;
      TH1D *tagpass2leadjetpt2nJet;
      TH1D *tagpass3leadjetpt2nJet;
      TH1D *tagpass4leadjetpt2nJet;
      TH1D *tagpass5leadjetpt2nJet;
      TH1D *tagpass6leadjetpt2nJet;
      TH1D *tagpass7leadjetpt2nJet;
      TH1D *tagpass8leadjetpt2nJet;
      TH1D *tagfail0leadjetpt2nJet;
      TH1D *tagfail1leadjetpt2nJet;
      TH1D *tagfail2leadjetpt2nJet;
      TH1D *tagfail3leadjetpt2nJet;
      TH1D *tagfail4leadjetpt2nJet; 
      TH1D *tagfail5leadjetpt2nJet; 
      TH1D *tagfail6leadjetpt2nJet; 
      TH1D *tagfail7leadjetpt2nJet; 
      TH1D *tagfail8leadjetpt2nJet; 

      TH1D *tagpass0leadjetpt3nJet;
      TH1D *tagpass1leadjetpt3nJet;
      TH1D *tagpass2leadjetpt3nJet;
      TH1D *tagpass3leadjetpt3nJet;
      TH1D *tagpass4leadjetpt3nJet;
      TH1D *tagpass5leadjetpt3nJet;
      TH1D *tagpass6leadjetpt3nJet;
      TH1D *tagpass7leadjetpt3nJet;
      TH1D *tagpass8leadjetpt3nJet;
      TH1D *tagfail0leadjetpt3nJet;
      TH1D *tagfail1leadjetpt3nJet;
      TH1D *tagfail2leadjetpt3nJet;
      TH1D *tagfail3leadjetpt3nJet;
      TH1D *tagfail4leadjetpt3nJet; 
      TH1D *tagfail5leadjetpt3nJet; 
      TH1D *tagfail6leadjetpt3nJet; 
      TH1D *tagfail7leadjetpt3nJet; 
      TH1D *tagfail8leadjetpt3nJet; 

      TH1D *tagpass0leadjetpt4nJet;
      TH1D *tagpass1leadjetpt4nJet;
      TH1D *tagpass2leadjetpt4nJet;
      TH1D *tagpass3leadjetpt4nJet;
      TH1D *tagpass4leadjetpt4nJet;
      TH1D *tagpass5leadjetpt4nJet;
      TH1D *tagpass6leadjetpt4nJet;
      TH1D *tagpass7leadjetpt4nJet;
      TH1D *tagpass8leadjetpt4nJet;
      TH1D *tagfail0leadjetpt4nJet;
      TH1D *tagfail1leadjetpt4nJet;
      TH1D *tagfail2leadjetpt4nJet;
      TH1D *tagfail3leadjetpt4nJet;
      TH1D *tagfail4leadjetpt4nJet; 
      TH1D *tagfail5leadjetpt4nJet; 
      TH1D *tagfail6leadjetpt4nJet; 
      TH1D *tagfail7leadjetpt4nJet; 
      TH1D *tagfail8leadjetpt4nJet; 

      TH1D *tagpass0leadjetpt5nJet;
      TH1D *tagpass1leadjetpt5nJet;
      TH1D *tagpass2leadjetpt5nJet;
      TH1D *tagpass3leadjetpt5nJet;
      TH1D *tagpass4leadjetpt5nJet;
      TH1D *tagpass5leadjetpt5nJet;
      TH1D *tagpass6leadjetpt5nJet;
      TH1D *tagpass7leadjetpt5nJet;
      TH1D *tagpass8leadjetpt5nJet;
      TH1D *tagfail0leadjetpt5nJet;
      TH1D *tagfail1leadjetpt5nJet;
      TH1D *tagfail2leadjetpt5nJet;
      TH1D *tagfail3leadjetpt5nJet;
      TH1D *tagfail4leadjetpt5nJet; 
      TH1D *tagfail5leadjetpt5nJet; 
      TH1D *tagfail6leadjetpt5nJet; 
      TH1D *tagfail7leadjetpt5nJet; 
      TH1D *tagfail8leadjetpt5nJet; 

      TH1D *HLTnumberOfMatches_PASS;
      TH1D *HLTnumberOfMatches_FAIL;
      TH1D *HLTnumberOfMatches_TOTALELE;
};

//
// constructors and destructor
//
EfficiencyFilter::EfficiencyFilter (const edm::ParameterSet & parameters)
{
  theElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("electronCollection");
  superClusterCollection_EB_ = parameters.getParameter < edm::InputTag > ("superClusterCollection_EB");
  superClusterCollection_EE_ = parameters.getParameter < edm::InputTag > ("superClusterCollection_EE");
  VertexCollectionTag_  = parameters.getParameter<edm::InputTag>("VertexCollectionTag");
  std::string outputfile_D = parameters.getUntrackedParameter<std::string>("filename");
  triggerCollection_=parameters.getUntrackedParameter<edm::InputTag>("triggerCollectionTag");
  useCombinedPrescales_ = parameters.getParameter<bool>("UseCombinedPrescales");
  triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");
  useAllTriggers_       = (triggerNames_.size()==0);
  removePU_             = parameters.getParameter<bool>("removePU");
  WP80_efficiency_                  = parameters.getParameter<bool>("WP80_efficiency");
  HLTele17_efficiency_              = parameters.getParameter<bool>("HLTele17_efficiency");
  HLTele8_efficiency_       = parameters.getParameter<bool>("HLTele8_efficiency");
  RECO_efficiency_                  = parameters.getParameter<bool>("RECO_efficiency");
  New_HE_                           = parameters.getParameter<bool>("New_HE");
  electronIsolatedProducer_ = parameters.getParameter< edm::InputTag > ("electronIsolatedProducer");
  candTag_ = parameters.getParameter< edm::InputTag > ("candTag");
  theJetCollectionLabel_       = parameters.getParameter<edm::InputTag>("JetCollectionLabel");
  triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");


 
  //Initializations...
  edm::Service<TFileService> fs;

  int nInvMassBins = 40;
  double InvMassLowLimit = 71.0;
  double InvMassHighLimit = 111.0;

  probeall_mee  = fs->make<TH1D>("probeall_mee","Invariant mass when probe fails or passes", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail_mee = fs->make<TH1D>("probefail_mee","Invariant mass when probe fails", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass_mee = fs->make<TH1D>("probepass_mee","Invariant mass when probe passes", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0jet = fs->make<TH1D>("probepass0Jet","Invariant mass when probe passes + no Jet", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1jet = fs->make<TH1D>("probepass1Jet","Invariant mass when probe passes + 1 Jet", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2jet = fs->make<TH1D>("probepass2Jet","Invariant mass when probe passes + 2 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3jet = fs->make<TH1D>("probepass3Jet","Invariant mass when probe passes + 3 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4jet = fs->make<TH1D>("probepass4Jet","Invariant mass when probe passes + 4 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5jet = fs->make<TH1D>("probepass5Jet","Invariant mass when probe passes + 5 or more Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0jet = fs->make<TH1D>("probefail0Jet","Invariant mass when probe fail + no Jet", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1jet = fs->make<TH1D>("probefail1Jet","Invariant mass when probe fail + 1 Jet", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2jet = fs->make<TH1D>("probefail2Jet","Invariant mass when probe fail + 2 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3jet = fs->make<TH1D>("probefail3Jet","Invariant mass when probe fail + 3 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4jet = fs->make<TH1D>("probefail4Jet","Invariant mass when probe fail + 4 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5jet = fs->make<TH1D>("probefail5Jet","Invariant mass when probe fail + 5 or more Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0pu = fs->make<TH1D>("probepass0pu","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1pu = fs->make<TH1D>("probepass1pu","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2pu = fs->make<TH1D>("probepass2pu","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3pu = fs->make<TH1D>("probepass3pu","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4pu = fs->make<TH1D>("probepass4pu","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0pu = fs->make<TH1D>("probefail0pu","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1pu = fs->make<TH1D>("probefail1pu","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2pu = fs->make<TH1D>("probefail2pu","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3pu = fs->make<TH1D>("probefail3pu","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4pu = fs->make<TH1D>("probefail4pu","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0eta = fs->make<TH1D>("probepass0eta","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1eta = fs->make<TH1D>("probepass1eta","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2eta = fs->make<TH1D>("probepass2eta","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3eta = fs->make<TH1D>("probepass3eta","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4eta = fs->make<TH1D>("probepass4eta","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0eta = fs->make<TH1D>("probefail0eta","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1eta = fs->make<TH1D>("probefail1eta","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2eta = fs->make<TH1D>("probefail2eta","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3eta = fs->make<TH1D>("probefail3eta","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4eta = fs->make<TH1D>("probefail4eta","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0leadjetpt = fs->make<TH1D>("probepass0leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1leadjetpt = fs->make<TH1D>("probepass1leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2leadjetpt = fs->make<TH1D>("probepass2leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3leadjetpt = fs->make<TH1D>("probepass3leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4leadjetpt = fs->make<TH1D>("probepass4leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5leadjetpt = fs->make<TH1D>("probepass5leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass6leadjetpt = fs->make<TH1D>("probepass6leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass7leadjetpt = fs->make<TH1D>("probepass7leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass8leadjetpt = fs->make<TH1D>("probepass8leadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0leadjetpt = fs->make<TH1D>("probefail0leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1leadjetpt = fs->make<TH1D>("probefail1leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2leadjetpt = fs->make<TH1D>("probefail2leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3leadjetpt = fs->make<TH1D>("probefail3leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4leadjetpt = fs->make<TH1D>("probefail4leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5leadjetpt = fs->make<TH1D>("probefail5leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail6leadjetpt = fs->make<TH1D>("probefail6leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail7leadjetpt = fs->make<TH1D>("probefail7leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail8leadjetpt = fs->make<TH1D>("probefail8leadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0leadjetpt0nJet = fs->make<TH1D>("probepass0leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1leadjetpt0nJet = fs->make<TH1D>("probepass1leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2leadjetpt0nJet = fs->make<TH1D>("probepass2leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3leadjetpt0nJet = fs->make<TH1D>("probepass3leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4leadjetpt0nJet = fs->make<TH1D>("probepass4leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5leadjetpt0nJet = fs->make<TH1D>("probepass5leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass6leadjetpt0nJet = fs->make<TH1D>("probepass6leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass7leadjetpt0nJet = fs->make<TH1D>("probepass7leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass8leadjetpt0nJet = fs->make<TH1D>("probepass8leadjetpt0nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0leadjetpt0nJet = fs->make<TH1D>("probefail0leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1leadjetpt0nJet = fs->make<TH1D>("probefail1leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2leadjetpt0nJet = fs->make<TH1D>("probefail2leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3leadjetpt0nJet = fs->make<TH1D>("probefail3leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4leadjetpt0nJet = fs->make<TH1D>("probefail4leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5leadjetpt0nJet = fs->make<TH1D>("probefail5leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail6leadjetpt0nJet = fs->make<TH1D>("probefail6leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail7leadjetpt0nJet = fs->make<TH1D>("probefail7leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail8leadjetpt0nJet = fs->make<TH1D>("probefail8leadjetpt0nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0leadjetpt1nJet = fs->make<TH1D>("probepass0leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1leadjetpt1nJet = fs->make<TH1D>("probepass1leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2leadjetpt1nJet = fs->make<TH1D>("probepass2leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3leadjetpt1nJet = fs->make<TH1D>("probepass3leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4leadjetpt1nJet = fs->make<TH1D>("probepass4leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5leadjetpt1nJet = fs->make<TH1D>("probepass5leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass6leadjetpt1nJet = fs->make<TH1D>("probepass6leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass7leadjetpt1nJet = fs->make<TH1D>("probepass7leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass8leadjetpt1nJet = fs->make<TH1D>("probepass8leadjetpt1nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0leadjetpt1nJet = fs->make<TH1D>("probefail0leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1leadjetpt1nJet = fs->make<TH1D>("probefail1leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2leadjetpt1nJet = fs->make<TH1D>("probefail2leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3leadjetpt1nJet = fs->make<TH1D>("probefail3leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4leadjetpt1nJet = fs->make<TH1D>("probefail4leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5leadjetpt1nJet = fs->make<TH1D>("probefail5leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail6leadjetpt1nJet = fs->make<TH1D>("probefail6leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail7leadjetpt1nJet = fs->make<TH1D>("probefail7leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail8leadjetpt1nJet = fs->make<TH1D>("probefail8leadjetpt1nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0leadjetpt2nJet = fs->make<TH1D>("probepass0leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1leadjetpt2nJet = fs->make<TH1D>("probepass1leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2leadjetpt2nJet = fs->make<TH1D>("probepass2leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3leadjetpt2nJet = fs->make<TH1D>("probepass3leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4leadjetpt2nJet = fs->make<TH1D>("probepass4leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5leadjetpt2nJet = fs->make<TH1D>("probepass5leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass6leadjetpt2nJet = fs->make<TH1D>("probepass6leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass7leadjetpt2nJet = fs->make<TH1D>("probepass7leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass8leadjetpt2nJet = fs->make<TH1D>("probepass8leadjetpt2nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0leadjetpt2nJet = fs->make<TH1D>("probefail0leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1leadjetpt2nJet = fs->make<TH1D>("probefail1leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2leadjetpt2nJet = fs->make<TH1D>("probefail2leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3leadjetpt2nJet = fs->make<TH1D>("probefail3leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4leadjetpt2nJet = fs->make<TH1D>("probefail4leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5leadjetpt2nJet = fs->make<TH1D>("probefail5leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail6leadjetpt2nJet = fs->make<TH1D>("probefail6leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail7leadjetpt2nJet = fs->make<TH1D>("probefail7leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail8leadjetpt2nJet = fs->make<TH1D>("probefail8leadjetpt2nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0leadjetpt3nJet = fs->make<TH1D>("probepass0leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1leadjetpt3nJet = fs->make<TH1D>("probepass1leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2leadjetpt3nJet = fs->make<TH1D>("probepass2leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3leadjetpt3nJet = fs->make<TH1D>("probepass3leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4leadjetpt3nJet = fs->make<TH1D>("probepass4leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5leadjetpt3nJet = fs->make<TH1D>("probepass5leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass6leadjetpt3nJet = fs->make<TH1D>("probepass6leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass7leadjetpt3nJet = fs->make<TH1D>("probepass7leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass8leadjetpt3nJet = fs->make<TH1D>("probepass8leadjetpt3nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0leadjetpt3nJet = fs->make<TH1D>("probefail0leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1leadjetpt3nJet = fs->make<TH1D>("probefail1leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2leadjetpt3nJet = fs->make<TH1D>("probefail2leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3leadjetpt3nJet = fs->make<TH1D>("probefail3leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4leadjetpt3nJet = fs->make<TH1D>("probefail4leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5leadjetpt3nJet = fs->make<TH1D>("probefail5leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail6leadjetpt3nJet = fs->make<TH1D>("probefail6leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail7leadjetpt3nJet = fs->make<TH1D>("probefail7leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail8leadjetpt3nJet = fs->make<TH1D>("probefail8leadjetpt3nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0leadjetpt4nJet = fs->make<TH1D>("probepass0leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1leadjetpt4nJet = fs->make<TH1D>("probepass1leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2leadjetpt4nJet = fs->make<TH1D>("probepass2leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3leadjetpt4nJet = fs->make<TH1D>("probepass3leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4leadjetpt4nJet = fs->make<TH1D>("probepass4leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5leadjetpt4nJet = fs->make<TH1D>("probepass5leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass6leadjetpt4nJet = fs->make<TH1D>("probepass6leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass7leadjetpt4nJet = fs->make<TH1D>("probepass7leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass8leadjetpt4nJet = fs->make<TH1D>("probepass8leadjetpt4nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0leadjetpt4nJet = fs->make<TH1D>("probefail0leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1leadjetpt4nJet = fs->make<TH1D>("probefail1leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2leadjetpt4nJet = fs->make<TH1D>("probefail2leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3leadjetpt4nJet = fs->make<TH1D>("probefail3leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4leadjetpt4nJet = fs->make<TH1D>("probefail4leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5leadjetpt4nJet = fs->make<TH1D>("probefail5leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail6leadjetpt4nJet = fs->make<TH1D>("probefail6leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail7leadjetpt4nJet = fs->make<TH1D>("probefail7leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail8leadjetpt4nJet = fs->make<TH1D>("probefail8leadjetpt4nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0leadjetpt5nJet = fs->make<TH1D>("probepass0leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1leadjetpt5nJet = fs->make<TH1D>("probepass1leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2leadjetpt5nJet = fs->make<TH1D>("probepass2leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3leadjetpt5nJet = fs->make<TH1D>("probepass3leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4leadjetpt5nJet = fs->make<TH1D>("probepass4leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5leadjetpt5nJet = fs->make<TH1D>("probepass5leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass6leadjetpt5nJet = fs->make<TH1D>("probepass6leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass7leadjetpt5nJet = fs->make<TH1D>("probepass7leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass8leadjetpt5nJet = fs->make<TH1D>("probepass8leadjetpt5nJet","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0leadjetpt5nJet = fs->make<TH1D>("probefail0leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1leadjetpt5nJet = fs->make<TH1D>("probefail1leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2leadjetpt5nJet = fs->make<TH1D>("probefail2leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3leadjetpt5nJet = fs->make<TH1D>("probefail3leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4leadjetpt5nJet = fs->make<TH1D>("probefail4leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5leadjetpt5nJet = fs->make<TH1D>("probefail5leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail6leadjetpt5nJet = fs->make<TH1D>("probefail6leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail7leadjetpt5nJet = fs->make<TH1D>("probefail7leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail8leadjetpt5nJet = fs->make<TH1D>("probefail8leadjetpt5nJet","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagall_mee  = fs->make<TH1D>("tagall_mee","Invariant mass when tag fails or passes", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail_mee = fs->make<TH1D>("tagfail_mee","Invariant mass when tag fails WP80", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass_mee = fs->make<TH1D>("tagpass_mee","Invariant mass when tag passes WP80", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0jet = fs->make<TH1D>("tagpass0Jet","Invariant mass when tag passes + no Jet", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1jet = fs->make<TH1D>("tagpass1Jet","Invariant mass when tag passes + 1 Jet", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2jet = fs->make<TH1D>("tagpass2Jet","Invariant mass when tag passes + 2 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3jet = fs->make<TH1D>("tagpass3Jet","Invariant mass when tag passes + 3 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4jet = fs->make<TH1D>("tagpass4Jet","Invariant mass when tag passes + 4 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5jet = fs->make<TH1D>("tagpass5Jet","Invariant mass when tag passes + 5 or more Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0jet = fs->make<TH1D>("tagfail0Jet","Invariant mass when tag fail + no Jet", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1jet = fs->make<TH1D>("tagfail1Jet","Invariant mass when tag fail + 1 Jet", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2jet = fs->make<TH1D>("tagfail2Jet","Invariant mass when tag fail + 2 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3jet = fs->make<TH1D>("tagfail3Jet","Invariant mass when tag fail + 3 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4jet = fs->make<TH1D>("tagfail4Jet","Invariant mass when tag fail + 4 Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5jet = fs->make<TH1D>("tagfail5Jet","Invariant mass when tag fail + 5 or more Jets", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0pu = fs->make<TH1D>("tagpass0pu","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1pu = fs->make<TH1D>("tagpass1pu","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2pu = fs->make<TH1D>("tagpass2pu","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3pu = fs->make<TH1D>("tagpass3pu","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4pu = fs->make<TH1D>("tagpass4pu","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0pu = fs->make<TH1D>("tagfail0pu","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1pu = fs->make<TH1D>("tagfail1pu","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2pu = fs->make<TH1D>("tagfail2pu","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3pu = fs->make<TH1D>("tagfail3pu","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4pu = fs->make<TH1D>("tagfail4pu","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0eta = fs->make<TH1D>("tagpass0eta","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1eta = fs->make<TH1D>("tagpass1eta","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2eta = fs->make<TH1D>("tagpass2eta","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3eta = fs->make<TH1D>("tagpass3eta","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4eta = fs->make<TH1D>("tagpass4eta","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0eta = fs->make<TH1D>("tagfail0eta","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1eta = fs->make<TH1D>("tagfail1eta","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2eta = fs->make<TH1D>("tagfail2eta","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3eta = fs->make<TH1D>("tagfail3eta","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4eta = fs->make<TH1D>("tagfail4eta","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0leadjetpt = fs->make<TH1D>("tagpass0leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1leadjetpt = fs->make<TH1D>("tagpass1leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2leadjetpt = fs->make<TH1D>("tagpass2leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3leadjetpt = fs->make<TH1D>("tagpass3leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4leadjetpt = fs->make<TH1D>("tagpass4leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5leadjetpt = fs->make<TH1D>("tagpass5leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass6leadjetpt = fs->make<TH1D>("tagpass6leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass7leadjetpt = fs->make<TH1D>("tagpass7leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass8leadjetpt = fs->make<TH1D>("tagpass8leadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0leadjetpt = fs->make<TH1D>("tagfail0leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1leadjetpt = fs->make<TH1D>("tagfail1leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2leadjetpt = fs->make<TH1D>("tagfail2leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3leadjetpt = fs->make<TH1D>("tagfail3leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4leadjetpt = fs->make<TH1D>("tagfail4leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5leadjetpt = fs->make<TH1D>("tagfail5leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail6leadjetpt = fs->make<TH1D>("tagfail6leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail7leadjetpt = fs->make<TH1D>("tagfail7leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail8leadjetpt = fs->make<TH1D>("tagfail8leadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0leadjetpt0nJet = fs->make<TH1D>("tagpass0leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1leadjetpt0nJet = fs->make<TH1D>("tagpass1leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2leadjetpt0nJet = fs->make<TH1D>("tagpass2leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3leadjetpt0nJet = fs->make<TH1D>("tagpass3leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4leadjetpt0nJet = fs->make<TH1D>("tagpass4leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5leadjetpt0nJet = fs->make<TH1D>("tagpass5leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass6leadjetpt0nJet = fs->make<TH1D>("tagpass6leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass7leadjetpt0nJet = fs->make<TH1D>("tagpass7leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass8leadjetpt0nJet = fs->make<TH1D>("tagpass8leadjetpt0nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0leadjetpt0nJet = fs->make<TH1D>("tagfail0leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1leadjetpt0nJet = fs->make<TH1D>("tagfail1leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2leadjetpt0nJet = fs->make<TH1D>("tagfail2leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3leadjetpt0nJet = fs->make<TH1D>("tagfail3leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4leadjetpt0nJet = fs->make<TH1D>("tagfail4leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5leadjetpt0nJet = fs->make<TH1D>("tagfail5leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail6leadjetpt0nJet = fs->make<TH1D>("tagfail6leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail7leadjetpt0nJet = fs->make<TH1D>("tagfail7leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail8leadjetpt0nJet = fs->make<TH1D>("tagfail8leadjetpt0nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0leadjetpt1nJet = fs->make<TH1D>("tagpass0leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1leadjetpt1nJet = fs->make<TH1D>("tagpass1leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2leadjetpt1nJet = fs->make<TH1D>("tagpass2leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3leadjetpt1nJet = fs->make<TH1D>("tagpass3leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4leadjetpt1nJet = fs->make<TH1D>("tagpass4leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5leadjetpt1nJet = fs->make<TH1D>("tagpass5leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass6leadjetpt1nJet = fs->make<TH1D>("tagpass6leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass7leadjetpt1nJet = fs->make<TH1D>("tagpass7leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass8leadjetpt1nJet = fs->make<TH1D>("tagpass8leadjetpt1nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0leadjetpt1nJet = fs->make<TH1D>("tagfail0leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1leadjetpt1nJet = fs->make<TH1D>("tagfail1leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2leadjetpt1nJet = fs->make<TH1D>("tagfail2leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3leadjetpt1nJet = fs->make<TH1D>("tagfail3leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4leadjetpt1nJet = fs->make<TH1D>("tagfail4leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5leadjetpt1nJet = fs->make<TH1D>("tagfail5leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail6leadjetpt1nJet = fs->make<TH1D>("tagfail6leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail7leadjetpt1nJet = fs->make<TH1D>("tagfail7leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail8leadjetpt1nJet = fs->make<TH1D>("tagfail8leadjetpt1nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0leadjetpt2nJet = fs->make<TH1D>("tagpass0leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1leadjetpt2nJet = fs->make<TH1D>("tagpass1leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2leadjetpt2nJet = fs->make<TH1D>("tagpass2leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3leadjetpt2nJet = fs->make<TH1D>("tagpass3leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4leadjetpt2nJet = fs->make<TH1D>("tagpass4leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5leadjetpt2nJet = fs->make<TH1D>("tagpass5leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass6leadjetpt2nJet = fs->make<TH1D>("tagpass6leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass7leadjetpt2nJet = fs->make<TH1D>("tagpass7leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass8leadjetpt2nJet = fs->make<TH1D>("tagpass8leadjetpt2nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0leadjetpt2nJet = fs->make<TH1D>("tagfail0leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1leadjetpt2nJet = fs->make<TH1D>("tagfail1leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2leadjetpt2nJet = fs->make<TH1D>("tagfail2leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3leadjetpt2nJet = fs->make<TH1D>("tagfail3leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4leadjetpt2nJet = fs->make<TH1D>("tagfail4leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5leadjetpt2nJet = fs->make<TH1D>("tagfail5leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail6leadjetpt2nJet = fs->make<TH1D>("tagfail6leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail7leadjetpt2nJet = fs->make<TH1D>("tagfail7leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail8leadjetpt2nJet = fs->make<TH1D>("tagfail8leadjetpt2nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0leadjetpt3nJet = fs->make<TH1D>("tagpass0leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1leadjetpt3nJet = fs->make<TH1D>("tagpass1leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2leadjetpt3nJet = fs->make<TH1D>("tagpass2leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3leadjetpt3nJet = fs->make<TH1D>("tagpass3leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4leadjetpt3nJet = fs->make<TH1D>("tagpass4leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5leadjetpt3nJet = fs->make<TH1D>("tagpass5leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass6leadjetpt3nJet = fs->make<TH1D>("tagpass6leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass7leadjetpt3nJet = fs->make<TH1D>("tagpass7leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass8leadjetpt3nJet = fs->make<TH1D>("tagpass8leadjetpt3nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0leadjetpt3nJet = fs->make<TH1D>("tagfail0leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1leadjetpt3nJet = fs->make<TH1D>("tagfail1leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2leadjetpt3nJet = fs->make<TH1D>("tagfail2leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3leadjetpt3nJet = fs->make<TH1D>("tagfail3leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4leadjetpt3nJet = fs->make<TH1D>("tagfail4leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5leadjetpt3nJet = fs->make<TH1D>("tagfail5leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail6leadjetpt3nJet = fs->make<TH1D>("tagfail6leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail7leadjetpt3nJet = fs->make<TH1D>("tagfail7leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail8leadjetpt3nJet = fs->make<TH1D>("tagfail8leadjetpt3nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0leadjetpt4nJet = fs->make<TH1D>("tagpass0leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1leadjetpt4nJet = fs->make<TH1D>("tagpass1leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2leadjetpt4nJet = fs->make<TH1D>("tagpass2leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3leadjetpt4nJet = fs->make<TH1D>("tagpass3leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4leadjetpt4nJet = fs->make<TH1D>("tagpass4leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5leadjetpt4nJet = fs->make<TH1D>("tagpass5leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass6leadjetpt4nJet = fs->make<TH1D>("tagpass6leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass7leadjetpt4nJet = fs->make<TH1D>("tagpass7leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass8leadjetpt4nJet = fs->make<TH1D>("tagpass8leadjetpt4nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0leadjetpt4nJet = fs->make<TH1D>("tagfail0leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1leadjetpt4nJet = fs->make<TH1D>("tagfail1leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2leadjetpt4nJet = fs->make<TH1D>("tagfail2leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3leadjetpt4nJet = fs->make<TH1D>("tagfail3leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4leadjetpt4nJet = fs->make<TH1D>("tagfail4leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5leadjetpt4nJet = fs->make<TH1D>("tagfail5leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail6leadjetpt4nJet = fs->make<TH1D>("tagfail6leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail7leadjetpt4nJet = fs->make<TH1D>("tagfail7leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail8leadjetpt4nJet = fs->make<TH1D>("tagfail8leadjetpt4nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0leadjetpt5nJet = fs->make<TH1D>("tagpass0leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1leadjetpt5nJet = fs->make<TH1D>("tagpass1leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2leadjetpt5nJet = fs->make<TH1D>("tagpass2leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3leadjetpt5nJet = fs->make<TH1D>("tagpass3leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4leadjetpt5nJet = fs->make<TH1D>("tagpass4leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5leadjetpt5nJet = fs->make<TH1D>("tagpass5leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass6leadjetpt5nJet = fs->make<TH1D>("tagpass6leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass7leadjetpt5nJet = fs->make<TH1D>("tagpass7leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass8leadjetpt5nJet = fs->make<TH1D>("tagpass8leadjetpt5nJet","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0leadjetpt5nJet = fs->make<TH1D>("tagfail0leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1leadjetpt5nJet = fs->make<TH1D>("tagfail1leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2leadjetpt5nJet = fs->make<TH1D>("tagfail2leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3leadjetpt5nJet = fs->make<TH1D>("tagfail3leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4leadjetpt5nJet = fs->make<TH1D>("tagfail4leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5leadjetpt5nJet = fs->make<TH1D>("tagfail5leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail6leadjetpt5nJet = fs->make<TH1D>("tagfail6leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail7leadjetpt5nJet = fs->make<TH1D>("tagfail7leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail8leadjetpt5nJet = fs->make<TH1D>("tagfail8leadjetpt5nJet","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probeall_leadjetpt= fs->make<TH1D>("probeall_leadjetpt","Pt of the leading jet", 200, 0, 200.0);
  tagall_leadjetpt= fs->make<TH1D>("tagall_leadjetpt","Pt of the leading jet", 200, 0, 200.0);

  probeall_pt= fs->make<TH1D>("probeall_pt","Pt of the electron probe", 120, 0, 120.0);
  probepass_pt= fs->make<TH1D>("probepass_pt","Pt of the electron probe when passing WP80", 120, 0, 120.0);
  probefail_pt= fs->make<TH1D>("probefail_pt","Pt of the electron probe when failing WP80", 120, 0, 120.0);
  tagall_pt= fs->make<TH1D>("tagall_pt","Pt of the electron tag", 120, 0, 120.0);
  tagpass_pt= fs->make<TH1D>("tagpass_pt","Pt of the electron tag when passing WP80", 120, 0, 120.0);
  tagfail_pt= fs->make<TH1D>("tagfail_pt","Pt of the electron tag when failing WP80", 120, 0, 120.0);

  probeall_eta= fs->make<TH1D>("probeall_eta","Eta of the electron probe", 60, -3.0, 3.0);
  probepass_eta= fs->make<TH1D>("probepass_eta","Eta of the electron probe when passing WP80", 60, -3.0, 3.0);
  probefail_eta= fs->make<TH1D>("probefail_eta","Eta of the electron probe when failing WP80", 60, -3.0, 3.0);
  tagall_eta= fs->make<TH1D>("tagall_eta","Eta of the electron tag", 60, -3.0, 3.0);
  tagpass_eta= fs->make<TH1D>("tagpass_eta","Eta of the electron tag when passing WP80", 60, -3.0, 3.0);
  tagfail_eta= fs->make<TH1D>("tagfail_eta","Eta of the electron tag when failing WP80", 60, -3.0, 3.0);

  //  scNumber_per_event    = fs->make<TH1D>("scNumber_per_event","Total # of SC > 5 GeV",10,0,10);
  //  eleNumber_scMatch    = fs->make<TH1D>("eleNumber_scMatch","Total # of SC > 5 GeV",10,0,10);

  HLTnumberOfMatches_PASS    = fs->make<TH1D>("HLTnumberOfMatches_PASS","Total # of electron matching the HLT",10,0,10);
  HLTnumberOfMatches_FAIL    = fs->make<TH1D>("HLTnumberOfMatches_FAIL","Total # of electron NOT matching the HLT",10,0,10);
  HLTnumberOfMatches_TOTALELE= fs->make<TH1D>("HLTnumberOfMatches_TOTALELE","Total # of electrons per event",10,0,10);
}



EfficiencyFilter::~EfficiencyFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}



