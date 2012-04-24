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
      edm::InputTag theHLTElectronCollectionLabel;
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

      TH1D *probepass0th2f;
      TH1D *probepass1th2f;
      TH1D *probepass2th2f;
      TH1D *probepass3th2f;
      TH1D *probepass4th2f;
      TH1D *probepass5th2f;
      TH1D *probepass6th2f;
      TH1D *probepass7th2f;
      TH1D *probepass8th2f;
      TH1D *probepass9th2f;
      TH1D *probepass10th2f;
      TH1D *probepass11th2f;
      TH1D *probepass12th2f;
      TH1D *probefail0th2f;
      TH1D *probefail1th2f;
      TH1D *probefail2th2f;
      TH1D *probefail3th2f;
      TH1D *probefail4th2f; 
      TH1D *probefail5th2f; 
      TH1D *probefail6th2f; 
      TH1D *probefail7th2f; 
      TH1D *probefail8th2f; 
      TH1D *probefail9th2f; 
      TH1D *probefail10th2f; 
      TH1D *probefail11th2f; 
      TH1D *probefail12th2f; 

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

      TH1D *tagpass0th2f;
      TH1D *tagpass1th2f;
      TH1D *tagpass2th2f;
      TH1D *tagpass3th2f;
      TH1D *tagpass4th2f;
      TH1D *tagpass5th2f;
      TH1D *tagpass6th2f;
      TH1D *tagpass7th2f;
      TH1D *tagpass8th2f;
      TH1D *tagpass9th2f;
      TH1D *tagpass10th2f;
      TH1D *tagpass11th2f;
      TH1D *tagpass12th2f;
      TH1D *tagfail0th2f;
      TH1D *tagfail1th2f;
      TH1D *tagfail2th2f;
      TH1D *tagfail3th2f;
      TH1D *tagfail4th2f; 
      TH1D *tagfail5th2f; 
      TH1D *tagfail6th2f; 
      TH1D *tagfail7th2f; 
      TH1D *tagfail8th2f; 
      TH1D *tagfail9th2f; 
      TH1D *tagfail10th2f; 
      TH1D *tagfail11th2f; 
      TH1D *tagfail12th2f; 

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
  theHLTElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("HLTelectronCollection");
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

  probepass0th2f = fs->make<TH1D>("probepass0th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1th2f = fs->make<TH1D>("probepass1th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2th2f = fs->make<TH1D>("probepass2th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3th2f = fs->make<TH1D>("probepass3th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4th2f = fs->make<TH1D>("probepass4th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5th2f = fs->make<TH1D>("probepass5th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass6th2f = fs->make<TH1D>("probepass6th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass7th2f = fs->make<TH1D>("probepass7th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass8th2f = fs->make<TH1D>("probepass8th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass9th2f = fs->make<TH1D>("probepass9th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass10th2f = fs->make<TH1D>("probepass10th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass11th2f = fs->make<TH1D>("probepass11th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass12th2f = fs->make<TH1D>("probepass12th2f","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0th2f = fs->make<TH1D>("probefail0th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1th2f = fs->make<TH1D>("probefail1th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2th2f = fs->make<TH1D>("probefail2th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3th2f = fs->make<TH1D>("probefail3th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4th2f = fs->make<TH1D>("probefail4th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5th2f = fs->make<TH1D>("probefail5th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail6th2f = fs->make<TH1D>("probefail6th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail7th2f = fs->make<TH1D>("probefail7th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail8th2f = fs->make<TH1D>("probefail8th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail9th2f = fs->make<TH1D>("probefail9th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail10th2f = fs->make<TH1D>("probefail10th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail11th2f = fs->make<TH1D>("probefail11th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail12th2f = fs->make<TH1D>("probefail12th2f","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

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

  tagpass0th2f = fs->make<TH1D>("tagpass0th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1th2f = fs->make<TH1D>("tagpass1th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2th2f = fs->make<TH1D>("tagpass2th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3th2f = fs->make<TH1D>("tagpass3th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4th2f = fs->make<TH1D>("tagpass4th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5th2f = fs->make<TH1D>("tagpass5th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass6th2f = fs->make<TH1D>("tagpass6th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass7th2f = fs->make<TH1D>("tagpass7th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass8th2f = fs->make<TH1D>("tagpass8th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass9th2f = fs->make<TH1D>("tagpass9th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass10th2f = fs->make<TH1D>("tagpass10th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass11th2f = fs->make<TH1D>("tagpass11th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass12th2f = fs->make<TH1D>("tagpass12th2f","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0th2f = fs->make<TH1D>("tagfail0th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1th2f = fs->make<TH1D>("tagfail1th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2th2f = fs->make<TH1D>("tagfail2th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3th2f = fs->make<TH1D>("tagfail3th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4th2f = fs->make<TH1D>("tagfail4th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5th2f = fs->make<TH1D>("tagfail5th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail6th2f = fs->make<TH1D>("tagfail6th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail7th2f = fs->make<TH1D>("tagfail7th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail8th2f = fs->make<TH1D>("tagfail8th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail9th2f = fs->make<TH1D>("tagfail9th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail10th2f = fs->make<TH1D>("tagfail10th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail11th2f = fs->make<TH1D>("tagfail11th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail12th2f = fs->make<TH1D>("tagfail12th2f","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

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



