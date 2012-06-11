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

      bool matchMC_;
      edm::InputTag genParticleCollection_;
      edm::InputTag theElectronCollectionLabel;
      edm::InputTag theTagHLTElectronCollectionLabel;
      edm::InputTag theProbeHLTElectronCollectionLabel;
      edm::InputTag superClusterCollection_EB_;
      edm::InputTag superClusterCollection_EE_;
      edm::InputTag triggerCollection_; 
      edm::InputTag electronIsolatedProducer_;
      edm::InputTag candTag_ ;
      edm::InputTag theJetCollectionLabel_;
      edm::InputTag VertexCollectionTag_;
      std::vector<edm::InputTag>  isoValInputTags_;
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

      TH1D *probepass0subleadjetpt;
      TH1D *probepass1subleadjetpt;
      TH1D *probepass2subleadjetpt;
      TH1D *probepass3subleadjetpt;
      TH1D *probepass4subleadjetpt;
      TH1D *probepass5subleadjetpt;
      TH1D *probefail0subleadjetpt ;
      TH1D *probefail1subleadjetpt ;
      TH1D *probefail2subleadjetpt ;
      TH1D *probefail3subleadjetpt ;
      TH1D *probefail4subleadjetpt; 
      TH1D *probefail5subleadjetpt; 

      TH1D *probepass0subsubleadjetpt;
      TH1D *probepass1subsubleadjetpt;
      TH1D *probepass2subsubleadjetpt;
      TH1D *probefail0subsubleadjetpt ;
      TH1D *probefail1subsubleadjetpt ;
      TH1D *probefail2subsubleadjetpt ;

      TH1D *probepass0subsubsubleadjetpt;
      TH1D *probefail0subsubsubleadjetpt ;

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

      TH1D *tagpass0subleadjetpt;
      TH1D *tagpass1subleadjetpt;
      TH1D *tagpass2subleadjetpt;
      TH1D *tagpass3subleadjetpt;
      TH1D *tagpass4subleadjetpt;
      TH1D *tagpass5subleadjetpt;
      TH1D *tagfail0subleadjetpt ;
      TH1D *tagfail1subleadjetpt ;
      TH1D *tagfail2subleadjetpt ;
      TH1D *tagfail3subleadjetpt ;
      TH1D *tagfail4subleadjetpt; 
      TH1D *tagfail5subleadjetpt; 

      TH1D *tagpass0subsubleadjetpt;
      TH1D *tagpass1subsubleadjetpt;
      TH1D *tagpass2subsubleadjetpt;
      TH1D *tagfail0subsubleadjetpt ;
      TH1D *tagfail1subsubleadjetpt ;
      TH1D *tagfail2subsubleadjetpt ;

      TH1D *tagpass0subsubsubleadjetpt;
      TH1D *tagfail0subsubsubleadjetpt ;
};

//
// constructors and destructor
//
EfficiencyFilter::EfficiencyFilter (const edm::ParameterSet & parameters)
{
  matchMC_ = parameters.getParameter<bool>("matchMC");
  genParticleCollection_ = parameters.getUntrackedParameter<edm::InputTag>("genParticleCollection", edm::InputTag("genParticles"));
  theElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("electronCollection");
  superClusterCollection_EB_ = parameters.getParameter < edm::InputTag > ("superClusterCollection_EB");
  superClusterCollection_EE_ = parameters.getParameter < edm::InputTag > ("superClusterCollection_EE");
  theTagHLTElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("TagHLTelectronCollection");
  theProbeHLTElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("ProbeHLTelectronCollection");
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
  electronIsolatedProducer_ = parameters.getParameter< edm::InputTag > ("electronIsolatedProducer");
  candTag_ = parameters.getParameter< edm::InputTag > ("candTag");
  theJetCollectionLabel_       = parameters.getParameter<edm::InputTag>("JetCollectionLabel");
  triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");
  isoValInputTags_      = parameters.getParameter<std::vector<edm::InputTag> >("isoValInputTags");


 
  //Initializations...
  edm::Service<TFileService> fs;

  int nInvMassBins = 60;
  double InvMassLowLimit = 60.0;
  double InvMassHighLimit = 120.0;

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

  probepass0subleadjetpt = fs->make<TH1D>("probepass0subleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1subleadjetpt = fs->make<TH1D>("probepass1subleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2subleadjetpt = fs->make<TH1D>("probepass2subleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3subleadjetpt = fs->make<TH1D>("probepass3subleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4subleadjetpt = fs->make<TH1D>("probepass4subleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass5subleadjetpt = fs->make<TH1D>("probepass5subleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0subleadjetpt = fs->make<TH1D>("probefail0subleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1subleadjetpt = fs->make<TH1D>("probefail1subleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2subleadjetpt = fs->make<TH1D>("probefail2subleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3subleadjetpt = fs->make<TH1D>("probefail3subleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4subleadjetpt = fs->make<TH1D>("probefail4subleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail5subleadjetpt = fs->make<TH1D>("probefail5subleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0subsubleadjetpt = fs->make<TH1D>("probepass0subsubleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1subsubleadjetpt = fs->make<TH1D>("probepass1subsubleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2subsubleadjetpt = fs->make<TH1D>("probepass2subsubleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail0subsubleadjetpt = fs->make<TH1D>("probefail0subsubleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1subsubleadjetpt = fs->make<TH1D>("probefail1subsubleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2subsubleadjetpt = fs->make<TH1D>("probefail2subsubleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probepass0subsubsubleadjetpt = fs->make<TH1D>("probepass0subsubsubleadjetpt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail0subsubsubleadjetpt = fs->make<TH1D>("probefail0subsubsubleadjetpt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

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

  tagpass0subleadjetpt = fs->make<TH1D>("tagpass0subleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1subleadjetpt = fs->make<TH1D>("tagpass1subleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2subleadjetpt = fs->make<TH1D>("tagpass2subleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass3subleadjetpt = fs->make<TH1D>("tagpass3subleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass4subleadjetpt = fs->make<TH1D>("tagpass4subleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass5subleadjetpt = fs->make<TH1D>("tagpass5subleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0subleadjetpt = fs->make<TH1D>("tagfail0subleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1subleadjetpt = fs->make<TH1D>("tagfail1subleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2subleadjetpt = fs->make<TH1D>("tagfail2subleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail3subleadjetpt = fs->make<TH1D>("tagfail3subleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail4subleadjetpt = fs->make<TH1D>("tagfail4subleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail5subleadjetpt = fs->make<TH1D>("tagfail5subleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0subsubleadjetpt = fs->make<TH1D>("tagpass0subsubleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass1subsubleadjetpt = fs->make<TH1D>("tagpass1subsubleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagpass2subsubleadjetpt = fs->make<TH1D>("tagpass2subsubleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagfail0subsubleadjetpt = fs->make<TH1D>("tagfail0subsubleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail1subsubleadjetpt = fs->make<TH1D>("tagfail1subsubleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail2subsubleadjetpt = fs->make<TH1D>("tagfail2subsubleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  tagpass0subsubsubleadjetpt = fs->make<TH1D>("tagpass0subsubsubleadjetpt","Invariant mass when tag pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  tagfail0subsubsubleadjetpt = fs->make<TH1D>("tagfail0subsubsubleadjetpt","Invariant mass when tag fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

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
}



EfficiencyFilter::~EfficiencyFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}



