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
#include "VecBosSelection/Selection/interface/SelectionUtils.h"

class EfficiencyPtEtaFilter : public edm::EDFilter, public SelectionUtils {
   public:
      explicit EfficiencyPtEtaFilter(const edm::ParameterSet &);
      ~EfficiencyPtEtaFilter();

       virtual void beginJob();
       virtual bool beginRun(edm::Run &, edm::EventSetup const&);

   private:
      virtual bool filter(edm::Event&, edm::EventSetup const&);
      virtual void endJob() ;

      // ----------member data ---------------------------

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

      TH1D *probepass1eta1pt;
      TH1D *probepass1eta2pt;
      TH1D *probepass1eta3pt;
      TH1D *probepass1eta4pt;
      TH1D *probepass2eta1pt;
      TH1D *probepass2eta2pt;
      TH1D *probepass2eta3pt;
      TH1D *probepass2eta4pt;
      TH1D *probepass3eta1pt;
      TH1D *probepass3eta2pt;
      TH1D *probepass3eta3pt;
      TH1D *probepass3eta4pt;
      TH1D *probepass4eta1pt;
      TH1D *probepass4eta2pt;
      TH1D *probepass4eta3pt;
      TH1D *probepass4eta4pt;

      TH1D *probefail1eta1pt;
      TH1D *probefail1eta2pt;
      TH1D *probefail1eta3pt;
      TH1D *probefail1eta4pt;
      TH1D *probefail2eta1pt;
      TH1D *probefail2eta2pt;
      TH1D *probefail2eta3pt;
      TH1D *probefail2eta4pt;
      TH1D *probefail3eta1pt;
      TH1D *probefail3eta2pt;
      TH1D *probefail3eta3pt;
      TH1D *probefail3eta4pt;
      TH1D *probefail4eta1pt;
      TH1D *probefail4eta2pt;
      TH1D *probefail4eta3pt;
      TH1D *probefail4eta4pt;

      TH1D *probeall_mee;
      TH1D *probepass_mee;
      TH1D *probefail_mee;

      TH1D *probeall_pt;
      TH1D *probepass_pt;
      TH1D *probefail_pt;

      TH1D *probeall_eta;
      TH1D *probepass_eta;
      TH1D *probefail_eta;

      TH1D *tagall_pt;
      TH1D *tagall_eta;

};

//
// constructors and destructor
//
EfficiencyPtEtaFilter::EfficiencyPtEtaFilter (const edm::ParameterSet & parameters)
{
  theElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("electronCollection");
  theTagHLTElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("TagHLTelectronCollection");
  theProbeHLTElectronCollectionLabel = parameters.getParameter < edm::InputTag > ("ProbeHLTelectronCollection");
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
  electronIsolatedProducer_ = parameters.getParameter< edm::InputTag > ("electronIsolatedProducer");
  candTag_ = parameters.getParameter< edm::InputTag > ("candTag");
  theJetCollectionLabel_       = parameters.getParameter<edm::InputTag>("JetCollectionLabel");
  triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");


 
  //Initializations...
  edm::Service<TFileService> fs;

  int nInvMassBins = 60;
  double InvMassLowLimit = 60.0;
  double InvMassHighLimit = 120.0;

  probepass1eta1pt = fs->make<TH1D>("probepass1eta1pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1eta2pt = fs->make<TH1D>("probepass1eta2pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1eta3pt = fs->make<TH1D>("probepass1eta3pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass1eta4pt = fs->make<TH1D>("probepass1eta4pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2eta1pt = fs->make<TH1D>("probepass2eta1pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2eta2pt = fs->make<TH1D>("probepass2eta2pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2eta3pt = fs->make<TH1D>("probepass2eta3pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass2eta4pt = fs->make<TH1D>("probepass2eta4pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3eta1pt = fs->make<TH1D>("probepass3eta1pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3eta2pt = fs->make<TH1D>("probepass3eta2pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3eta3pt = fs->make<TH1D>("probepass3eta3pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass3eta4pt = fs->make<TH1D>("probepass3eta4pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4eta1pt = fs->make<TH1D>("probepass4eta1pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4eta2pt = fs->make<TH1D>("probepass4eta2pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4eta3pt = fs->make<TH1D>("probepass4eta3pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass4eta4pt = fs->make<TH1D>("probepass4eta4pt","Invariant mass when probe pass", nInvMassBins, InvMassLowLimit, InvMassHighLimit);

  probefail1eta1pt = fs->make<TH1D>("probefail1eta1pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1eta2pt = fs->make<TH1D>("probefail1eta2pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1eta3pt = fs->make<TH1D>("probefail1eta3pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail1eta4pt = fs->make<TH1D>("probefail1eta4pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2eta1pt = fs->make<TH1D>("probefail2eta1pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2eta2pt = fs->make<TH1D>("probefail2eta2pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2eta3pt = fs->make<TH1D>("probefail2eta3pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail2eta4pt = fs->make<TH1D>("probefail2eta4pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3eta1pt = fs->make<TH1D>("probefail3eta1pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3eta2pt = fs->make<TH1D>("probefail3eta2pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3eta3pt = fs->make<TH1D>("probefail3eta3pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail3eta4pt = fs->make<TH1D>("probefail3eta4pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4eta1pt = fs->make<TH1D>("probefail4eta1pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4eta2pt = fs->make<TH1D>("probefail4eta2pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4eta3pt = fs->make<TH1D>("probefail4eta3pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail4eta4pt = fs->make<TH1D>("probefail4eta4pt","Invariant mass when probe fail", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
 
  probeall_mee  = fs->make<TH1D>("probeall_mee","Invariant mass when probe fails or passes", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probefail_mee = fs->make<TH1D>("probefail_mee","Invariant mass when probe fails", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
  probepass_mee = fs->make<TH1D>("probepass_mee","Invariant mass when probe passes", nInvMassBins, InvMassLowLimit, InvMassHighLimit);
 
  probeall_pt= fs->make<TH1D>("probeall_pt","Pt of the electron probe", 120, 0, 120.0);
  probepass_pt= fs->make<TH1D>("probepass_pt","Pt of the electron probe when passing WP80", 120, 0, 120.0);
  probefail_pt= fs->make<TH1D>("probefail_pt","Pt of the electron probe when failing WP80", 120, 0, 120.0);

  probeall_eta= fs->make<TH1D>("probeall_eta","Eta of the electron probe", 60, -3.0, 3.0);
  probepass_eta= fs->make<TH1D>("probepass_eta","Eta of the electron probe when passing WP80", 60, -3.0, 3.0);
  probefail_eta= fs->make<TH1D>("probefail_eta","Eta of the electron probe when failing WP80", 60, -3.0, 3.0);

  tagall_eta= fs->make<TH1D>("tagall_eta","Eta of the electron tag", 60, -3.0, 3.0);
  tagall_pt= fs->make<TH1D>("tagall_pt","Pt of the electron tag", 120, 0, 120.0);
}



EfficiencyPtEtaFilter::~EfficiencyPtEtaFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}



