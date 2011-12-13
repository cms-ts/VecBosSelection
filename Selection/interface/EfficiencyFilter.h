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
      TH1D *probefailed0jet ;
      TH1D *probefailed1jet ;
      TH1D *probefailed2jet ;
      TH1D *probefailed3jet ;
      TH1D *probefailed4jet; 
      TH1D *probepass1jetEE;
      TH1D *probepass2jetEE;
      TH1D *probepass3jetEE;
      TH1D *probepass4jetEE;
      TH1D *probepass0jetEE;
      TH1D *probefailed0jetEE;
      TH1D *probefailed1jetEE;
      TH1D *probefailed2jetEE;
      TH1D *probefailed3jetEE;
      TH1D *probefailed4jetEE; 
      TH1D *probepass1jetEB;
      TH1D *probepass2jetEB;
      TH1D *probepass3jetEB;
      TH1D *probepass4jetEB;
      TH1D *probepass0jetEB;
      TH1D *probefailed0jetEB;
      TH1D *probefailed1jetEB;
      TH1D *probefailed2jetEB;
      TH1D *probefailed3jetEB;
      TH1D *probefailed4jetEB; 
      TH1D *probeinEB ;
      TH1D *probeinEE ;
      TH1D *probept;
      TH1D *tagpt;
      TH1D *probept_0_35;
      TH1D *probept_35_45;
      TH1D *probept_45_55;
      TH1D *probept_55_inf;
      //HLT*
      TH1D *HLTEfficiency;
      TH1D *HLTEfficiency0Jet;
      TH1D *HLTEfficiency1Jet;
      TH1D *HLTEfficiency2Jet;
      TH1D *HLTEfficiency3Jet;
      TH1D *HLTEfficiency4Jet;
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

  probepass0jetEE = fs->make<TH1D>("Probepass0JetEE","Invariant mass when probe passes + no Jet and it is contained in EE", 60, 60.0, 120.0);
  probepass1jetEE = fs->make<TH1D>("Probepass1JetEE","Invariant mass when probe passes + 1 Jet and it is contained in EE", 60, 60.0, 120.0);
  probepass2jetEE = fs->make<TH1D>("Probepass2JetEE","Invariant mass when probe passes + 2 Jets and it is contained in EE", 60, 60.0, 120.0);
  probepass3jetEE = fs->make<TH1D>("Probepass3JetEE","Invariant mass when probe passes + 3 Jets and it is contained in EE", 60, 60.0, 120.0);
  probepass4jetEE = fs->make<TH1D>("Probepass4JetEE","Invariant mass when probe passes + 4 Jets and it is contained in EE", 60, 60.0, 120.0);
  probepass0jetEB = fs->make<TH1D>("Probepass0JetEB","Invariant mass when probe passes + no Jet and it is contained in EB", 60, 60.0, 120.0);
  probepass1jetEB = fs->make<TH1D>("Probepass1JetEB","Invariant mass when probe passes + 1 Jet and it is contained in EB", 60, 60.0, 120.0);
  probepass2jetEB = fs->make<TH1D>("Probepass2JetEB","Invariant mass when probe passes + 2 Jets and it is contained in EB", 60, 60.0, 120.0);
  probepass3jetEB = fs->make<TH1D>("Probepass3JetEB","Invariant mass when probe passes + 3 Jets and it is contained in EB", 60, 60.0, 120.0);
  probepass4jetEB = fs->make<TH1D>("Probepass4JetEB","Invariant mass when probe passes + 4 Jets and it is contained in EB", 60, 60.0, 120.0);

  probefailed0jet = fs->make<TH1D>("Probefailed0Jet","Invariant mass when probe failed + no Jet", 60, 60.0, 120.0);
  probefailed1jet = fs->make<TH1D>("Probefailed1Jet","Invariant mass when probe failed + 1 Jet", 60, 60.0, 120.0);
  probefailed2jet = fs->make<TH1D>("Probefailed2Jet","Invariant mass when probe failed + 2 Jets", 60, 60.0, 120.0);
  probefailed3jet = fs->make<TH1D>("Probefailed3Jet","Invariant mass when probe failed + 3 Jets", 60, 60.0, 120.0);
  probefailed4jet = fs->make<TH1D>("Probefailed4Jet","Invariant mass when probe failed + 4 Jets", 60, 60.0, 120.0);

  probefailed0jetEE = fs->make<TH1D>("Probefailed0JetEE","Invariant mass when probe failed + no Jet and it is contained in EE", 60, 60.0, 120.0);
  probefailed1jetEE = fs->make<TH1D>("Probefailed1JetEE","Invariant mass when probe failed + 1 Jet and it is contained in EE", 60, 60.0, 120.0);
  probefailed2jetEE = fs->make<TH1D>("Probefailed2JetEE","Invariant mass when probe failed + 2 Jets and it is contained in EE", 60, 60.0, 120.0);
  probefailed3jetEE = fs->make<TH1D>("Probefailed3JetEE","Invariant mass when probe failed + 3 Jets and it is contained in EE", 60, 60.0, 120.0);
  probefailed4jetEE = fs->make<TH1D>("Probefailed4JetEE","Invariant mass when probe failed + 4 Jets and it is contained in EE", 60, 60.0, 120.0);
  probefailed0jetEB = fs->make<TH1D>("Probefailed0JetEB","Invariant mass when probe failed + no Jet and it is contained in EB", 60, 60.0, 120.0);
  probefailed1jetEB = fs->make<TH1D>("Probefailed1JetEB","Invariant mass when probe failed + 1 Jet and it is contained in EB", 60, 60.0, 120.0);
  probefailed2jetEB = fs->make<TH1D>("Probefailed2JetEB","Invariant mass when probe failed + 2 Jets and it is contained in EB", 60, 60.0, 120.0);
  probefailed3jetEB = fs->make<TH1D>("Probefailed3JetEB","Invariant mass when probe failed + 3 Jets and it is contained in EB", 60, 60.0, 120.0);
  probefailed4jetEB = fs->make<TH1D>("Probefailed4JetEB","Invariant mass when probe failed + 4 Jets and it is contained in EB", 60, 60.0, 120.0);

  probeinEB = fs->make<TH1D>("ProbeinEB","Invariant mass when probe is in EB", 60, 60.0, 120.0);
  probeinEE = fs->make<TH1D>("ProbeinEE","Invariant mass when probe is in EE", 60, 60.0, 120.0);
  probept= fs->make<TH1D>("Probept","Pt of the electron probe", 120, 0, 120.0);
  tagpt= fs->make<TH1D>("tagpt","Pt of the electron tag", 120, 0, 120.0);
  probept_0_35= fs->make<TH1D>("probept_0_35","Invariant Z mass, when probe pt < 35 GeV",60,60,120);
  probept_35_45= fs->make<TH1D>("probept_35_45","Invariant Z mass, when probe pt > 35 && pt < 45 GeV",60,60,120);
  probept_45_55= fs->make<TH1D>("probept_45_55","Invariant Z mass, when probe pt > 45 && pt < 55 GeV",60,60,120);
  probept_55_inf= fs->make<TH1D>("probept_55_inf","Invariant Z mass, when probe pt > 55",60,60,120);

  //Plot HLT efficiency...
  HLTEfficiency= fs->make<TH1D>("HLTEfficiency","HLT Efficiency: # of electron matching the HLT",3,0,3);
  HLTEfficiency0Jet= fs->make<TH1D>("HLTEfficiency0Jet","HLT Efficiency: # of electron matching the HLT + 0 Jet",3,0,3);
  HLTEfficiency1Jet= fs->make<TH1D>("HLTEfficiency1Jet","HLT Efficiency: # of electron matching the HLT + 1 Jet",3,0,3);
  HLTEfficiency2Jet= fs->make<TH1D>("HLTEfficiency2Jet","HLT Efficiency: # of electron matching the HLT + 2 Jet",3,0,3);
  HLTEfficiency3Jet= fs->make<TH1D>("HLTEfficiency3Jet","HLT Efficiency: # of electron matching the HLT + 3 Jet",3,0,3);
  HLTEfficiency4Jet= fs->make<TH1D>("HLTEfficiency4Jet","HLT Efficiency: # of electron matching the HLT + 4 Jet",3,0,3);
}



EfficiencyFilter::~EfficiencyFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}



