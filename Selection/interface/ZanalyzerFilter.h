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

bool debug=false; //Activate with true if you wonna have verbosity for debug

class ZanalyzerFilter : public edm::EDFilter {
	public:
		explicit ZanalyzerFilter(const edm::ParameterSet &);
		~ZanalyzerFilter();

		virtual void beginJob();
		virtual bool beginRun(edm::Run &, edm::EventSetup const&);

	private:
		virtual bool filter(edm::Event&, edm::EventSetup const&);
		virtual void endJob() ;

		// ----------member data ---------------------------

		edm::InputTag theElectronCollectionLabel;
		edm::InputTag triggerCollection_; 
		bool useCombinedPrescales_; // switch between HLT only and L1*HLT prescales
		bool useAllTriggers_; // if no trigger names are provided, use all triggers to find event weight
		HLTConfigProvider hltConfig_;        // to get configuration for L1s/Pre
		std::vector<std::string> triggerNames_; // name of the algorithms selected by our analysis
		std::vector<unsigned int> triggerIndices_; // index of the algorithms selected by our analysis
		bool doTheHLTAnalysis_;
		bool removePU_;


		//  std::string outputFile_;
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


};



//
// constructors and destructor
//
ZanalyzerFilter::ZanalyzerFilter (const edm::ParameterSet & parameters)
{
	theElectronCollectionLabel = parameters.getParameter <edm::InputTag> ("electronCollection");
	triggerCollection_	= parameters.getParameter<edm::InputTag>("triggerCollectionTag");
	useCombinedPrescales_ = parameters.getParameter<bool>("UseCombinedPrescales");
	triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");
	useAllTriggers_       = (triggerNames_.size()==0);
	doTheHLTAnalysis_     = parameters.getParameter<bool>("doTheHLTAnalysis");
	removePU_             = parameters.getParameter<bool>("removePU");



  //Initializations...
  edm::Service<TFileService> fs;

  eventAccept= fs->make<TH1D>("eventAccept","Good Event Multiplicity", 20, 0, 20);
  h_invMass = fs->make<TH1F>("Z peak - WP80","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassEE =  fs->make<TH1F>("Z peak - WP80 Endcap-Endcap","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassEB = fs->make<TH1F>("Z peak - WP80 Endcap-Barrel","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassBB = fs->make<TH1F>("Z peak - WP80 Barrel-Barrel","Z peak;InvMass (Gev)", 140, 0.0, 140.0);




}



ZanalyzerFilter::~ZanalyzerFilter ()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


