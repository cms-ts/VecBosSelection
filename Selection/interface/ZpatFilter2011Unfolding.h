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
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "VecBosSelection/Selection/interface/SelectionUtils.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class ZpatFilter2011Unfolding : public edm::EDFilter, public SelectionUtils {
	public:
		explicit ZpatFilter2011Unfolding(const edm::ParameterSet &);
		~ZpatFilter2011Unfolding();

		virtual void beginJob();
		virtual bool beginRun(edm::Run &, edm::EventSetup const&);

	private:
		virtual bool filter(edm::Event&, edm::EventSetup const&);
		virtual void endJob() ;

		// ----------member data ---------------------------

		edm::InputTag theElectronCollectionLabel;
		edm::InputTag triggerCollection_; 
		bool useCombinedPrescales_; // switch between HLT only and L1*HLT prescales
		edm::InputTag               conversionsInputTag_;
		edm::InputTag               beamSpotInputTag_;
		edm::InputTag               primaryVertexInputTag_;
		std::vector<edm::InputTag>  isoValInputTags_;
		edm::InputTag genParticleCollection_;

		bool useAllTriggers_; // if no trigger names are provided, use all triggers to find event weight
		HLTConfigProvider hltConfig_;        // to get configuration for L1s/Pre
		std::vector<std::string> triggerNames_; // name of the algorithms selected by our analysis
		std::vector<unsigned int> triggerIndices_; // index of the algorithms selected by our analysis
		bool doTheHLTAnalysis_;
		//bool removePU_;
		bool useNewID_;
		bool doIsolation_;
		bool doID_;
		bool doWP90_;
		double secondEleEnThrhold_; 
		double firstEleEnThrhold_;
		double lowZmassLimit_;
		double highZmassLimit_;
		double maxEtaForElectron_;
		bool isUnfolding_;

		//  std::string outputFile_;
		TH1F* gsfelEt;
		TH1D* Conversion;
		TH1D* Isolation;
		TH1D* Identification;
		TH1D* Selected;
		TH1F* h_invMass;
		TH1F* h_invMassEE;
		TH1F* h_invMassEB;
		TH1F* h_invMassBB;
		TH1F* h_zPt_3e;
		TH1I* passIDEleCriteria;
		TH1I*  eleSelStepByStep;
		TH1I*  numberOfEleAfterHLTTrigger;
};



//
// constructors and destructor
//
ZpatFilter2011Unfolding::ZpatFilter2011Unfolding (const edm::ParameterSet & parameters)
{
	theElectronCollectionLabel = parameters.getParameter <edm::InputTag> ("electronCollection");
	triggerCollection_	= parameters.getParameter<edm::InputTag>("triggerCollectionTag");
	useCombinedPrescales_ = parameters.getParameter<bool>("UseCombinedPrescales");
	triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");
	useAllTriggers_       = (triggerNames_.size()==0);
	doTheHLTAnalysis_     = parameters.getParameter<bool>("doTheHLTAnalysis");
	//removePU_             = parameters.getParameter<bool>("removePU");
	useNewID_             = parameters.getParameter<bool>("useNewID");
	doIsolation_          = parameters.getUntrackedParameter<bool>("doIsolation",true);
	doID_                 = parameters.getUntrackedParameter<bool>("doID",true);
	doWP90_                 = parameters.getUntrackedParameter<bool>("doWP90",false);
	secondEleEnThrhold_   = parameters.getParameter<double>("secondEleEnThrhold");
	firstEleEnThrhold_    = parameters.getParameter<double>("firstEleEnThrhold");
	lowZmassLimit_        = parameters.getParameter<double>("lowZmassLimit");
	highZmassLimit_       = parameters.getParameter<double>("highZmassLimit");
	maxEtaForElectron_    = parameters.getParameter<double>("maxEtaForElectron");
	conversionsInputTag_  = parameters.getParameter<edm::InputTag>("conversionsInputTag");
        beamSpotInputTag_     = parameters.getParameter<edm::InputTag>("beamSpotInputTag");
        primaryVertexInputTag_= parameters.getParameter<edm::InputTag>("primaryVertexInputTag");
        isoValInputTags_      = parameters.getParameter<std::vector<edm::InputTag> >("isoValInputTags");
	isUnfolding_          = parameters.getUntrackedParameter<bool>("isUnfolding",true);
	genParticleCollection_= parameters.getUntrackedParameter<edm::InputTag>("genParticleCollection",edm::InputTag("genParticles"));
	

  //Initializations...
  edm::Service<TFileService> fs;

  h_invMass = fs->make<TH1F>("Z peak - WP80","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassEE =  fs->make<TH1F>("Z peak - WP80 Endcap-Endcap","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassEB = fs->make<TH1F>("Z peak - WP80 Endcap-Barrel","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassBB = fs->make<TH1F>("Z peak - WP80 Barrel-Barrel","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_zPt_3e = fs->make<TH1F>("h_zPt_3e","zPt_3e", 200, 0.0, 200.0);
  passIDEleCriteria = fs->make<TH1I>("passIDEleCriteria","Ele Id pass/not pass... 3 entries each ele", 4, 0, 4);
  eleSelStepByStep = fs->make<TH1I>("eleSelStepByStep","History of selected/rejected ele", 13, 0, 13);
  numberOfEleAfterHLTTrigger = fs->make<TH1I>("numberOfEleAfterHLTTrigger","Size of ele collection, after HLT cuts,", 3, 0, 3);
}



ZpatFilter2011Unfolding::~ZpatFilter2011Unfolding ()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


