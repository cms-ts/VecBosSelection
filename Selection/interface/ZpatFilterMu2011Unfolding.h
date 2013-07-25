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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

class ZpatFilterMu2011Unfolding : public edm::EDFilter, public SelectionUtils {
	public:
		explicit ZpatFilterMu2011Unfolding(const edm::ParameterSet &);
		~ZpatFilterMu2011Unfolding();

		virtual void beginJob();

	private:
		virtual bool filter(edm::Event&, edm::EventSetup const&);
		virtual void endJob() ;

		// ----------member data ---------------------------

		edm::InputTag theMuCollectionLabel;
		edm::InputTag genParticleCollection_;
		bool isUnfolding_;
		
		double lowZmassLimit_;
		double highZmassLimit_;

		//  std::string outputFile_;
		TH1F* h_invMass;
		TH1F* h_invMassEE;
		TH1F* h_invMassEB;
		TH1F* h_invMassBB;
		TH1F* h_zPt_3mu;
		TH1I*  muSelStepByStep;
};



//
// constructors and destructor
//
ZpatFilterMu2011Unfolding::ZpatFilterMu2011Unfolding (const edm::ParameterSet & parameters)
{
	theMuCollectionLabel = parameters.getParameter <edm::InputTag> ("muonCollection");
	lowZmassLimit_        = parameters.getParameter<double>("lowZmassLimit");
	highZmassLimit_       = parameters.getParameter<double>("highZmassLimit");
	isUnfolding_          = parameters.getUntrackedParameter<bool>("isUnfolding",true);
	genParticleCollection_= parameters.getUntrackedParameter<edm::InputTag>("genParticleCollection",edm::InputTag("genParticles"));
	
  //Initializations...
  edm::Service<TFileService> fs;

  h_invMass = fs->make<TH1F>("Z peak - WP80","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassEE =  fs->make<TH1F>("Z peak - Endcap-Endcap","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassEB = fs->make<TH1F>("Z peak - Endcap-Barrel","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_invMassBB = fs->make<TH1F>("Z peak - Barrel-Barrel","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
  h_zPt_3mu = fs->make<TH1F>("h_zPt_3mu","zPt_3mu", 200, 0.0, 200.0);
  muSelStepByStep = fs->make<TH1I>("muSelStepByStep","History of selected/rejected mu", 13, 0, 13);
}



ZpatFilterMu2011Unfolding::~ZpatFilterMu2011Unfolding ()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


