#ifndef analyzerMuCuts_h
#define analyzerMuCuts_h

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

class analyzerMuCuts : public edm::EDAnalyzer{
	public:
		explicit analyzerMuCuts(const edm::ParameterSet &);
		~analyzerMuCuts();

		virtual void beginJob();

	private:
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		// ----------member data ---------------------------

		edm::InputTag theMuCollectionLabel;
		edm::InputTag theMuMatchedCollectionLabel;
		TH1I*  muSelStepByStep;
		//  std::string outputFile_;
		//TH1F* h_invMass;
};



//
// constructors and destructor
//
analyzerMuCuts::analyzerMuCuts (const edm::ParameterSet & parameters)
{
	theMuCollectionLabel = parameters.getParameter <edm::InputTag> ("muonCollection");
	theMuMatchedCollectionLabel = parameters.getParameter <edm::InputTag> ("muonMatchedCollection");

  //Initializations...
  edm::Service<TFileService> fs;
  muSelStepByStep = fs->make<TH1I>("muSelStepByStep","History of selected/rejected mu", 13, 0, 13);
  //h_invMass = fs->make<TH1F>("Z peak - WP80","Z peak;InvMass (Gev)", 140, 0.0, 140.0);
}



analyzerMuCuts::~analyzerMuCuts ()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

#endif
