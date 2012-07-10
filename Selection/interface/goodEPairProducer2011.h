#ifndef goodEPairProducer2011_h_
#define goodEPairProducer2011_h_

#include <memory>
#include <string>
#include <cmath>
#include <iostream>

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "VecBosSelection/Selection/interface/SelectionUtils.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1.h"


class goodEPairProducer2011 : public edm::EDProducer, public SelectionUtils {
	public:
		explicit goodEPairProducer2011(const edm::ParameterSet &);
		~goodEPairProducer2011();

		virtual void produce(edm::Event&, edm::EventSetup const&);
		virtual void beginJob();
		//virtual bool beginRun(edm::Run &, edm::EventSetup const&);

	private:
		virtual void endJob() ;

		// ----------member data ---------------------------

		edm::InputTag theElectronCollectionLabel;
		edm::InputTag pflowEleCollection_;
		edm::InputTag pflowMuCollection_;
		edm::InputTag               conversionsInputTag_;
		edm::InputTag               beamSpotInputTag_;
		edm::InputTag               primaryVertexInputTag_;
		std::vector<edm::InputTag>  isoValInputTags_;
		edm::InputTag               ZmumuCandidates_;
		bool useNewID_;
		bool doIsolation_;
		bool doID_;
		bool doWP90_;
		bool isElectron_;

		double secondEleEnThrhold_; 
		double firstEleEnThrhold_;
		double lowZmassLimit_;
		double highZmassLimit_;
		double maxEtaForElectron_;

		bool Debug2; //Activate with true if you wonna have verbosity for debug
		TH1I* passIDEleCriteria;
		TH1I* eleSelStepByStep;
		TH1D* eventMultip;
		TH1F* gsfelEt;
		TH1D* Conversion;
		TH1D* Isolation;
		TH1D* Identification;
		TH1D* Selected;

		TH1D* eventAccept;
		TH1D* h_electronEn;
		TH1D* h_electronEta;
		TH1D* h_electronPt;
		TH1D* h_electronInvMass;
		TH1D* h_electronInvMassPass;
		TH1D* h_zPt_3e;
		//TH1D* h_zPt1j_3e;
		bool passSelection;
 
};



//
// constructors and destructor
//
goodEPairProducer2011::goodEPairProducer2011 (const edm::ParameterSet & parameters)
{
   Debug2 = false;
   passSelection=true;
   theElectronCollectionLabel = parameters.getParameter <edm::InputTag> ("electronCollection");
   pflowEleCollection_ = parameters.getUntrackedParameter<edm::InputTag>("pflowEleCollection",edm::InputTag("pfNoPileUp"));
   pflowMuCollection_ = parameters.getUntrackedParameter<edm::InputTag>("pflowMuCollection",edm::InputTag("pfNoPileUp"));
 
  ZmumuCandidates_ = parameters.getUntrackedParameter<edm::InputTag>("ZmumuCandidates");

   useNewID_ = parameters.getParameter<bool>("useNewID");
   doIsolation_          = parameters.getUntrackedParameter<bool>("doIsolation",true);
   doID_                 = parameters.getUntrackedParameter<bool>("doID",true);
   doWP90_               = parameters.getUntrackedParameter<bool>("doWP90",false);
   secondEleEnThrhold_   = parameters.getParameter<double>("secondEleEnThrhold");
   firstEleEnThrhold_    = parameters.getParameter<double>("firstEleEnThrhold");
   lowZmassLimit_        = parameters.getParameter<double>("lowZmassLimit");
   highZmassLimit_       = parameters.getParameter<double>("highZmassLimit");
   maxEtaForElectron_    = parameters.getParameter<double>("maxEtaForElectron");
   conversionsInputTag_  = parameters.getParameter<edm::InputTag>("conversionsInputTag");
   beamSpotInputTag_     = parameters.getParameter<edm::InputTag>("beamSpotInputTag");
   primaryVertexInputTag_= parameters.getParameter<edm::InputTag>("primaryVertexInputTag");
   isoValInputTags_      = parameters.getParameter<std::vector<edm::InputTag> >("isoValInputTags");
   isElectron_            = parameters.getUntrackedParameter<bool>("isElectron");

   produces<reco::PFCandidateCollection>();
   

   //Initializations...
   edm::Service<TFileService> fs;

   passIDEleCriteria = fs->make<TH1I>("passIDEleCriteria","Ele Id pass/not pass... 3 entries each ele", 4, 0, 4);
   eleSelStepByStep = fs->make<TH1I>("eleSelStepByStep","History of selected/rejected ele", 13, 0, 13);
   
   eventAccept= fs->make<TH1D>("eventAccept","There is corresponding PF", 20, 0, 20);
   h_electronEn= fs->make<TH1D>("h_electronEn","electronEn", 200, 0, 200);
   h_electronEta= fs->make<TH1D>("h_electronEta","electronEta", 200, -2.5, 2.5);
   h_electronPt= fs->make<TH1D>("h_electronPt","electronPt", 200, 0, 200);
   h_electronInvMass= fs->make<TH1D>("h_electronInvMass","electronInvMass", 60, 60, 120);
   h_electronInvMassPass= fs->make<TH1D>("h_electronInvMassPass","electronInvMassPass", 60, 60, 120);
   h_zPt_3e = fs->make<TH1D>("h_zPt_3e","zPt_3e", 200, 0, 200);
   //h_zPt1j_3e = fs->make<TH1D>("h_zPt1j_3e","zPt1j_3e", 200, 0, 200);
}



goodEPairProducer2011::~goodEPairProducer2011 ()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

#endif


