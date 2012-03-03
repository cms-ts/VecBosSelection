#ifndef ZAnalyzerProducer_h_
#define ZAnalyzerProducer_h_

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


class ZanalyzerProducer : public edm::EDProducer, public SelectionUtils {
	public:
		explicit ZanalyzerProducer(const edm::ParameterSet &);
		~ZanalyzerProducer();

		virtual void produce(edm::Event&, edm::EventSetup const&);
		virtual void beginJob();
		//virtual bool beginRun(edm::Run &, edm::EventSetup const&);

	private:
		virtual void endJob() ;

		// ----------member data ---------------------------

		edm::InputTag theElectronCollectionLabel;
		bool removePU_;
		bool Debug2; //Activate with true if you wonna have verbosity for debug
		TH1I* passIDEleCriteria;
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
		bool passSelection;
 
};



//
// constructors and destructor
//
ZanalyzerProducer::ZanalyzerProducer (const edm::ParameterSet & parameters)
{
   Debug2 = false;
   passSelection=true;
   theElectronCollectionLabel = parameters.getParameter <edm::InputTag> ("electronCollection");
   removePU_ = parameters.getParameter<bool>("removePU");
   produces<reco::GsfElectronCollection>();
   

   //Initializations...
   edm::Service<TFileService> fs;

   passIDEleCriteria = fs->make<TH1I>("passIDEleCriteria","Ele Id pass/not pass... 3 entries each ele", 3, 0, 3);
   
   eventAccept= fs->make<TH1D>("eventAccept","There is corresponding GSF", 20, 0, 20);
   h_electronEn= fs->make<TH1D>("h_electronEn","electronEn", 200, 0, 200);
   h_electronEta= fs->make<TH1D>("h_electronEta","electronEta", 200, -2.5, 2.5);
   h_electronPt= fs->make<TH1D>("h_electronPt","electronPt", 200, 0, 200);
   h_electronInvMass= fs->make<TH1D>("h_electronInvMass","electronInvMass", 60, 60, 120);
   h_electronInvMassPass= fs->make<TH1D>("h_electronInvMassPass","electronInvMassPass", 60, 60, 120);

}



ZanalyzerProducer::~ZanalyzerProducer ()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

#endif

