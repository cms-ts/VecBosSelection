#ifndef photonRemoval_h_
#define photonRemoval_h_

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


class photonRemoval : public edm::EDProducer{
	public:
		explicit photonRemoval(const edm::ParameterSet &);
		~photonRemoval();

		virtual void produce(edm::Event&, edm::EventSetup const&);
		virtual void beginJob();
		//virtual bool beginRun(edm::Run &, edm::EventSetup const&);

	private:
		virtual void endJob() ;		
		template <class T> bool searchHistory(T const&);

		// ----------member data ---------------------------

		edm::InputTag particleCollectionLabel;
		edm::InputTag bottomCollectionLabel;

		double barrelRCone;
		double endcapRCone;
		TLorentzVector e1_test;
		double minZMass;
		double maxZMass;

		TH1F * gammaRemovedPt;
		TH1F * gammaRemovedEta;
		TH1F * eRemovedPt;
		TH1F * eRemovedEta;
		TH1F * eNotRemovedMass;
		TH1F * eNotRemovedPt;
		TH1F * eNotRemovedEta;
 
};



//
// constructors and destructor
//
photonRemoval::photonRemoval (const edm::ParameterSet & parameters)
{
   particleCollectionLabel = parameters.getParameter <edm::InputTag> ("particleCollection");
   bottomCollectionLabel   = parameters.getParameter <edm::InputTag> ("bottomCollection");
   barrelRCone             = parameters.getUntrackedParameter<double>("barrelRCone",0.05);
   endcapRCone             = parameters.getUntrackedParameter<double>("endcapRCone",0.07);
   produces<reco::GenParticleCollection>();
   
   minZMass = -1;
   maxZMass = 9999;
   
   //Initializations...
   edm::Service<TFileService> fs;

   gammaRemovedPt = fs->make<TH1F>("gammaRemovedPt","gammaRemovedPt", 100, 0, 50);
   gammaRemovedEta = fs->make<TH1F>("gammaRemovedEta","gammaRemovedEta", 100, -3.0, 3.0);
   eRemovedPt = fs->make<TH1F>("eRemovedPt","eRemovedPt", 100, 0, 100);
   eRemovedEta = fs->make<TH1F>("eRemovedEta","eRemovedEta", 100, -3.0, 3.0);
   eNotRemovedMass = fs->make<TH1F>("eNotRemovedMass","eNotRemovedMass",200,0,200);
   eNotRemovedPt = fs->make<TH1F>("eNotRemovedPt","eNotRemovedPt",200,0,200);
   eNotRemovedEta = fs->make<TH1F>("eNotRemovedEta","eNotRemovedEta", 100, -3.0, 3.0);
}



photonRemoval::~photonRemoval ()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

template <class T>
bool photonRemoval::searchHistory(T const& itSearch)
{
   bool pass = false;
   if (fabs(itSearch->pdgId()) ==11 && itSearch->status()==3 && itSearch->mother()->pdgId()==23) {
      e1_test.SetPtEtaPhiM(itSearch->pt(),itSearch->eta(),
                           itSearch->phi(),itSearch->mass());

      pass=true;
      return pass;
   } else if (itSearch->status()==3 ){
      pass=false;
      return pass;
   } else {
      pass=searchHistory(itSearch->mother());
   }
   return pass;
}


#endif


