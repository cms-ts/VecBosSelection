// -*- C++ -*-
//
// Package:    Zanalyzer
// Class:      Zanalyzer
// 
/**\class Zanalyzer Zanalyzer.cc Zmonitoring/Zanalyzer/src/Zanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vieri Candelise
//         Created:  Wed May 11 14:53:26 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/ZanalyzerFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include <string>
#include <cmath>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"


using namespace std;
using namespace edm;
using namespace reco;


//
// constructors and destructor
//
ZanalyzerFilter::ZanalyzerFilter (const edm::ParameterSet & parameters)
{
  theElectronCollectionLabel =
    parameters.getParameter < InputTag > ("electronCollection");
  std::string outputFile_D = parameters.getUntrackedParameter<std::string>("filename");
  outputFile_ = parameters.getUntrackedParameter<std::string>("outputFile", outputFile_D);
  triggerCollectionTag_=parameters.getUntrackedParameter<edm::InputTag>("triggerCollectionTag");
}



ZanalyzerFilter::~ZanalyzerFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}




//
// member functions
//

// ------------ method called for each event  ------------
bool
ZanalyzerFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{

  using namespace edm;
  Handle < GsfElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ())
    return false;

  //Match The HLT Trigger
  edm::Handle<edm::TriggerResults> HLTResults;
  iEvent.getByLabel(triggerCollectionTag_, HLTResults);
  if (!HLTResults.isValid ())
    return false;

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTResults);   
  int trigger_size = HLTResults->size();
  bool flag=false;

  string path[10]={"HLT_Ele10_LW_L1R","HLT_Ele15_SW_L1R","HLT_Ele15_SW_CaloEleId_L1R","HLT_Ele17_SW_CaloEleId_L1R","HLT_Ele17_SW_TightEleId_L1R","HLT_Ele17_SW_TightEleId_L1R_v2","HLT_Ele17_SW_TightEleId_L1R_v3","HLT_Photon10_L1R","HLT_Photon15_L1R","HTL_Photon15_Cleaned_L1R"};
  
  for (int i=0; i<10;i++){
    int pos=(int)triggerNames.triggerIndex(path[i]);
    if(pos<trigger_size){
      path[i]=(int) HLTResults->accept(pos);
      if (path[i]==1) {
	flag=true;
	cout<<"Matched "<<path[i]<<endl;
      }
    }
  }

  if (!flag) 
    {
      //cout<<"It does not match with HLT triggers, sorry"<<endl;
      //return false;
    }
 
  bool isBarrelElectrons;
  bool isEndcapElectrons;
  bool isIsolatedBarrel;
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIsolatedEndcap;
  bool isIDEndcap;
  bool isConvertedEndcap;
  int elIsAccepted=0;
  int elIsAcceptedEB=0;
  int elIsAcceptedEE=0;

  std::vector<TLorentzVector> LV;

  //Check the numb of electrons in the event
  //eventMultip->Fill(electronCollection->size());

  for (reco::GsfElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    // if (electronCollection->size()==1) cout<<"One electron event has energy of "<<recoElectron->et()<<endl;

    if (recoElectron->et () <= 25)  continue;
 
    // Define Isolation variables
    double IsoTrk = (recoElectron->dr03TkSumPt () / recoElectron->et ());
    double IsoEcal = (recoElectron->dr03EcalRecHitSumEt () / recoElectron->et ());
    double IsoHcal = (recoElectron->dr03HcalTowerSumEt () / recoElectron->et ());
    double HE = (recoElectron->hcalOverEcal ());

    //Define ID variables

    float DeltaPhiTkClu = recoElectron->deltaPhiSuperClusterTrackAtVtx ();
    float DeltaEtaTkClu = recoElectron->deltaEtaSuperClusterTrackAtVtx ();
    float sigmaIeIe = recoElectron->sigmaIetaIeta ();

    //Define Conversion Rejection Variables

    float Dcot = recoElectron->convDcot ();
    float Dist = recoElectron->convDist ();
    int NumberOfExpectedInnerHits = recoElectron->gsfTrack ()->trackerExpectedHitsInner ().numberOfHits ();

    //quality flags

    isBarrelElectrons = false;
    isEndcapElectrons = false;
    isIsolatedBarrel = false;
    isIDBarrel = false;
    isConvertedBarrel = false;
    isIsolatedEndcap = false;
    isIDEndcap = false;
    isConvertedEndcap = false;
 
  /***** Barrel WP80 Cuts *****/

    if (fabs (recoElectron->eta ()) <= 1.4442) {

      /* Isolation */
      if (IsoTrk < 0.09 && IsoEcal < 0.07 && IsoHcal < 0.10) {
	isIsolatedBarrel = true;
      }

      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.004 && fabs (DeltaPhiTkClu) < 0.06
	  && sigmaIeIe < 0.01 && HE < 0.04) {
	isIDBarrel = true;
      }

      /* Conversion Rejection */
      if ((fabs (Dist) >= 0.02 || fabs (Dcot) >= 0.02)
	  && NumberOfExpectedInnerHits == 0) {
	isConvertedBarrel = true;
      }

    }

    if (isIsolatedBarrel && isIDBarrel && isConvertedBarrel) {
      elIsAccepted++;
      elIsAcceptedEB++;
      TLorentzVector b_e2(recoElectron->momentum ().x (),recoElectron->momentum ().y (),recoElectron->momentum ().z (), recoElectron->p ());
      LV.push_back(b_e2);
    }

  /***** Endcap WP80 Cuts *****/

    if (fabs (recoElectron->eta ()) >= 1.5660
	&& fabs (recoElectron->eta ()) <= 2.5000) {

      /* Isolation */
      if (IsoTrk < 0.04 && IsoEcal < 0.05 && IsoHcal < 0.025) {
	isIsolatedEndcap = true;
      }

      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.007 && fabs (DeltaPhiTkClu) < 0.03
	  && sigmaIeIe < 0.03 && HE < 0.15) {
	isIDEndcap = true;
      }

      /* Conversion Rejection */
      if ((fabs (Dcot) > 0.02 || fabs (Dist) > 0.02)
	  && NumberOfExpectedInnerHits == 0) {
	isConvertedEndcap = true;
      }
    }

    if (isIsolatedEndcap && isIDEndcap && isConvertedEndcap) {
      elIsAccepted++;
      elIsAcceptedEE++;
       TLorentzVector e_e2(recoElectron->momentum ().x (),recoElectron->momentum ().y (),recoElectron->momentum ().z (), recoElectron->p ());
      LV.push_back(e_e2);
    }

  }
  if (elIsAcceptedEB<=1)    return false;
  TLorentzVector e_pair = LV[0] + LV[1];
  double e_ee_invMass = e_pair.M ();
  if (LV.size()==2) h_invMass->Fill(e_ee_invMass);  
  if (elIsAccepted>2) cout<<"In this events we have more than two electrons accpeted!!!!!!!"<<endl;
  eventAccept->Fill(elIsAccepted);
  
  return true;
  
  // end of reco electron loop
}


// ------------ method called once each job just before starting event loop  ------------
void
ZanalyzerFilter::beginJob (){
  //const float pi = 3.14159265;
      //h_ee_invMass_EB = 0;
      //h_ee_invMass_EE = 0;
      //h_ee_invMass_BB = 0;
  fOFile = new TFile("ZAnalysisFilter.root","RECREATE");
  //eventMultip=new TH1D("eventMultip","Event Multiplicity", 20, 0, 20);
  eventAccept=new TH1D("eventAccept","Good Event Multiplicity", 20, 0, 20);
  //gsfelEt=new TH1F("gsfelEt","gsf electron tEnergy ", 100, 0, 100);
  //Conversion=new TH1D("Conversion","NumberOfElectron Not Converted", 3, 0, 3);
  //Isolation=new TH1D("Isolation","NumberOfElectron Isolated", 3, 0, 3);
  //Identification=new TH1D("Identification","NumberOfElectron good ID", 3, 0, 3);
  //Selected=new TH1D("Selected","Number Of Selected electrons", 3, 0, 3);
  //h_bb_invMass_BB=new TH1F("Z peak - WP80 EB-EB","Z peak - WP80 EB-EB;InvMass (Gev)", 60, 60.0, 120.0);
  //h_ee_invMass_EB =new TH1F("Z peak - WP80 EB-EE","Z peak - WP80 EB-EE;InvMass (GeV)", 60, 60.0, 120.0);
  //h_ee_invMass_EE =new TH1F("Z peak - WP80 EE-EE","Z peak - WP80 EE-EE;InvMass (Gev)", 60, 60.0, 120.0);
  h_invMass =new TH1F("Z peak - WP80","Z peak;InvMass (Gev)", 60, 60.0, 120.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZanalyzerFilter::endJob ()
{
  fOFile->cd();
  //eventMultip->Write();
  eventAccept->Write();
  //gsfelEt->Write();
  //Conversion->Write();
  //Isolation->Write();
  //Identification->Write();
  //Selected->Write();
  //h_bb_invMass_BB->Write();
  //h_ee_invMass_EB->Write();
  //h_ee_invMass_EE->Write();
  h_invMass->Write();

  fOFile->Write() ;
  fOFile->Close() ;
}

DEFINE_FWK_MODULE (ZanalyzerFilter);
