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
// Original Author:  Vieri Candelise & Matteo Marone
//         Created:  Wed May 11 14:53:26 CEST 2011
// $Id: ZanalyzerFilter.cc,v 1.12 2011/11/09 10:18:43 marone Exp $
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

#include "TString.h"

using namespace std;
using namespace edm;
using namespace reco;


bool hltispresent=1;
bool davdebug=0;


//
// member functions
//

// ------------ method called for each event  ------------
bool
ZanalyzerFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
	if (debug) cout<<"------- NEW Event -----"<<endl;
	using namespace edm;
	Handle < GsfElectronCollection > electronCollection;
	iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
	if (!electronCollection.isValid ())
		return false;


	//Match The HLT Trigger
	using edm::TriggerResults;
	Handle<TriggerResults> HLTResults;
	iEvent.getByLabel(triggerCollection_, HLTResults);
	const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTResults);
	bool flag=false;

	if (HLTResults.isValid() && doTheHLTAnalysis_) {
		/// Storing the Prescale information: loop over the triggers and record prescale

		// Matching the HLT information event per event if no hlt info is present in iRun
		if(hltispresent==false){
			std::vector<std::string>  hlNames;
			hlNames.clear();
			hlNames=triggerNames.triggerNames();
			triggerIndices_.clear();


			for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
				if(find(hlNames.begin(),hlNames.end(),triggerNames_[itrig])!=hlNames.end()){
					triggerIndices_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
				}
				else{
					triggerIndices_.push_back(2048);
				}
			}


		}

		unsigned int minimalPrescale(10000);
		unsigned int prescale(0);
		bool bit(true);
		std::pair<int,int> prescalepair;
		std::vector<int>  triggerSubset;
		for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
			if(triggerIndices_[itrig]!=2048) {

				// check trigger response
				bit = HLTResults->accept(triggerIndices_[itrig]);
				triggerSubset.push_back(bit);

				if(bit) {
				  flag=true;
cout << "il valore della flag " << flag << " should be quite always 1 \n";
				  //If path is accepted, then together with its prescale it is stored in a map.
						int prescaleset = hltConfig_.prescaleSet(iEvent,iSetup);
						if(prescaleset!=-1) {
							prescalepair = hltConfig_.prescaleValues(iEvent,iSetup,triggerNames_[itrig]);
							if (debug) cout<<"prescale.first "<<prescalepair.first<<" prescalepair.second "<<prescalepair.second<<endl;
							//getting prescale info
							prescale = useCombinedPrescales_ ? prescalepair.first*prescalepair.second : prescalepair.second;
							if((useCombinedPrescales_ && prescalepair.first<0) || prescalepair.second<0) {
								edm::LogWarning("MyAnalyzer") << " Unable to get prescale from event for trigger " << triggerNames.triggerName(itrig) << " :" << prescalepair.first << ", " << prescalepair.second;
								prescale = -999;
							}

							if(prescalepair.first<0 || prescalepair.second<0) { prescale = -999; }
						}

						minimalPrescale = minimalPrescale <  prescale ? minimalPrescale : prescale;
						if (debug) cout<<"prescale "<<prescale<<" minimal Prescale "<<minimalPrescale<<" for trigger "<<triggerNames.triggerName(itrig)<<endl;


				} //Chiusura del if(bit)
				else {
					//edm::LogError("HistoProducer") << " Unable to get prescale set from event. Check that L1 data products are present.";
				}
			}
			else {
				// that trigger is presently not in the menu
				triggerSubset.push_back(false);
			}
		} //chiusura for
	} //chiusura HLT studies

	if (!flag) 
	  {
	    if(!useAllTriggers_) return false;

	  }
 

  //===========================

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


  for (reco::GsfElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {

    if (recoElectron->et () <= 25)  continue;

    ///////////////

    // Define Isolation variables
    double IsoTrk = (recoElectron->dr03TkSumPt () / recoElectron->et ());
    double IsoEcal = (recoElectron->dr03EcalRecHitSumEt () / recoElectron->et ());
    double IsoHcal = (recoElectron->dr03HcalTowerSumEt () / recoElectron->et ());
    double HE = recoElectron->hadronicOverEm();

    //Define ID variables

    float DeltaPhiTkClu = recoElectron->deltaPhiSuperClusterTrackAtVtx ();
    float DeltaEtaTkClu = recoElectron->deltaEtaSuperClusterTrackAtVtx ();
    float sigmaIeIe = recoElectron->sigmaIetaIeta ();

    //Define Conversion Rejection Variables

    float Dcot = recoElectron->convDcot ();
    float Dist = recoElectron->convDist ();
    int NumberOfExpectedInnerHits = recoElectron->gsfTrack ()->trackerExpectedHitsInner ().numberOfHits ();

    if (removePU_){
      double lepIsoRho;
		  
      /////// Pileup density "rho" for lepton isolation subtraction /////
      edm::Handle<double> rhoLepIso;
      const edm::InputTag eventrhoLepIso("kt6PFJetsForIsolation", "rho");
      iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
      if( *rhoLepIso == *rhoLepIso) { 
	lepIsoRho = *rhoLepIso;
      }
      else { 
	lepIsoRho =  999999.9;
      }
		  
      //EB
      if (fabs (recoElectron->eta()) <= 1.4442) {      
	//
	IsoTrk = (recoElectron->dr03TkSumPt () - lepIsoRho*0) / recoElectron->et ();
	IsoEcal = (recoElectron->dr03EcalRecHitSumEt () - lepIsoRho*0.096) / recoElectron->et ();
	IsoHcal = (recoElectron->dr03HcalTowerSumEt ()  - lepIsoRho*0.020) / recoElectron->et ();
	if(IsoEcal<=0.) IsoEcal=0.;
	if(IsoHcal<=0.) IsoHcal=0.;
      }
      //EE
      if (fabs (recoElectron->eta()) >= 1.5660
	  && fabs (recoElectron->eta()) <= 2.5000) {
	//
	IsoTrk = (recoElectron->dr03TkSumPt () - lepIsoRho*0) / recoElectron->et ();
	IsoEcal = (recoElectron->dr03EcalRecHitSumEt () - lepIsoRho*0.044) / recoElectron->et ();
	IsoHcal = (recoElectron->dr03HcalTowerSumEt ()  - lepIsoRho*0.041) / recoElectron->et ();
	if(IsoEcal<=0.) IsoEcal=0.;
	if(IsoHcal<=0.) IsoHcal=0.;
      }
    }
		

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
      //cout<<"isIsolatedBarrel "<<isIsolatedBarrel<<" isIDBarrel "<< isIDBarrel<<" isConvertedBarrel "<<isConvertedBarrel<<endl;
    }

    if (isIsolatedBarrel && isIDBarrel && isConvertedBarrel) {
      elIsAccepted++;
      elIsAcceptedEB++;
      TLorentzVector b_e2;
      b_e2.SetPtEtaPhiM(recoElectron->pt(),recoElectron->eta(),recoElectron->phi(), 0.0);
      //      TLorentzVector b_e2(recoElectron->momentum ().x (),recoElectron->momentum ().y (),recoElectron->momentum ().z (), recoElectron->p ());
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
      //cout<<"isIsolatedEndcap "<<isIsolatedEndcap<<" isIDEndcap "<< isIDEndcap<<" isConvertedEndcap "<<isConvertedEndcap<<endl;
    }

    if (isIsolatedEndcap && isIDEndcap && isConvertedEndcap) {
      elIsAccepted++;
      elIsAcceptedEE++;
      TLorentzVector e_e2;
      e_e2.SetPtEtaPhiM(recoElectron->pt(),recoElectron->eta(),recoElectron->phi(), 0.0);
      LV.push_back(e_e2);
    }

  }
  if (elIsAccepted<=1)    return false;
  double e_ee_invMass=0; 
  if (elIsAccepted>2) cout<<"WARNING: In this events we have more than two electrons accpeted!!!!!!!"<<endl;
  if (LV.size()==2){
    TLorentzVector e_pair = LV[0] + LV[1];
    e_ee_invMass = e_pair.M ();
    h_invMass->Fill(e_ee_invMass);
  }  

  if (elIsAcceptedEB==2){
    h_invMassBB->Fill(e_ee_invMass);
  }
  if (elIsAcceptedEE==2){
    h_invMassEE->Fill(e_ee_invMass);
  }
  if (elIsAcceptedEB==1 && elIsAcceptedEE==1){
    h_invMassEB->Fill(e_ee_invMass);
  }

  eventAccept->Fill(elIsAccepted);
  LV.clear();

  return true;

}


// ------------ method called once each job just before starting event loop  ------------
void
ZanalyzerFilter::beginJob (){

  //beginJob

}


// ------------ method called when starting to processes a run  ------------
bool
ZanalyzerFilter::beginRun(edm::Run &iRun, edm::EventSetup const& iSetup)
{

hltispresent=true;
	//HLT names
	std::vector<std::string>  hlNames;
	hlNames.clear();
	bool changed (true);
	if (hltConfig_.init(iRun,iSetup,triggerCollection_.process(),changed)) {
		if (changed) {
			hlNames = hltConfig_.triggerNames();
		}
	} else {
		edm::LogError("MyAnalyzer") << " HLT config extraction failure with process name " << triggerCollection_.process();
	}

	//debug dav
	if(hlNames.size()==0) { 
		hltispresent=false;
	}

	if (debug) cout<<"useAllTriggers? "<<useAllTriggers_<<endl;
	if(useAllTriggers_) triggerNames_ = hlNames;
	//triggerNames_ = hlNames;
	//HLT indices
	triggerIndices_.clear();


unsigned int myflag=0;

	for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
		if(find(hlNames.begin(),hlNames.end(),triggerNames_[itrig])!=hlNames.end()){
			triggerIndices_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
		}
		else{
			triggerIndices_.push_back(2048);
			myflag++;
		}
	}
	return true;

	//qui sarebbe piÃ¹ intelligente tornare false se nella lista dei path hlt non c'Ã¨ quello che ci interessa cosÃ¬ eviti di andare comunque ogni evento in cerca di lui (ovvio devi verificare che doTheHLTanalysis sia true qui!)
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZanalyzerFilter::endJob ()
{

  //endJob

}

DEFINE_FWK_MODULE (ZanalyzerFilter);
