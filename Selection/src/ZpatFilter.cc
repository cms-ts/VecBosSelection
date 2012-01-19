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
// Original Author:  Vieri Candelise, Matteo Marone & Davide Scaini
//         Created:  Thu Dec 11 10:46:26 CEST 2011
// $Id: ZpatFilter.cc,v 1.3 2012/01/18 10:08:00 schizzi Exp $
//
//


// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/ZpatFilter.h"
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


bool hltispresent2=1; //Necessary to correctly match the HLT
bool Debug2=false;    //Activate with true if you wonna have verbosity for Debug

  using namespace edm;

//
// member functions
//

// ------------ method called for each event  ------------
bool
ZpatFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  
  if (Debug2) cout<<"------- NEW Event -----"<<endl;
	
	//  Match The HLT Trigger
	/// to study fired HLT paths
	using edm::TriggerResults;
	Handle<TriggerResults> HLTResults;
	iEvent.getByLabel(triggerCollection_, HLTResults);
	const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTResults);
	bool flag=false;

	if (HLTResults.isValid() && doTheHLTAnalysis_) {
		/// Storing the Prescale information: loop over the triggers and record prescale

		// Matching the HLT information event per event if no hlt info is present in iRun
		if(hltispresent2==false){
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
				  //If path is accepted, then together with its prescale it is stored in a map.
						int prescaleset = hltConfig_.prescaleSet(iEvent,iSetup);
						if(prescaleset!=-1) {
							prescalepair = hltConfig_.prescaleValues(iEvent,iSetup,triggerNames_[itrig]);
							if (Debug2) cout<<"prescale.first "<<prescalepair.first<<" prescalepair.second "<<prescalepair.second<<endl;
							//getting prescale info
							prescale = useCombinedPrescales_ ? prescalepair.first*prescalepair.second : prescalepair.second;
							if((useCombinedPrescales_ && prescalepair.first<0) || prescalepair.second<0) {
								edm::LogWarning("MyAnalyzer") << " Unable to get prescale from event for trigger " << triggerNames.triggerName(itrig) << " :" << prescalepair.first << ", " << prescalepair.second;
								prescale = -999;
							}

							if(prescalepair.first<0 || prescalepair.second<0) { prescale = -999; }
						}

						minimalPrescale = minimalPrescale <  prescale ? minimalPrescale : prescale;
						if (Debug2) cout<<"prescale "<<prescale<<" minimal Prescale "<<minimalPrescale<<" for trigger "<<triggerNames.triggerName(itrig)<<endl;


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


		if (!flag) 
		{
			if(!useAllTriggers_) return false;

		}

	} //chiusura HLT studies

	
  //===========================




  //////////////////////////////////////////////////////////////////////////////////////////////////
  //------ NEW
  //Get Pattuple... Yes, we gave up....
  Handle < pat::ElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ()) return false;


  pat::ElectronCollection::const_iterator highestptele;
  pat::ElectronCollection::const_iterator secondptele;

  int i=0;
  if (electronCollection->size()<=1) return false;
  bool protection=false;


  protection=false;
  /// NEW DS
  // Cutting on WP80
  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    if ( SelectionUtils::DoWP80(recoElectron,iEvent) && SelectionUtils::DoHLTMatch(recoElectron,iEvent) && recoElectron->pt()>20.0){
      if (Debug2) cout<<"Tag is a WP80 electron..."<<endl;
      //Sort in Pt
      if (Debug2) cout<<"Electron pt value ->"<<recoElectron->pt()<<endl;
      if (Debug2) cout<<" MMMM ele trigger size "<<recoElectron->triggerObjectMatches().size()<<endl;
      if (i==0) {
	highestptele=recoElectron;
	i++;
      }
      if (i>0) {
	if (highestptele->pt()<recoElectron->pt()) {
	  highestptele=recoElectron;
	  i++;
	}
      }
    } else {
      if (Debug2) cout<<"Tag IS NOT a WP80 electron...Exit"<<endl;
    }
  }
  
  if (!protection) {
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }

  int j=0;
  protection=false;
  /// NEW DS
  // Cutting on WP80
  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    if ( SelectionUtils::DoWP80(recoElectron,iEvent) && SelectionUtils::DoHLTMatch(recoElectron,iEvent) && recoElectron->pt()>10.0){
      if (Debug2) cout<<"Probe is a WP80 electron..."<<endl;
      //Sort in Pt
      if (Debug2) cout<<"Electron pt value ->"<<recoElectron->pt()<<endl;
      if (Debug2) cout<<" MMMM ele trigger size "<<recoElectron->triggerObjectMatches().size()<<endl;
      if (j==0 && highestptele->pt()>recoElectron->pt()) {
	secondptele=recoElectron;
	j++;
      }
      if (j>0) {
	if (secondptele->pt()<recoElectron->pt() && highestptele->pt()>recoElectron->pt()) {
	  secondptele=recoElectron;
	  j++;
	}
      }
    } else {
      if (Debug2) cout<<"Probe IS NOT a WP80 electron...Exit"<<endl;
    }
  }

  if (!protection) {
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }

  if(i<1 || j<1) return false; //you NEED at least two electrons :)

  //--------------
  // Match the HLT
  pat::ElectronCollection::const_iterator tag;
  pat::ElectronCollection::const_iterator probe;
	  
	  tag=highestptele; //Ã¨ solo un rinominare le cose... (non e' necessario solo retaggio di codice copiato)
	  probe=secondptele;

  //Calculating Invariant Mass
  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(tag->pt(),tag->eta(),tag->phi(), 0.0);
  TLorentzVector probev;
  probev.SetPtEtaPhiM(probe->pt(),probe->eta(),probe->phi(), 0.0);

  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();


  //Cut on the tag and probe mass...
  if (e_ee_invMass>120 || e_ee_invMass<60) return false;

  //Filling Histograms
  h_invMass->Fill(e_ee_invMass);

  //Where is probe in eta
  if ((fabs(probe->eta()) <=1.44) && (fabs(tag->eta()) <=1.44)) h_invMassBB->Fill(e_ee_invMass);
  if ((fabs(probe->eta()) >1.44) && (fabs(tag->eta()) >1.44)) h_invMassEE->Fill(e_ee_invMass);
  if (((fabs(probe->eta()) <=1.44) && (fabs(tag->eta()) >1.44)) || ((fabs(probe->eta()) >1.44) && (fabs(tag->eta()) <=1.44))) h_invMassEB->Fill(e_ee_invMass);

  
  eventAccept->Fill(i);
  //------ END NEW
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////


  return true;

}


// ------------ method called once each job just before starting event loop  ------------
void
ZpatFilter::beginJob (){

  //beginJob

}


// ------------ method called when starting to processes a run  ------------
bool
ZpatFilter::beginRun(edm::Run &iRun, edm::EventSetup const& iSetup)
{

	hltispresent2=true;
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
		hltispresent2=false;
	}

	if (Debug2) cout<<"useAllTriggers? "<<useAllTriggers_<<endl;
	if(useAllTriggers_) triggerNames_ = hlNames;
	//triggerNames_ = hlNames;
	//HLT indices
	triggerIndices_.clear();



	for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
		if(find(hlNames.begin(),hlNames.end(),triggerNames_[itrig])!=hlNames.end()){
			triggerIndices_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
		}
		else{
			triggerIndices_.push_back(2048);
		}
	}
	return true;

	//qui sarebbe piÃ¹ intelligente tornare false se nella lista dei path hlt non c'Ã¨ quello che ci interessa cosÃ¬ eviti di andare comunque ogni evento in cerca di lui (ovvio devi verificare che doTheHLTanalysis sia true qui!)
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZpatFilter::endJob ()
{

  //endJob

}

DEFINE_FWK_MODULE (ZpatFilter);
