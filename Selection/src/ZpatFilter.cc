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
  eleSelStepByStep->SetBinContent(1,eleSelStepByStep->GetBinContent(1)+1); //Number of events in which hlt has fired (1)
	
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

	
  eleSelStepByStep->SetBinContent(2,eleSelStepByStep->GetBinContent(2)+1); //Number of events in which hlt has fired (1)

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //------ NEW
  //Get Pattuple... Yes, we gave up....
  Handle < pat::ElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ())  return false;

  pat::ElectronCollection::const_iterator highestptele;
  pat::ElectronCollection::const_iterator secondptele;

  double secondEleEnThrhold=secondEleEnThrhold_;
  double firstEleEnThrhold=firstEleEnThrhold_;  
  double lowZmassLimit=lowZmassLimit_;
  double highZmassLimit=highZmassLimit_;
  double maxEtaForElectron=maxEtaForElectron_;

  //Set of counters to follow the Z selection history  and form the plot labels..
  int twoEleGoodEtaCount=0;
  int twoEleHLTCount=0;
  int WP80Count=0;
  int lowThrholdCount=0;
  int isID=0; 
  int isIso=0;
  int isConv=0;
  
  int i=0;

  if (electronCollection->size()<=1) {
    numberOfEleAfterHLTTrigger->SetBinContent(1,numberOfEleAfterHLTTrigger->GetBinContent(1)+1);
    numberOfEleAfterHLTTrigger->SetBinContent(3,numberOfEleAfterHLTTrigger->GetBinContent(3)+1);
    cout<<"Events with DoubleEle trigger but eleCollection size <2!!! suspicious! -> "<<theElectronCollectionLabel<<endl;
    return false;
  }
  numberOfEleAfterHLTTrigger->SetBinContent(1,numberOfEleAfterHLTTrigger->GetBinContent(1)+1);
  numberOfEleAfterHLTTrigger->SetBinContent(2,numberOfEleAfterHLTTrigger->GetBinContent(2)+1);
  bool protection=false;
  int jj=0;
  int sizePat=electronCollection->size();
  eleSelStepByStep->SetBinContent(3,eleSelStepByStep->GetBinContent(3)+1); // (1) + at least 2 ele in the event (2)

  protection=false;
  /// NEW DS
  // Cutting on WP80
  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    jj++;
    //Check whether the electron is within the acceptance
    if (fabs(recoElectron ->superCluster()->eta()) > maxEtaForElectron) {
      continue;
    }
    twoEleGoodEtaCount++;
    if (twoEleGoodEtaCount==2) {
      eleSelStepByStep->SetBinContent(4,eleSelStepByStep->GetBinContent(4)+1); //(2) + 2 ele whithin eta acceptance (3)
    }
    if (SelectionUtils::DoHLTMatch(recoElectron,iEvent)) twoEleHLTCount++;
    if (twoEleHLTCount==2) {
      eleSelStepByStep->SetBinContent(5,eleSelStepByStep->GetBinContent(5)+1); //(3) + 2 ele HLT matched (4)
      twoEleHLTCount++; //To avoid multiple insertion...
    }
    //Perform checks on each ele ID criteria
    // Here you get a plot full of information. Each electron contributes with one entry (so total numer of entries = 3* #electrons)
    // To have the "%", each bin value shold be divided by total numer of entries/3
    std::vector<bool> result=SelectionUtils::MakeEleIDAnalysis(recoElectron,iEvent,removePU_); 
    passIDEleCriteria->SetBinContent(1,passIDEleCriteria->GetBinContent(1)+1);
    if (result[0]) {
      passIDEleCriteria->SetBinContent(2,passIDEleCriteria->GetBinContent(2)+1);
      if (isIso<2) isIso++;
    }
    if (result[1]) {
      passIDEleCriteria->SetBinContent(3,passIDEleCriteria->GetBinContent(3)+1);
      if (isID<2) isID++;
    }
    if (result[2]) {
      passIDEleCriteria->SetBinContent(4,passIDEleCriteria->GetBinContent(4)+1);
      if (isConv<2) isConv++;
    }
    if (isID==2) eleSelStepByStep->SetBinContent(6,eleSelStepByStep->GetBinContent(6)+1); // (5) + 2 ele pt > lowTh (6)
    if (isID==2 && isIso==2) eleSelStepByStep->SetBinContent(7,eleSelStepByStep->GetBinContent(7)+1); // (5) + 2 ele pt > lowTh (6)
    if (isID==2 && isIso==2 && isConv==2) eleSelStepByStep->SetBinContent(8,eleSelStepByStep->GetBinContent(8)+1); // (5) + 2 ele pt > lowTh (6)

    if (result[0] && result[1] && result[2]) WP80Count++;
    if (WP80Count==2) eleSelStepByStep->SetBinContent(9,eleSelStepByStep->GetBinContent(9)+1); //(4) + 2 ele WP80 (5)

    if ( SelectionUtils::DoWP80(recoElectron,iEvent,removePU_) && SelectionUtils::DoHLTMatch(recoElectron,iEvent) && recoElectron->pt()>secondEleEnThrhold){
       lowThrholdCount++;
      if (lowThrholdCount==2)eleSelStepByStep->SetBinContent(10,eleSelStepByStep->GetBinContent(10)+1); // (5) + 2 ele pt > lowTh (6)

      if (Debug2) cout<<"Tag is a WP80 electron..."<<endl;
      //Sort in Pt
      if (Debug2) cout<<"Electron pt value ->"<<recoElectron->pt()<<endl;
      if (Debug2) cout<<" MMMM ele trigger size "<<recoElectron->triggerObjectMatches().size()<<endl;
      if (i==0) highestptele=recoElectron;
      if (i==1){
	if (highestptele->pt()<recoElectron->pt()){
	  secondptele=highestptele;
	  highestptele=recoElectron;
	}
	else{
	  secondptele=recoElectron;
	}
      }
      if (i>1){
	if (highestptele->pt()<recoElectron->pt()){
	  secondptele=highestptele;
	  highestptele=recoElectron;
	}
	else{
	  if (secondptele->pt()<recoElectron->pt()) secondptele=recoElectron;
	}
      }
      i++;
    } else {
      if (Debug2) cout<<"Tag IS NOT a WP80 electron...Exit"<<endl;
    }
  }

  if (!protection) {
    cout<<"size pat is "<<sizePat<<" while jj is "<<jj<<" and protection "<<protection<<endl;
    cout<<"problems with PAT collection, in ZpatFilter.cc-->... Please check..."<<endl;    
    return false;
  }

  if(i<2 || highestptele->pt()<firstEleEnThrhold) return false; //you NEED at least two electrons :)
  eleSelStepByStep->SetBinContent(11,eleSelStepByStep->GetBinContent(11)+1); // (8) + 1 ele pt > LowPt and + 1 ele pt > highPt (7)

  //Check if the charge are opposite..
  if(highestptele->charge() == secondptele->charge()) return false;
   eleSelStepByStep->SetBinContent(12,eleSelStepByStep->GetBinContent(12)+1); // (9) + Opposite Charge

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
  if (e_ee_invMass>highZmassLimit || e_ee_invMass<lowZmassLimit) return false;
  eleSelStepByStep->SetBinContent(13,eleSelStepByStep->GetBinContent(13)+1);

  //Filling Histograms
  h_invMass->Fill(e_ee_invMass);

  //Where is probe in eta
  if ((fabs(probe->eta()) <=1.44) && (fabs(tag->eta()) <=1.44)) h_invMassBB->Fill(e_ee_invMass);
  if ((fabs(probe->eta()) >1.44) && (fabs(tag->eta()) >1.44)) h_invMassEE->Fill(e_ee_invMass);
  if (((fabs(probe->eta()) <=1.44) && (fabs(tag->eta()) >1.44)) || ((fabs(probe->eta()) >1.44) && (fabs(tag->eta()) <=1.44))) h_invMassEB->Fill(e_ee_invMass);

  //------ END NEW
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////


  return true;

}


// ------------ method called once each job just before starting event loop  ------------
void
ZpatFilter::beginJob (){
  cout<<endl;
  cout<<"##############################"<<endl;
  cout<<"#   Z Selection Parameters   #"<<endl;
  cout<<"##############################"<<endl;
  cout<<endl; 
  cout<<"Transverse Energy cut on first Ele="<<firstEleEnThrhold_<<"GeV, and "<<secondEleEnThrhold_<<"GeV on the second"<<endl;
  cout<<"Z invariant mass limit: low="<<lowZmassLimit_<<"GeV, high="<<highZmassLimit_<<"GeV"<<endl;
  cout<<"Electron max acceptance="<<maxEtaForElectron_<<endl;
  cout<<"List of triggers being included:"<<endl;
  for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
    cout<<triggerNames_[itrig]<<" ";
  }
  cout<<endl;
  eleSelStepByStep->GetXaxis()->SetBinLabel(1,"Total # of Events");
  eleSelStepByStep->GetXaxis()->SetBinLabel(2,"event HLT Fired");
  eleSelStepByStep->GetXaxis()->SetBinLabel(3,">= 2 ele");
  eleSelStepByStep->GetXaxis()->SetBinLabel(4,"2 ele <= eta Acceptance");
  eleSelStepByStep->GetXaxis()->SetBinLabel(5,"2 ele HLT matched");
  eleSelStepByStep->GetXaxis()->SetBinLabel(6,"2 ele ID");
  eleSelStepByStep->GetXaxis()->SetBinLabel(7,"2 ele Isolated");
  eleSelStepByStep->GetXaxis()->SetBinLabel(8,"2 ele !Converted");
  eleSelStepByStep->GetXaxis()->SetBinLabel(9,"2 ele WP80");
  eleSelStepByStep->GetXaxis()->SetBinLabel(10,"2 ele pt > lowPt");
  eleSelStepByStep->GetXaxis()->SetBinLabel(11,"ele pt > lowPt + pt > HighPt");
  eleSelStepByStep->GetXaxis()->SetBinLabel(12,"Opposite Charge");
  eleSelStepByStep->GetXaxis()->SetBinLabel(13,"Within Z window mass");
  passIDEleCriteria->GetXaxis()->SetBinLabel(1,"TotEle");
  passIDEleCriteria->GetXaxis()->SetBinLabel(2,"Isolated");
  passIDEleCriteria->GetXaxis()->SetBinLabel(3,"ID");
  passIDEleCriteria->GetXaxis()->SetBinLabel(4,"NotConverted");
  numberOfEleAfterHLTTrigger->GetXaxis()->SetBinLabel(1,"Total");
  numberOfEleAfterHLTTrigger->GetXaxis()->SetBinLabel(2,"Size >=2");
  numberOfEleAfterHLTTrigger->GetXaxis()->SetBinLabel(3,"Size <2");

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
