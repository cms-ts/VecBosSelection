// -*- C++ -*-

// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/analyzerMuCuts.h"
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

//
// member functions
//

// ------------ method called for each event  ------------
void
analyzerMuCuts::analyze (const edm::Event & iEvent, const edm::EventSetup & iSetup)
{
 
  muSelStepByStep->SetBinContent(1,muSelStepByStep->GetBinContent(1)+1); //Total Number of events

  Handle < pat::MuonCollection > muonCollection;
  iEvent.getByLabel (theMuCollectionLabel, muonCollection);

  Handle < pat::MuonCollection > muonMatchedCollection;
  iEvent.getByLabel (theMuMatchedCollectionLabel, muonMatchedCollection);

  double highZmassLimit = 111;
  double lowZmassLimit = 71;

  //a bunch of counters...
  int inAcceptance=0;
  int inPt=0;
  int isGlobalAndTracker=0;
  int isDB=0;
  int isChi2=0;
  int inHLT=0;
  int isPixelAndMuonAndTracker=0;

  if (muonCollection.isValid ()) {
    // 3-> 2 muons
    if (muonCollection->size()<=1) {
      return;
    }
    muSelStepByStep->SetBinContent(2,muSelStepByStep->GetBinContent(2)+1); //AT least 1 muon
    
    for (pat::MuonCollection::const_iterator recoMu = muonCollection->begin (); recoMu != muonCollection->end (); recoMu++) {
      
      if (recoMu->triggerObjectMatches().size()>0) inHLT++;
      if (fabs(recoMu->eta())<=2.4) inAcceptance++;
      bool isGlobal=recoMu->isGlobalMuon();
      bool isTracker=recoMu->isTrackerMuon();
      if (isGlobal==true && isTracker==true) isGlobalAndTracker++;
      double muonpt=recoMu->pt();
      if (muonpt>20) inPt++;
      //7 -> is DB<0.02
      if (recoMu->dB()<0.02) isDB++;
      if (isGlobal) {
	if (recoMu->globalTrack()->normalizedChi2()<10) isChi2++;
      }
      if (isGlobal==true && isTracker==true){
	if (recoMu->innerTrack()->hitPattern().numberOfValidPixelHits()>0 
	    && recoMu->globalTrack()->hitPattern().numberOfValidMuonHits()>0 
	    && recoMu->innerTrack()->hitPattern().trackerLayersWithMeasurement()>8) isPixelAndMuonAndTracker++;
      }
    }
    
    // 3-> 2 HLT machted mu    
    if (inHLT>1) {
      muSelStepByStep->SetBinContent(3,muSelStepByStep->GetBinContent(3)+1); //In acceptance
    }
    else {return;}
    
    // 4-> 2 within acceptance    
    if (inAcceptance>1) {
      muSelStepByStep->SetBinContent(4,muSelStepByStep->GetBinContent(4)+1); //In acceptance
    }
    else {return;}
    
    if (isGlobalAndTracker>1) {
      muSelStepByStep->SetBinContent(5,muSelStepByStep->GetBinContent(5)+1); //Global and Tracker
    }
    else{return;}
    
    if (isDB>1 && isChi2>1){
      muSelStepByStep->SetBinContent(6,muSelStepByStep->GetBinContent(6)+1); // DB & Chi2
    }
    else{return;}

    if (inPt>1){
      muSelStepByStep->SetBinContent(7,muSelStepByStep->GetBinContent(7)+1); // Pt>20
    }
    else{return;}

    if (isPixelAndMuonAndTracker>1){
      muSelStepByStep->SetBinContent(8,muSelStepByStep->GetBinContent(8)+1); // Pt>20
    }
    else{return;}
  }
  else{
    cout<<"Muon Collection in analyzerMuCuts is not Valid"<<endl;
    return;
  }
  
  //Importare MacthedMuons, e replicare codice controllo Carica e Massa Invariante
  if (!muonMatchedCollection.isValid ()){
     cout<<"Muon Matched Collection in analyzerMuCuts is not Valid"<<endl;
     return;
  }

  if ( muonMatchedCollection->size()<2) return;
  muSelStepByStep->SetBinContent(9,muSelStepByStep->GetBinContent(9)+1); // 2 mu in Selezioni matched

  // search the two highest pt muons
  pat::MuonCollection::const_iterator highestptmu;
  pat::MuonCollection::const_iterator secondptmu;
  int i=0;
  for (pat::MuonCollection::const_iterator itmu = muonMatchedCollection->begin (); 
       itmu != muonMatchedCollection->end (); itmu++) {

    if (i==0) highestptmu=itmu;
    if (i==1){
       if (highestptmu->pt()<itmu->pt()){
          secondptmu=highestptmu;
          highestptmu=itmu;
       }
       else{
          secondptmu=itmu;
       }
    }
    if (i>1){
       if (highestptmu->pt()<itmu->pt()){
          secondptmu=highestptmu;
          highestptmu=itmu;
       }
       else{
          if (secondptmu->pt()<itmu->pt()) secondptmu=itmu;
       }
    }
    i++;
  }
  //Check if the charge are opposite
  if(highestptmu->charge() == secondptmu->charge()) return;  
  muSelStepByStep->SetBinContent(10,muSelStepByStep->GetBinContent(10)+1); // 2 mu with opposite charge 

  //Calculating Invariant Mass
  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(highestptmu->pt(),highestptmu->eta(),highestptmu->phi(), 0.0);
  TLorentzVector probev;
  probev.SetPtEtaPhiM(secondptmu->pt(),secondptmu->eta(),secondptmu->phi(), 0.0);
  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();
  //Cut on the tag and probe mass...
  if (e_ee_invMass>highZmassLimit || e_ee_invMass<lowZmassLimit) return ;
  muSelStepByStep->SetBinContent(11,muSelStepByStep->GetBinContent(11)+1); // Invariant mass [71,111]

  return;
}


// ------------ method called once each job just before starting event loop  ------------
void
analyzerMuCuts::beginJob (){

  muSelStepByStep->GetXaxis()->SetBinLabel(1,"Total # of Events");
  muSelStepByStep->GetXaxis()->SetBinLabel(2,">= 2 muons");
  muSelStepByStep->GetXaxis()->SetBinLabel(3,">= 2 muons HLT fired");
  muSelStepByStep->GetXaxis()->SetBinLabel(4,">=2 muons in acceptance");
  muSelStepByStep->GetXaxis()->SetBinLabel(5,">=2 Global & Tracker");
  muSelStepByStep->GetXaxis()->SetBinLabel(6,">=2 DB<0.02 and Chi2<10");
  muSelStepByStep->GetXaxis()->SetBinLabel(7,">=2 pt>20");
  muSelStepByStep->GetXaxis()->SetBinLabel(8,">=2 Hits in Tracker & Pixel & Muon");
  muSelStepByStep->GetXaxis()->SetBinLabel(9,">=1 Zmumu Candadates");
  muSelStepByStep->GetXaxis()->SetBinLabel(10,"Opposite Charge");
  muSelStepByStep->GetXaxis()->SetBinLabel(11,"InvariantMass");
}


// ------------ method called once each job just after ending the event loop  ------------
void
analyzerMuCuts::endJob ()
{

  //endJob

}

DEFINE_FWK_MODULE (analyzerMuCuts);
