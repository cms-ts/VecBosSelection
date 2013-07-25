// -*- C++ -*-
//
// Package:    Zanalyzer
// Class:      Zanalyzer
// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/ZpatFilterMu2011Unfolding.h"
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

bool DebugMu111=false;    //Activate with true if you wonna have verbosity for Debug

  using namespace edm;

//
// member functions
//

// ------------ method called for each event  ------------
bool
ZpatFilterMu2011Unfolding::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{


  bool valueToBeReturnedReco=true;
  bool valueToBeReturnedGen=false;

  //Evaluate if the events would have been selected at the generator level...
  if (isUnfolding_){
    edm::Handle<reco::GenParticleCollection> genPart;
    iEvent.getByLabel (genParticleCollection_,genPart);
    //Find the ele of the Z  
    int leptonId=11;
    TLorentzVector l1, l2, l_pair;
    std::vector<TLorentzVector> leptonContainer;
    double zInvMass = 0;  

    double whichlepton = 13;
    
    //CAMBINARE COLLEZIONE IN MODO CHE LooPPI SUI ELETTRONI CON I GAMMA SOMMATI!!!!!!!!!!!!!

    for(reco::GenParticleCollection::const_iterator itgen=genPart->begin();itgen!=genPart->end();itgen++){
      if ( fabs(itgen->pdgId())==whichlepton && itgen->mother()->pdgId()==23){ //itgen->status()==1 &&
	l1.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());
	leptonContainer.push_back(l1);
	if (DebugMu111) cout<<"lepton id->"<<itgen->pdgId()<<" eta"<<itgen->eta()<<" et->"<<itgen->et()<<" mass->"<<itgen->mass()<<endl;
      } 
    }
    
    if (leptonContainer.size()>=2) {
      l_pair = leptonContainer[0] + leptonContainer[1];
      zInvMass = l_pair.M();
      if (DebugMu111) cout<<" invMass->"<<zInvMass<<endl;  
      
      //Check If it is a valid generated event or not...
      if ( (fabs(leptonContainer[0].Eta())<=2.4 && fabs(leptonContainer[1].Eta())<=2.4) && (fabs(leptonContainer[0].Pt())>=20 && fabs(leptonContainer[1].Pt())>=20) && (zInvMass>=71 && zInvMass<=111) ){
	if (DebugMu111) cout<<"This is a well generated Z Boson that you sohuld have reconstructed..."<<endl;
	valueToBeReturnedGen=true;
      }
      else{
	valueToBeReturnedGen=false;
	if (DebugMu111) cout<<"This is not a Z Boson inside the acceptance"<<endl;
      }
    }
    else{
      valueToBeReturnedGen=false;
    }
    leptonContainer.clear();
  }

if (DebugMu111)  cout<<"What gen level returns->"<<valueToBeReturnedGen<<endl;

  //=============================================================//
  //=============================================================//
  //=============================================================//

  muSelStepByStep->SetBinContent(1,muSelStepByStep->GetBinContent(1)+1); //Total Number of events, 1 HLT fired + Selection (1)  

  if (DebugMu111) cout<<"------- NEW Event -----"<<endl;
   
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //------ NEW
  //Get Pattuple... Yes, we gave up....
  Handle < pat::MuonCollection > muonCollection;
  iEvent.getByLabel (theMuCollectionLabel, muonCollection);
  if (!muonCollection.isValid ())  return valueToBeReturnedGen;
  
  pat::MuonCollection::const_iterator highestptmu;
  pat::MuonCollection::const_iterator secondptmu;
  
  double lowZmassLimit=lowZmassLimit_;
  double highZmassLimit=highZmassLimit_;

  int i=0;

  if (muonCollection->size()<=1) return valueToBeReturnedGen;
 
  bool protection=false;
  int jj=0;
  int sizePat=muonCollection->size();

 /// NEW DS
  // Cutting on WP80
  for (pat::MuonCollection::const_iterator recoElectron = muonCollection->begin (); recoElectron != muonCollection->end (); recoElectron++) {
    protection=true;
    jj++;

    if (i==0) highestptmu=recoElectron;
    if (i==1){
       if (highestptmu->pt()<recoElectron->pt()){
	  secondptmu=highestptmu;
	  highestptmu=recoElectron;
       }
       else{
	  secondptmu=recoElectron;
       }
    }
    if (i>1){
       if (highestptmu->pt()<recoElectron->pt()){
	  secondptmu=highestptmu;
	  highestptmu=recoElectron;
       }
       else{
	  if (secondptmu->pt()<recoElectron->pt()) secondptmu=recoElectron;
       }
    }
    i++;
  }

  if (!protection) {
    cout<<"size pat is "<<sizePat<<" while jj is "<<jj<<" and protection "<<protection<<endl;
    cout<<"problems with PAT collection, in ZpatFilter.cc-->... Please check..."<<endl;    
    return false;
  }

  if(i<2) return valueToBeReturnedGen; //you NEED at least two electrons :)
  
  //Check if the charge are opposite..
  if(highestptmu->charge() == secondptmu->charge()) return valueToBeReturnedGen;

  //--------------
  // Match the HLT
  pat::MuonCollection::const_iterator tag;
  pat::MuonCollection::const_iterator probe;
	  
	  tag=highestptmu; // solo un rinominare le cose... (non e' necessario solo retaggio di codice copiato)
	  probe=secondptmu;

  //Calculating Invariant Mass
  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(tag->pt(),tag->eta(),tag->phi(), 0.0);
  TLorentzVector probev;
  probev.SetPtEtaPhiM(probe->pt(),probe->eta(),probe->phi(), 0.0);

  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();


  //Cut on the tag and probe mass...
  if (e_ee_invMass>highZmassLimit || e_ee_invMass<lowZmassLimit) return valueToBeReturnedGen;

  //Pippo -> Number Of Events having more than 2 electrons and eta < 2.4 & 1 HLT+WP80 with > 10 GeV & 1 HLT+WP80 with > 20 GeV and Withihn the window energy mass
  if (i> 2 ) h_zPt_3mu->Fill(e_pair.Pt());
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
ZpatFilterMu2011Unfolding::beginJob (){
  cout<<endl;
  cout<<"##############################"<<endl;
  cout<<"#   Z Selection Parameters   #"<<endl;
  cout<<"##############################"<<endl;
  cout<<endl; 

  cout<<"Z invariant mass limit: low="<<lowZmassLimit_<<"GeV, high="<<highZmassLimit_<<"GeV"<<endl;
  cout<<endl;

  muSelStepByStep->GetXaxis()->SetBinLabel(1,"Total # of Events");

}


// ------------ method called once each job just after ending the event loop  ------------
void
ZpatFilterMu2011Unfolding::endJob ()
{

  //endJob

}

DEFINE_FWK_MODULE (ZpatFilterMu2011Unfolding);
