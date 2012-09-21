// -*- C++ -*-
//
// Package:    Zanalyzer
// Class:      Zanalyzer
// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/ZpatFilterMu2011.h"
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

bool DebugMu11=false;    //Activate with true if you wonna have verbosity for Debug

  using namespace edm;

//
// member functions
//

// ------------ method called for each event  ------------
bool
ZpatFilterMu2011::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{

  muSelStepByStep->SetBinContent(1,muSelStepByStep->GetBinContent(1)+1); //Total Number of events, HLT fired + Selection (1)  

  if (DebugMu11) cout<<"------- NEW Event -----"<<endl;
   
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //------ NEW
  //Get Pattuple... Yes, we gave up....
  Handle < pat::MuonCollection > muonCollection;
  iEvent.getByLabel (theMuCollectionLabel, muonCollection);
  if (!muonCollection.isValid ())  return false;
  
  pat::MuonCollection::const_iterator highestptmu;
  pat::MuonCollection::const_iterator secondptmu;
  
  double lowZmassLimit=lowZmassLimit_;
  double highZmassLimit=highZmassLimit_;

  int i=0;

  if (muonCollection->size()<=1) return false;
 
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

  if(i<2) return false; //you NEED at least two electrons :)
  
  //Check if the charge are opposite..
  if(highestptmu->charge() == secondptmu->charge()) return false;

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
  if (e_ee_invMass>highZmassLimit || e_ee_invMass<lowZmassLimit) return false;

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
ZpatFilterMu2011::beginJob (){
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
ZpatFilterMu2011::endJob ()
{

  //endJob

}

DEFINE_FWK_MODULE (ZpatFilterMu2011);
