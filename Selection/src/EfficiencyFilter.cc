// -*- C++ -*-
//
// Package:    Efficiency
// Class:      Efficiency
// 
/**\class Efficiency Efficiency.cc Zmonitoring/Efficiency/src/Efficiency.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Vieri Candelise & Matteo Marone
//         Created:  Wed May 11 14:53:26 CESDo2011
// $Id: EfficiencyFilter.cc,v 1.17 2012/03/05 14:20:48 schizzi Exp $



// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/EfficiencyFilter.h"
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
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "HLTrigger/HLTcore/interface/TriggerSummaryAnalyzerAOD.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


using namespace std;
using namespace edm;
using namespace reco;

bool Debug=false; //Activate with true if you wonna have verbosity for Debug

//
// member functions
//

// ------------ method called for each event  ------------
bool
EfficiencyFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  if (Debug) cout<<"------- NEW Event -----"<<endl;
  using namespace edm;

  Handle<reco::SuperClusterCollection> superClusters_EB_h;
  iEvent.getByLabel(superClusterCollection_EB_,superClusters_EB_h );
  if ( ! superClusters_EB_h.isValid() ) return false;

  Handle<reco::SuperClusterCollection> superClusters_EE_h;
  iEvent.getByLabel(superClusterCollection_EE_,superClusters_EE_h );
  if ( ! superClusters_EE_h.isValid() ) return false;
  
  reco::SuperClusterCollection::const_iterator highestenergy_SC;
  bool highestenergy_SC_isPassingProbe=false;
  
  int l=0;
  //  int sc_counter=0;
  
  //EB superclusters:
  for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EB_h->begin(); superCluster != superClusters_EB_h->end(); superCluster++) {
    if (superCluster->energy()>20.0 && fabs(superCluster->eta())<=1.4442) {
      //      sc_counter++;
      if (l==0) highestenergy_SC=superCluster;
      if (l>0){
	if (highestenergy_SC->energy() < superCluster->energy()){
	  highestenergy_SC=superCluster;
	}
      }
      l++;
    }
  }
  //EE superclusters:  
  for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EE_h->begin(); superCluster != superClusters_EE_h->end(); superCluster++) {
    if (superCluster->energy()>20.0 && (fabs(superCluster->eta())>=1.5660 && fabs(superCluster->eta())<=2.5000)) {
      //  sc_counter++;
      if (l==0) highestenergy_SC=superCluster;
      if (l>0){
	if (highestenergy_SC->energy() < superCluster->energy()){
	  highestenergy_SC=superCluster;
	}
      }
      l++;
    }
  }

  if (RECO_efficiency_ && l<1) return false;

  //  scNumber_per_event->Fill(sc_counter);

  //Get Pattuple... Yes, we gave up....
  Handle < pat::ElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ()) return false;


  pat::ElectronCollection::const_iterator highestptele;
  pat::ElectronCollection::const_iterator secondptele;

  //  int HLTmatches_PASS=0;
  //  int HLTmatches_TOTAL=0;
  //  int HLTmatches_FAIL=0;

  int i=0;
  if (electronCollection->size()<=1) return false;
  bool protection=false;
  int jj=0;
  //  int n_scEle_match=0;
  int sizePat=electronCollection->size();

  
  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    if (Debug) cout<<"Electron pt value ->"<<recoElectron->pt()<<endl;
    if (Debug) cout<<" MMMM ele trigger size "<<recoElectron->triggerObjectMatches().size()<<endl;
    //    HLTmatches_TOTAL++;
    //    if ( SelectionUtils::DoHLTMatch(recoElectron,iEvent) ) {HLTmatches_PASS++;} else {HLTmatches_FAIL++;}
    if (recoElectron->pt()>10.0 && ((WP80_efficiency_ && SelectionUtils::DoHLTMatch(recoElectron,iEvent)) || HLTele17_efficiency_ || HLTele8NOTele17_efficiency_)) {
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
    }
    if (recoElectron->pt()>10.0 && RECO_efficiency_) {
      if (sqrt((highestenergy_SC->eta()-recoElectron->eta())*(highestenergy_SC->eta()-recoElectron->eta())+(highestenergy_SC->phi()-recoElectron->phi())*(highestenergy_SC->phi()-recoElectron->phi())) < 0.3) {
	highestenergy_SC_isPassingProbe=true;
	if (Debug) cout<<"Ele and SC are matched!"<<endl;
	//	n_scEle_match++;
	continue;
      }
      if (i==0) highestptele=recoElectron;
      if (i>0){
	if (highestptele->pt()<recoElectron->pt()){
	  highestptele=recoElectron;
	}
      }
      i++;
    }
  }
  if (!protection) {
    cout<<"size pat is "<<sizePat<<" while jj is "<<jj<<" and protection "<<protection<<endl;
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }
  
  //  eleNumber_scMatch->Fill(n_scEle_match);

  if (Debug) cout<<"Out of the SC and RECOele loops!"<<endl;  

  //  HLTnumberOfMatches_PASS->Fill(HLTmatches_PASS);
  //  HLTnumberOfMatches_FAIL->Fill(HLTmatches_FAIL);
  //  HLTnumberOfMatches_TOTALELE->Fill(HLTmatches_TOTAL);

  if((WP80_efficiency_ || HLTele17_efficiency_ || HLTele8NOTele17_efficiency_) && (i<2 || highestptele->pt()<20.0)) {return false;}
  if(RECO_efficiency_ && (i<1 || highestptele->pt()<20.0)) {if (Debug) cout << "No TP ele-SC for RECO. i = " << i << endl; return false;}

  //Calculating tag & probe stuff
  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(highestptele->pt(),highestptele->eta(),highestptele->phi(), 0.0);
  TLorentzVector probev;
  if (RECO_efficiency_)  {
    probev.SetXYZT(highestenergy_SC->energy()*cos(highestenergy_SC->phi())/cosh(highestenergy_SC->eta()),
		   highestenergy_SC->energy()*sin(highestenergy_SC->phi())/cosh(highestenergy_SC->eta()),
		   highestenergy_SC->energy()*tanh(highestenergy_SC->eta()),
		   highestenergy_SC->energy());
  } else {
    probev.SetPtEtaPhiM(secondptele->pt(),secondptele->eta(),secondptele->phi(), 0.0);
  }


  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();

  //cut on the tag and probe mass...
  if (e_ee_invMass>120 || e_ee_invMass<60) return false;


  //Estraggo il contenuto dei jets
  int nJet = 0;
  int jetIndex = 0;
  double deltaRCone=0.3;
  
  Handle<PFJetCollection> pfJets;
  iEvent.getByLabel(theJetCollectionLabel_, pfJets);
  if (pfJets.isValid()) {
    PFJetCollection::const_iterator jet = pfJets->begin ();
    
    for (; jet != pfJets->end (); jet++, jetIndex++) {
      // check if the jet is equal to one of the isolated electrons
      double deltaR1;
      if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8NOTele17_efficiency_) {
	deltaR1 = sqrt( pow(jet->eta()-secondptele->eta(),2)+pow(jet->phi()-secondptele->phi(),2) );
      } else {
	deltaR1 = sqrt( pow(jet->eta()-highestenergy_SC->eta(),2)+pow(jet->phi()-highestenergy_SC->phi(),2) );
      }
      double deltaR2= sqrt( pow(jet->eta()-highestptele->eta(),2)+pow(jet->phi()-highestptele->phi(),2) );
      if (deltaR1 > deltaRCone && deltaR2 > deltaRCone && fabs(jet->eta())<2.4 && jet->pt()>30) {
	nJet++;
	if (Debug) cout<<"Jet eta "<<jet->eta()<<" pt "<<jet->pt()<<endl;
      }
    }
  }
  else{
    cout<<"No valid Jets Collection"<<endl;
  }
  
  if (Debug) cout<<"This event has jets #->"<<nJet<<endl;
  

  // Retrieve Number of vertexes
  edm::Handle<reco::VertexCollection> Vertexes;
  iEvent.getByLabel(VertexCollectionTag_, Vertexes);
  int numberOfVertices = Vertexes->size();
  
  
  // Filling TAP distributions:

  // WP80:

  if (WP80_efficiency_ && !New_HE_) {
    // 1st leg WP80 efficiency
    if (SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_)){
      probeall_pt->Fill(secondptele->pt());
      probeall_eta->Fill(secondptele->eta());
      probeall_mee->Fill(e_ee_invMass);
      if ( SelectionUtils::DoWP80Pf(secondptele,iEvent,removePU_) ){
	probepass_pt->Fill(secondptele->pt());
	probepass_eta->Fill(secondptele->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (nJet==0) probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) probepass4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  probepassHighPU->Fill(e_ee_invMass);
	} else {
	  probepassLowPU->Fill(e_ee_invMass);
	}
      }
      else{
	probefail_pt->Fill(secondptele->pt());
	probefail_eta->Fill(secondptele->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (nJet==0) probefail0jet->Fill(e_ee_invMass);
	if (nJet==1) probefail1jet->Fill(e_ee_invMass);
	if (nJet==2) probefail2jet->Fill(e_ee_invMass);
	if (nJet==3) probefail3jet->Fill(e_ee_invMass);
	if (nJet==4) probefail4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  probefailHighPU->Fill(e_ee_invMass);
	} else {
	  probefailLowPU->Fill(e_ee_invMass);
	}
      }
    }  
    if (SelectionUtils::DoWP80Pf(secondptele,iEvent,removePU_)){
      tagall_pt->Fill(highestptele->pt());
      tagall_eta->Fill(highestptele->eta());
      tagall_mee->Fill(e_ee_invMass);
      if ( SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_) ){
	tagpass_pt->Fill(highestptele->pt());
	tagpass_eta->Fill(highestptele->eta());
	tagpass_mee->Fill(e_ee_invMass);
	if (nJet==0) tagpass0jet->Fill(e_ee_invMass);
	if (nJet==1) tagpass1jet->Fill(e_ee_invMass);
	if (nJet==2) tagpass2jet->Fill(e_ee_invMass);
	if (nJet==3) tagpass3jet->Fill(e_ee_invMass);
	if (nJet==4) tagpass4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  tagpassHighPU->Fill(e_ee_invMass);
	} else {
	  tagpassLowPU->Fill(e_ee_invMass);
	}
      }
      else{
	tagfail_pt->Fill(highestptele->pt());
	tagfail_eta->Fill(highestptele->eta());
	tagfail_mee->Fill(e_ee_invMass);
	if (nJet==0) tagfail0jet->Fill(e_ee_invMass);
	if (nJet==1) tagfail1jet->Fill(e_ee_invMass);
	if (nJet==2) tagfail2jet->Fill(e_ee_invMass);
	if (nJet==3) tagfail3jet->Fill(e_ee_invMass);
	if (nJet==4) tagfail4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  tagfailHighPU->Fill(e_ee_invMass);
	} else {
	  tagfailLowPU->Fill(e_ee_invMass);
	}
      }
    }
  }


  // WP80 (New HE!!!):

  if (WP80_efficiency_ && New_HE_) {
    // 1st leg WP80 efficiency
    if (SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_)){
      probeall_pt->Fill(secondptele->pt());
      probeall_eta->Fill(secondptele->eta());
      probeall_mee->Fill(e_ee_invMass);
      if ( SelectionUtils::DoWP80Pf_NewHE(secondptele,iEvent,removePU_) ){
	probepass_pt->Fill(secondptele->pt());
	probepass_eta->Fill(secondptele->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (nJet==0) probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) probepass4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  probepassHighPU->Fill(e_ee_invMass);
	} else {
	  probepassLowPU->Fill(e_ee_invMass);
	}
      }
      else{
	probefail_pt->Fill(secondptele->pt());
	probefail_eta->Fill(secondptele->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (nJet==0) probefail0jet->Fill(e_ee_invMass);
	if (nJet==1) probefail1jet->Fill(e_ee_invMass);
	if (nJet==2) probefail2jet->Fill(e_ee_invMass);
	if (nJet==3) probefail3jet->Fill(e_ee_invMass);
	if (nJet==4) probefail4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  probefailHighPU->Fill(e_ee_invMass);
	} else {
	  probefailLowPU->Fill(e_ee_invMass);
	}
      }
    }  
    if (SelectionUtils::DoWP80Pf_NewHE(secondptele,iEvent,removePU_)){
      tagall_pt->Fill(highestptele->pt());
      tagall_eta->Fill(highestptele->eta());
      tagall_mee->Fill(e_ee_invMass);
      if ( SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_) ){
	tagpass_pt->Fill(highestptele->pt());
	tagpass_eta->Fill(highestptele->eta());
	tagpass_mee->Fill(e_ee_invMass);
	if (nJet==0) tagpass0jet->Fill(e_ee_invMass);
	if (nJet==1) tagpass1jet->Fill(e_ee_invMass);
	if (nJet==2) tagpass2jet->Fill(e_ee_invMass);
	if (nJet==3) tagpass3jet->Fill(e_ee_invMass);
	if (nJet==4) tagpass4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  tagpassHighPU->Fill(e_ee_invMass);
	} else {
	  tagpassLowPU->Fill(e_ee_invMass);
	}
      }
      else{
	tagfail_pt->Fill(highestptele->pt());
	tagfail_eta->Fill(highestptele->eta());
	tagfail_mee->Fill(e_ee_invMass);
	if (nJet==0) tagfail0jet->Fill(e_ee_invMass);
	if (nJet==1) tagfail1jet->Fill(e_ee_invMass);
	if (nJet==2) tagfail2jet->Fill(e_ee_invMass);
	if (nJet==3) tagfail3jet->Fill(e_ee_invMass);
	if (nJet==4) tagfail4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  tagfailHighPU->Fill(e_ee_invMass);
	} else {
	  tagfailLowPU->Fill(e_ee_invMass);
	}
      }
    }
  }


  // HLT:

  if (HLTele17_efficiency_ || HLTele8NOTele17_efficiency_) {
    if ( SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_) ){
      probeall_pt->Fill(secondptele->pt());
      probeall_eta->Fill(secondptele->eta());
      probeall_mee->Fill(e_ee_invMass);
      if ( SelectionUtils::DoHLTMatch(secondptele,iEvent) ){
	probepass_pt->Fill(secondptele->pt());
	probepass_eta->Fill(secondptele->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (nJet==0) probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) probepass4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  probepassHighPU->Fill(e_ee_invMass);
	} else {
	  probepassLowPU->Fill(e_ee_invMass);
	}
      }
      else{
	probefail_pt->Fill(secondptele->pt());
	probefail_eta->Fill(secondptele->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (nJet==0) probefail0jet->Fill(e_ee_invMass);
	if (nJet==1) probefail1jet->Fill(e_ee_invMass);
	if (nJet==2) probefail2jet->Fill(e_ee_invMass);
	if (nJet==3) probefail3jet->Fill(e_ee_invMass);
	if (nJet==4) probefail4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  probefailHighPU->Fill(e_ee_invMass);
	} else {
	  probefailLowPU->Fill(e_ee_invMass);
	}
      }
    }
    if ( SelectionUtils::DoWP80Pf(secondptele,iEvent,removePU_) ){    
      if ( SelectionUtils::DoHLTMatch(highestptele,iEvent) ){
	tagpass_pt->Fill(highestptele->pt());
	tagpass_eta->Fill(highestptele->eta());
	tagpass_mee->Fill(e_ee_invMass);
	if (nJet==0) tagpass0jet->Fill(e_ee_invMass);
	if (nJet==1) tagpass1jet->Fill(e_ee_invMass);
	if (nJet==2) tagpass2jet->Fill(e_ee_invMass);
	if (nJet==3) tagpass3jet->Fill(e_ee_invMass);
	if (nJet==4) tagpass4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  tagpassHighPU->Fill(e_ee_invMass);
	} else {
	  tagpassLowPU->Fill(e_ee_invMass);
	}
      }
      else{
	tagfail_pt->Fill(highestptele->pt());
	tagfail_eta->Fill(highestptele->eta());
	tagfail_mee->Fill(e_ee_invMass);
	if (nJet==0) tagfail0jet->Fill(e_ee_invMass);
	if (nJet==1) tagfail1jet->Fill(e_ee_invMass);
	if (nJet==2) tagfail2jet->Fill(e_ee_invMass);
	if (nJet==3) tagfail3jet->Fill(e_ee_invMass);
	if (nJet==4) tagfail4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  tagfailHighPU->Fill(e_ee_invMass);
	} else {
	  tagfailLowPU->Fill(e_ee_invMass);
	}
      }
    }
  }

  // RECO:

  if (RECO_efficiency_) {
    if (SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_)){
      probeall_pt->Fill(highestenergy_SC->energy());
      probeall_eta->Fill(highestenergy_SC->eta());
      probeall_mee->Fill(e_ee_invMass);
      if (highestenergy_SC_isPassingProbe){
	probepass_pt->Fill(highestenergy_SC->energy());
	probepass_eta->Fill(highestenergy_SC->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (nJet==0) probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) probepass4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  probepassHighPU->Fill(e_ee_invMass);
	} else {
	  probepassLowPU->Fill(e_ee_invMass);
	}
      }
      else{
	probefail_pt->Fill(highestenergy_SC->energy());
	probefail_eta->Fill(highestenergy_SC->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (nJet==0) probefail0jet->Fill(e_ee_invMass);
	if (nJet==1) probefail1jet->Fill(e_ee_invMass);
	if (nJet==2) probefail2jet->Fill(e_ee_invMass);
	if (nJet==3) probefail3jet->Fill(e_ee_invMass);
	if (nJet==4) probefail4jet->Fill(e_ee_invMass);
	if (numberOfVertices > 10) {
	  probefailHighPU->Fill(e_ee_invMass);
	} else {
	  probefailLowPU->Fill(e_ee_invMass);
	}
      }
    }  
  }
  
  return true;
  
}

  
// ------------ method called once each job just before starting event loop  ------------
void
EfficiencyFilter::beginJob (){

}

// ------------ method called when starting to processes a run  ------------
bool
EfficiencyFilter::beginRun(edm::Run &iRun, edm::EventSetup const& iSetup)
{
//HLT names
   std::vector<std::string>  hlNames;
   bool changed (true);
   if (hltConfig_.init(iRun,iSetup,triggerCollection_.process(),changed)) {
     if (changed) {
       hlNames = hltConfig_.triggerNames();
     }
   } else {
     edm::LogError("MyAnalyzer") << " HLT config extraction failure with process name " << triggerCollection_.process();
   }
   if (Debug) cout<<"useAllTriggers?"<<useAllTriggers_<<endl;
   if(useAllTriggers_) triggerNames_ = hlNames;
   //triggerNames_ = hlNames;
   //HLT indices
   triggerIndices_.clear();
   for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
     if(find(hlNames.begin(),hlNames.end(),triggerNames_[itrig])!=hlNames.end())
       triggerIndices_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
     else
       triggerIndices_.push_back(2048);
   }

   if (Debug){
     // text (Debug) output
     int i=0;
     for(std::vector<std::string>::const_iterator it = triggerNames_.begin(); it<triggerNames_.end();++it) {
       std::cout << (i++) << " = " << (*it) << std::endl;
     } 
   }
   return true;
}



// ------------ method called once each job just after ending the event loop  ------------
void
EfficiencyFilter::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (EfficiencyFilter);
