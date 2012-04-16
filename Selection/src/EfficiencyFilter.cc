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
// $Id: EfficiencyFilter.cc,v 1.22 2012/03/29 15:07:12 schizzi Exp $



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

  int HLTmatches_PASS=0;
  int HLTmatches_TOTAL=0;
  int HLTmatches_FAIL=0;

  int i=0;
  if (electronCollection->size()<=1) return false;
  bool protection=false;
  
  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    if (Debug) cout<<"Electron pt value ->"<<recoElectron->pt()<<endl;
    if (Debug) cout<<" MMMM ele trigger size "<<recoElectron->triggerObjectMatches().size()<<endl;
    HLTmatches_TOTAL++;
    if ( SelectionUtils::DoHLTMatch(recoElectron,iEvent) ) {HLTmatches_PASS++;} else {HLTmatches_FAIL++;}
    if (recoElectron->pt()>20.0 && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_)) {
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
    if (recoElectron->pt()>20.0 && RECO_efficiency_) {
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
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }
  
  if (Debug) cout<<"Out of the SC and RECOele loops!"<<endl;  

  HLTnumberOfMatches_PASS->Fill(HLTmatches_PASS);
  HLTnumberOfMatches_FAIL->Fill(HLTmatches_FAIL);
  HLTnumberOfMatches_TOTALELE->Fill(HLTmatches_TOTAL);

  if((WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) && i<2) return false;
  if(i<1) return false;

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
  if (e_ee_invMass>111 || e_ee_invMass<71) return false;


  //Estraggo il contenuto dei jets
  int nJet = 0;
  double deltaRCone=0.3;
  double leadingJet_pT=0.0;
  
  Handle<PFJetCollection> pfJets;
  iEvent.getByLabel(theJetCollectionLabel_, pfJets);
  if (pfJets.isValid()) {
    for(reco::PFJetCollection::const_iterator jet = pfJets->begin();jet != pfJets->end();jet++) { 
      // check if the jet is equal to one of the isolated electrons
      double deltaR1;
      if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
	deltaR1 = sqrt( pow(jet->eta()-secondptele->eta(),2)+pow(jet->phi()-secondptele->phi(),2) );
      } else {
	deltaR1 = sqrt( pow(jet->eta()-highestenergy_SC->eta(),2)+pow(jet->phi()-highestenergy_SC->phi(),2) );
      }
      double deltaR2= sqrt( pow(jet->eta()-highestptele->eta(),2)+pow(jet->phi()-highestptele->phi(),2) );
      if (deltaR1 > deltaRCone && deltaR2 > deltaRCone 
	  && fabs(jet->eta())<2.4 
	  && jet->pt()>30
	  && jet->chargedEmEnergyFraction() < 0.99
	  && jet->neutralHadronEnergyFraction() < 0.99
	  && jet->neutralEmEnergyFraction() < 0.99
	  && jet->chargedHadronEnergyFraction() > 0.0
	  && jet->chargedMultiplicity() > 0
	  ) {
	nJet++;
	if (jet->pt() > leadingJet_pT) leadingJet_pT = jet->pt();
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

  // WP80 & HLT:

  if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
    // 1st leg WP80 efficiency
    if ((WP80_efficiency_ && SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_) && SelectionUtils::DoHLTMatch(highestptele,iEvent)) ||
	((HLTele17_efficiency_ || HLTele8_efficiency_) && SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_))
	){
      probeall_pt->Fill(secondptele->pt());
      probeall_eta->Fill(secondptele->eta());
      probeall_mee->Fill(e_ee_invMass);
      probeall_leadjetpt->Fill(leadingJet_pT);
      if ((WP80_efficiency_ && !New_HE_ && SelectionUtils::DoWP80Pf(secondptele,iEvent,removePU_)) || 
	  (WP80_efficiency_ && New_HE_ && SelectionUtils::DoWP80Pf_NewHE(secondptele,iEvent,removePU_)) ||
	  ((HLTele17_efficiency_ || HLTele8_efficiency_) && SelectionUtils::DoHLTMatch(secondptele,iEvent))
	  ){
	probepass_pt->Fill(secondptele->pt());
	probepass_eta->Fill(secondptele->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (nJet==0) probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) probepass4jet->Fill(e_ee_invMass);
	if (nJet>4)  probepass5jet->Fill(e_ee_invMass);
	if (numberOfVertices < 5) probepass0pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 5 && numberOfVertices < 10) probepass1pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 10 && numberOfVertices < 15) probepass2pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 15 && numberOfVertices < 20) probepass3pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 20) probepass4pu->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) < 0.5) probepass0eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 0.5 && fabs(secondptele->eta()) <1.0)  probepass1eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 1.0 && fabs(secondptele->eta()) <1.4) probepass2eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) > 1.6 && fabs(secondptele->eta()) <2.0)  probepass3eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 2.0 && fabs(secondptele->eta()) <3.0)  probepass4eta->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet == 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt0nJet->Fill(e_ee_invMass);
	}
	if (nJet == 1) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt1nJet->Fill(e_ee_invMass);
	}
	if (nJet == 2) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt2nJet->Fill(e_ee_invMass);
	}
	if (nJet == 3) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt3nJet->Fill(e_ee_invMass);
	}
	if (nJet == 4) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt4nJet->Fill(e_ee_invMass);
	}
	if (nJet >= 5) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt5nJet->Fill(e_ee_invMass);
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
	if (nJet>4)  probefail5jet->Fill(e_ee_invMass);
	if (numberOfVertices < 5) probefail0pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 5 && numberOfVertices < 10) probefail1pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 10 && numberOfVertices < 15) probefail2pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 15 && numberOfVertices < 20) probefail3pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 20) probefail4pu->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) < 0.5) probefail0eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 0.5 && fabs(secondptele->eta()) <1.0)  probefail1eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 1.0 && fabs(secondptele->eta()) <1.4) probefail2eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) > 1.6 && fabs(secondptele->eta()) <2.0)  probefail3eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 2.0 && fabs(secondptele->eta()) <3.0)  probefail4eta->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet == 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt0nJet->Fill(e_ee_invMass);
	}
	if (nJet == 1) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt1nJet->Fill(e_ee_invMass);
	}
	if (nJet == 2) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt2nJet->Fill(e_ee_invMass);
	}
	if (nJet == 3) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt3nJet->Fill(e_ee_invMass);
	}
	if (nJet == 4) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt4nJet->Fill(e_ee_invMass);
	}
	if (nJet >= 5) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt5nJet->Fill(e_ee_invMass);
	}
      }
    }  
    if ((WP80_efficiency_ && SelectionUtils::DoWP80Pf(secondptele,iEvent,removePU_) && SelectionUtils::DoHLTMatch(secondptele,iEvent)) ||
	((HLTele17_efficiency_ || HLTele8_efficiency_) && SelectionUtils::DoWP80Pf(secondptele,iEvent,removePU_))
	){
      tagall_pt->Fill(highestptele->pt());
      tagall_eta->Fill(highestptele->eta());
      tagall_mee->Fill(e_ee_invMass);
      tagall_leadjetpt->Fill(leadingJet_pT);
      if ( (WP80_efficiency_ && !New_HE_ && SelectionUtils::DoWP80Pf(highestptele,iEvent,removePU_)) || 
	   (WP80_efficiency_ && New_HE_ && SelectionUtils::DoWP80Pf_NewHE(highestptele,iEvent,removePU_)) ||
	   ((HLTele17_efficiency_ || HLTele8_efficiency_) && SelectionUtils::DoHLTMatch(highestptele,iEvent))
	   ){
	tagpass_pt->Fill(highestptele->pt());
	tagpass_eta->Fill(highestptele->eta());
	tagpass_mee->Fill(e_ee_invMass);
	if (nJet==0) tagpass0jet->Fill(e_ee_invMass);
	if (nJet==1) tagpass1jet->Fill(e_ee_invMass);
	if (nJet==2) tagpass2jet->Fill(e_ee_invMass);
	if (nJet==3) tagpass3jet->Fill(e_ee_invMass);
	if (nJet==4) tagpass4jet->Fill(e_ee_invMass);
	if (nJet>4)  tagpass5jet->Fill(e_ee_invMass);
	if (numberOfVertices < 5) tagpass0pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 5 && numberOfVertices < 10) tagpass1pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 10 && numberOfVertices < 15) tagpass2pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 15 && numberOfVertices < 20) tagpass3pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 20) tagpass4pu->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) < 0.5) tagpass0eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 0.5 && fabs(secondptele->eta()) <1.0)  tagpass1eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 1.0 && fabs(secondptele->eta()) <1.4) tagpass2eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) > 1.6 && fabs(secondptele->eta()) <2.0)  tagpass3eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 2.0 && fabs(secondptele->eta()) <3.0)  tagpass4eta->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagpass0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagpass1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagpass2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagpass3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagpass4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagpass5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagpass6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagpass7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagpass8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet == 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagpass0leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagpass1leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagpass2leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagpass3leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagpass4leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagpass5leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagpass6leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagpass7leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagpass8leadjetpt0nJet->Fill(e_ee_invMass);
	}
	if (nJet == 1) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagpass0leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagpass1leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagpass2leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagpass3leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagpass4leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagpass5leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagpass6leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagpass7leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagpass8leadjetpt1nJet->Fill(e_ee_invMass);
	}
	if (nJet == 2) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagpass0leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagpass1leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagpass2leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagpass3leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagpass4leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagpass5leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagpass6leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagpass7leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagpass8leadjetpt2nJet->Fill(e_ee_invMass);
	}
	if (nJet == 3) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagpass0leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagpass1leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagpass2leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagpass3leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagpass4leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagpass5leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagpass6leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagpass7leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagpass8leadjetpt3nJet->Fill(e_ee_invMass);
	}
	if (nJet == 4) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagpass0leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagpass1leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagpass2leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagpass3leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagpass4leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagpass5leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagpass6leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagpass7leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagpass8leadjetpt4nJet->Fill(e_ee_invMass);
	}
	if (nJet >= 5) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagpass0leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagpass1leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagpass2leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagpass3leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagpass4leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagpass5leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagpass6leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagpass7leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagpass8leadjetpt5nJet->Fill(e_ee_invMass);
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
	if (nJet>4)  tagfail5jet->Fill(e_ee_invMass);
	if (numberOfVertices < 5) tagfail0pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 5 && numberOfVertices < 10) tagfail1pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 10 && numberOfVertices < 15) tagfail2pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 15 && numberOfVertices < 20) tagfail3pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 20) tagfail4pu->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) < 0.5) tagfail0eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 0.5 && fabs(secondptele->eta()) <1.0)  tagfail1eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 1.0 && fabs(secondptele->eta()) <1.4) tagfail2eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) > 1.6 && fabs(secondptele->eta()) <2.0)  tagfail3eta->Fill(e_ee_invMass);
	if (fabs(secondptele->eta()) >= 2.0 && fabs(secondptele->eta()) <3.0)  tagfail4eta->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagfail0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagfail1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagfail2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagfail3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagfail4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagfail5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagfail6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagfail7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagfail8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet == 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagfail0leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagfail1leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagfail2leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagfail3leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagfail4leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagfail5leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagfail6leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagfail7leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagfail8leadjetpt0nJet->Fill(e_ee_invMass);
	}
	if (nJet == 1) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagfail0leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagfail1leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagfail2leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagfail3leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagfail4leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagfail5leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagfail6leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagfail7leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagfail8leadjetpt1nJet->Fill(e_ee_invMass);
	}
	if (nJet == 2) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagfail0leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagfail1leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagfail2leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagfail3leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagfail4leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagfail5leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagfail6leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagfail7leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagfail8leadjetpt2nJet->Fill(e_ee_invMass);
	}
	if (nJet == 3) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagfail0leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagfail1leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagfail2leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagfail3leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagfail4leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagfail5leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagfail6leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagfail7leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagfail8leadjetpt3nJet->Fill(e_ee_invMass);
	}
	if (nJet == 4) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagfail0leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagfail1leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagfail2leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagfail3leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagfail4leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagfail5leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagfail6leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagfail7leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagfail8leadjetpt4nJet->Fill(e_ee_invMass);
	}
	if (nJet >= 5) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagfail0leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagfail1leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagfail2leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagfail3leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagfail4leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagfail5leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagfail6leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagfail7leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagfail8leadjetpt5nJet->Fill(e_ee_invMass);
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
      probeall_leadjetpt->Fill(leadingJet_pT);
      if (highestenergy_SC_isPassingProbe){
	probepass_pt->Fill(highestenergy_SC->energy());
	probepass_eta->Fill(highestenergy_SC->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (nJet==0) probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) probepass4jet->Fill(e_ee_invMass);
	if (nJet>4)  probepass5jet->Fill(e_ee_invMass);
	if (numberOfVertices < 5) probepass0pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 5 && numberOfVertices < 10) probepass1pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 10 && numberOfVertices < 15) probepass2pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 15 && numberOfVertices < 20) probepass3pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 20) probepass4pu->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) < 0.5) probepass0eta->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) >= 0.5 && fabs(highestenergy_SC->eta()) <1.0)  probepass1eta->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) >= 1.0 && fabs(highestenergy_SC->eta()) <1.4) probepass2eta->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) > 1.6 && fabs(highestenergy_SC->eta()) <2.0)  probepass3eta->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) >= 2.0 && fabs(highestenergy_SC->eta()) <3.0)  probepass4eta->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet == 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt0nJet->Fill(e_ee_invMass);
	}
	if (nJet == 1) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt1nJet->Fill(e_ee_invMass);
	}
	if (nJet == 2) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt2nJet->Fill(e_ee_invMass);
	}
	if (nJet == 3) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt3nJet->Fill(e_ee_invMass);
	}
	if (nJet == 4) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt4nJet->Fill(e_ee_invMass);
	}
	if (nJet >= 5) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt5nJet->Fill(e_ee_invMass);
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
	if (nJet>4)  probefail5jet->Fill(e_ee_invMass);
	if (numberOfVertices < 5) probefail0pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 5 && numberOfVertices < 10) probefail1pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 10 && numberOfVertices < 15) probefail2pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 15 && numberOfVertices < 20) probefail3pu->Fill(e_ee_invMass);
	if (numberOfVertices >= 20) probefail4pu->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) < 0.5) probefail0eta->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) >= 0.5 && fabs(highestenergy_SC->eta()) <1.0)  probefail1eta->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) >= 1.0 && fabs(highestenergy_SC->eta()) <1.4) probefail2eta->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) > 1.6 && fabs(highestenergy_SC->eta()) <2.0)  probefail3eta->Fill(e_ee_invMass);
	if (fabs(highestenergy_SC->eta()) >= 2.0 && fabs(highestenergy_SC->eta()) <3.0)  probefail4eta->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet == 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt0nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt0nJet->Fill(e_ee_invMass);
	}
	if (nJet == 1) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt1nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt1nJet->Fill(e_ee_invMass);
	}
	if (nJet == 2) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt2nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt2nJet->Fill(e_ee_invMass);
	}
	if (nJet == 3) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt3nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt3nJet->Fill(e_ee_invMass);
	}
	if (nJet == 4) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt4nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt4nJet->Fill(e_ee_invMass);
	}
	if (nJet >= 5) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt5nJet->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt5nJet->Fill(e_ee_invMass);
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
