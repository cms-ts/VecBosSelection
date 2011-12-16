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
// $Id: EfficiencyFilter.cc,v 1.8 2011/12/13 10:33:55 marone Exp $



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

  //Get Pattuple... Yes, we gave up....
  Handle < pat::ElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ()) return false;


  pat::ElectronCollection::const_iterator highestptele;
  pat::ElectronCollection::const_iterator secondptele;

  int i=0;
  int HLTmatches_SUM=0;
  int HLTmatches_TOTALELE=0;
  int HLTmatches_FAIL=0;

  if (electronCollection->size()==1) return false;
  bool protection=false;

  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    if (Debug) cout<<"Electron pt value ->"<<recoElectron->pt()<<endl;
    if (Debug) cout<<" MMMM ele trigger size "<<recoElectron->triggerObjectMatches().size()<<endl;
    if (i==0) highestptele=recoElectron;
    if ( DoHLTMatch(recoElectron,iEvent) ) HLTmatches_SUM++;
    if ( !DoHLTMatch(recoElectron,iEvent) ) HLTmatches_FAIL++;
    HLTmatches_TOTALELE++;
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

    HLTnumberOfMatches_PASS->Fill(HLTmatches_SUM);
    HLTnumberOfMatches_FAIL->Fill(HLTmatches_TOTALELE);
    HLTnumberOfMatches_TOTALELE->Fill(HLTmatches_FAIL);

  if (!protection) {
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }

  if (Debug) cout<<"First electron "<<highestptele->pt()<<" Second electron "<<secondptele->pt()<<endl;

  //Calculating tag & probe stuff
  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(highestptele->pt(),highestptele->eta(),highestptele->phi(), 0.0);
  TLorentzVector probev;
  probev.SetPtEtaPhiM(secondptele->pt(),secondptele->eta(),secondptele->phi(), 0.0);
  
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
	double deltaR1= sqrt( pow(jet->eta()-secondptele->eta(),2)+pow(jet->phi()-secondptele->phi(),2) );
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
    
    //Filling TAP distributions:



    if (DoWP80(highestptele,iEvent)){
      probept->Fill(secondptele->pt());
      probeall->Fill(e_ee_invMass);
      if ( DoWP80(secondptele,iEvent) ){
	if (Debug) cout<<"Probe is a WP80 electron->.."<<endl;
	probept_passWP80->Fill(secondptele->pt());
	WP80_probepass->Fill(e_ee_invMass);
	if (nJet==0) WP80_probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) WP80_probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) WP80_probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) WP80_probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) WP80_probepass4jet->Fill(e_ee_invMass);
      }
      else{
	probept_failWP80->Fill(secondptele->pt());
	WP80_probefail->Fill(e_ee_invMass);
	if (nJet==0) WP80_probefail0jet->Fill(e_ee_invMass);
	if (nJet==1) WP80_probefail1jet->Fill(e_ee_invMass);
	if (nJet==2) WP80_probefail2jet->Fill(e_ee_invMass);
	if (nJet==3) WP80_probefail3jet->Fill(e_ee_invMass);
	if (nJet==4) WP80_probefail4jet->Fill(e_ee_invMass);
      }
      if ( DoHLTMatch(secondptele,iEvent) ){
	if (Debug) cout<<"Probe is a WP80 electron->.."<<endl;
	HLT_probepass->Fill(e_ee_invMass);
	if (nJet==0) HLT_probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) HLT_probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) HLT_probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) HLT_probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) HLT_probepass4jet->Fill(e_ee_invMass);
      }
      else{
	HLT_probefail->Fill(e_ee_invMass);
	if (nJet==0) HLT_probefail0jet->Fill(e_ee_invMass);
	if (nJet==1) HLT_probefail1jet->Fill(e_ee_invMass);
	if (nJet==2) HLT_probefail2jet->Fill(e_ee_invMass);
	if (nJet==3) HLT_probefail3jet->Fill(e_ee_invMass);
	if (nJet==4) HLT_probefail4jet->Fill(e_ee_invMass);
      }
    }


    if (DoWP80(secondptele,iEvent)){
      tagpt->Fill(highestptele->pt());
      tagall->Fill(e_ee_invMass);
      if ( DoWP80(highestptele,iEvent) ){
	if (Debug) cout<<"Probe is a WP80 electron->.."<<endl;
	tagpt_passWP80->Fill(highestptele->pt());
	WP80_tagpass->Fill(e_ee_invMass);
	if (nJet==0) WP80_tagpass0jet->Fill(e_ee_invMass);
	if (nJet==1) WP80_tagpass1jet->Fill(e_ee_invMass);
	if (nJet==2) WP80_tagpass2jet->Fill(e_ee_invMass);
	if (nJet==3) WP80_tagpass3jet->Fill(e_ee_invMass);
	if (nJet==4) WP80_tagpass4jet->Fill(e_ee_invMass);
      }
      else{
	tagpt_failWP80->Fill(highestptele->pt());
	WP80_tagfail->Fill(e_ee_invMass);
	if (nJet==0) WP80_tagfail0jet->Fill(e_ee_invMass);
	if (nJet==1) WP80_tagfail1jet->Fill(e_ee_invMass);
	if (nJet==2) WP80_tagfail2jet->Fill(e_ee_invMass);
	if (nJet==3) WP80_tagfail3jet->Fill(e_ee_invMass);
	if (nJet==4) WP80_tagfail4jet->Fill(e_ee_invMass);
      }
      if ( DoHLTMatch(highestptele,iEvent) ){
	if (Debug) cout<<"Probe is a WP80 electron->.."<<endl;
	HLT_tagpass->Fill(e_ee_invMass);
	if (nJet==0) HLT_tagpass0jet->Fill(e_ee_invMass);
	if (nJet==1) HLT_tagpass1jet->Fill(e_ee_invMass);
	if (nJet==2) HLT_tagpass2jet->Fill(e_ee_invMass);
	if (nJet==3) HLT_tagpass3jet->Fill(e_ee_invMass);
	if (nJet==4) HLT_tagpass4jet->Fill(e_ee_invMass);
      }
      else{
	HLT_tagfail->Fill(e_ee_invMass);
	if (nJet==0) HLT_tagfail0jet->Fill(e_ee_invMass);
	if (nJet==1) HLT_tagfail1jet->Fill(e_ee_invMass);
	if (nJet==2) HLT_tagfail2jet->Fill(e_ee_invMass);
	if (nJet==3) HLT_tagfail3jet->Fill(e_ee_invMass);
	if (nJet==4) HLT_tagfail4jet->Fill(e_ee_invMass);
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

//DO the WP80 analysis

bool EfficiencyFilter::DoWP80(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
  double IsoTrk = 0;
  double IsoEcal = 0;
  double IsoHcal = 0;
  double HE = 0;
  
  if (removePU_){
    double lepIsoRho;
    
    /////// Pileup density "rho" for lepton isolation subtraction /////
    
    edm::Handle<double> rhoLepIso;
    const edm::InputTag eventrhoLepIso("kt6PFJetsForIsolation", "rho");
    iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
    if( *rhoLepIso == *rhoLepIso)  lepIsoRho = *rhoLepIso;
    else  lepIsoRho =  -999999.9;      
    IsoEcal = (recoElectron->dr03EcalRecHitSumEt () - lepIsoRho*0.096) / recoElectron->et ();
    IsoTrk = (recoElectron->dr03TkSumPt () - lepIsoRho*0.096) / recoElectron->et ();
    IsoHcal = (recoElectron->dr03HcalTowerSumEt ()  - lepIsoRho*0.096) / recoElectron->et ();
      HE = recoElectron->hadronicOverEm();
  }
  else{
    // Define Isolation variables
    IsoTrk = (recoElectron->dr03TkSumPt () / recoElectron->et ());
    IsoEcal = (recoElectron->dr03EcalRecHitSumEt () / recoElectron->et ());
    IsoHcal = (recoElectron->dr03HcalTowerSumEt () / recoElectron->et ());
    HE = recoElectron->hadronicOverEm();
  }
  //Define ID variables
  
  float DeltaPhiTkClu = recoElectron->deltaPhiSuperClusterTrackAtVtx ();
  float DeltaEtaTkClu = recoElectron->deltaEtaSuperClusterTrackAtVtx ();
  float sigmaIeIe = recoElectron->sigmaIetaIeta ();
  
  //Define Conversion Rejection Variables
  
  float Dcot = recoElectron->convDcot ();
  float Dist = recoElectron->convDist ();
  int NumberOfExpectedInnerHits = recoElectron->gsfTrack ()->trackerExpectedHitsInner ().numberOfHits ();
  
  //quality flags

  bool isBarrelElectrons;
  bool isEndcapElectrons;
  bool isIsolatedBarrel;
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIsolatedEndcap;
  bool isIDEndcap;
  bool isConvertedEndcap;
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
    return true;
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
      return true;
    }
    return false;
}


bool EfficiencyFilter::DoHLTMatch(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
  bool match=false;
  //DOesnt work, after chatting with Marco M.
  //for(std::vector<std::string>::const_iterator it = triggerNames_.begin(); it<triggerNames_.end();++it) {
  //string stringa=(string) *it;
  //if (recoElectron->triggerObjectMatchesByPath(stringa).size()>0) match=true;
  if (Debug) cout<<"electron trigger Object Size ->"<<recoElectron->triggerObjectMatches().size()<<endl;
  //} 
  if (recoElectron->triggerObjectMatches().size()>0) match=true;

  return match;
}


// ------------ method called once each job just after ending the event loop  ------------
void
EfficiencyFilter::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (EfficiencyFilter);
