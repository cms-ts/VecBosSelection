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
//         Created:  Wed May 11 14:53:26 CEST 2011
// $Id: EfficiencyFilter.cc,v 1.3 2011/09/05 09:15:32 marone Exp $
//
//


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
  Handle < GsfElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ())
    return false;

  //////////////////////
  //Match The HLT Trigger
  //////////////////////

  edm::Handle<edm::TriggerResults> HLTResults;
  iEvent.getByLabel(triggerCollection_, HLTResults);
  if (!HLTResults.isValid ()) return false;
  
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTResults);   
  bool flag=false;
    
 if (HLTResults.isValid()) {
   /// Storing the Prescale information: loop over the triggers and record prescale
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
	 if (Debug) cout<<"Matched "<<triggerNames.triggerName(itrig)<<endl;
	 int prescaleset = hltConfig_.prescaleSet(iEvent,iSetup);
	 if(prescaleset!=-1) {
	   prescalepair = hltConfig_.prescaleValues(iEvent,iSetup,triggerNames_[itrig]);
	   if (Debug) cout<<"prescale.first "<<prescalepair.first<<" prescalepair.second "<<prescalepair.second<<endl;
	   if((useCombinedPrescales_ && prescalepair.first<0) || prescalepair.second<0) {
	     edm::LogWarning("ZEfficiencyFilter") << " Unable to get prescale from event for trigger " << triggerNames.triggerName(itrig) << " :" 
					   << prescalepair.first << ", " << prescalepair.second;
	   }
	   prescale = useCombinedPrescales_ ? prescalepair.first*prescalepair.second : prescalepair.second;
	   minimalPrescale = minimalPrescale <  prescale ? minimalPrescale : prescale;
	   if (Debug) cout<<"prescale "<<prescale<<" minimal Prescale "<<minimalPrescale<<" for trigger "<<triggerNames.triggerName(itrig)<<endl;
	 } 
       }
     }
     else {
       // that trigger is presently not in the menu
       triggerSubset.push_back(false);
     }
   }
 }
  if (!flag) 
    {
      if(!useAllTriggers_) return false;
    }

    // ordino gli elettroni in pt 
  reco::GsfElectronCollection::const_iterator highestptele;
  reco::GsfElectronCollection::const_iterator secondptele;

  int i=0;
  
  if (electronCollection->size()==1) return false;

  for (reco::GsfElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    cout<<"Electron pt value ->"<<recoElectron->pt()<<endl;
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

  cout<<"First electron "<<highestptele->pt()<<" Second electron "<<secondptele->pt()<<endl;


  reco::GsfElectronCollection::const_iterator tag;
  reco::GsfElectronCollection::const_iterator probe;
  bool matchHLT=false;

  if ( DoHLTMatch(highestptele,iEvent) ) {
    tag=highestptele;
    probe=secondptele;
    matchHLT=true;
    cout<<"highest ele is matched by HLT"<<endl;
  }
  else{
    if ( DoHLTMatch(secondptele,iEvent) ) {
      tag=secondptele;
      probe=highestptele;
      matchHLT=true;
      cout<<"second ele is matched by HLT"<<endl;
    }
  }

  if (!matchHLT){
    cout<<"In this events, none of the two highest electrons (in pt) have triggered...exit"<<endl;
    return false;
  }

  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(tag->pt(),tag->eta(),tag->phi(), 0.0);
  TLorentzVector probev;
  probev.SetPtEtaPhiM(probe->pt(),probe->eta(),probe->phi(), 0.0);
  
  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();


  //Estraggo il contenuto dei jets
    int nJet = 0;
    int jetIndex = 0;
    Handle<PFJetCollection> pfJets;
    iEvent.getByLabel(theJetCollectionLabel_, pfJets);
    if (pfJets.isValid()) {
      PFJetCollection::const_iterator jet = pfJets->begin ();
      for (; jet != pfJets->end (); jet++, jetIndex++) {
	if (fabs(jet->eta())<2.4 && jet->pt()>30 && ((jet->eta()-tag->eta())>0.1 || (jet->phi()-tag->phi())>0.1 ) && ( (jet->eta()-probe->eta())>0.1 || (jet->phi()-probe->phi()>0.1)) ) {
	  nJet++;
	  cout<<"Jet eta "<<jet->eta()<<" pt "<<jet->pt()<<endl;
	}
      }
    }
    else{
      cout<<"No valid Jets Collection"<<endl;
    }
    
    cout<<"This event has jets #->"<<nJet<<endl;
  
  if ( DoWP80(tag,iEvent) ){
    cout<<"Tag is a WP80 electron..."<<endl;
  }
  else{
    cout<<"Tag IS a NOT WP80 electron...Exit"<<endl;
    return false;
  }

  if ( DoWP80(probe,iEvent) ){
    cout<<"Probe is a WP80 electron..."<<endl;
    probepass->Fill(e_ee_invMass);
    if (nJet==0) probepass0jet->Fill(e_ee_invMass);
    if (nJet==1) probepass1jet->Fill(e_ee_invMass);
    if (nJet==2) probepass2jet->Fill(e_ee_invMass);
    if (nJet==3) probepass3jet->Fill(e_ee_invMass);
    if (nJet==4) probepass4jet->Fill(e_ee_invMass);
    probeall->Fill(e_ee_invMass);
  }
  else{
    probefail->Fill(e_ee_invMass);
    probeall->Fill(e_ee_invMass);
    cout<<"Probe IS NOT a WP80 electron..."<<endl;
  }



  return true;
 
}


// ------------ method called once each job just before starting event loop  ------------
void
EfficiencyFilter::beginJob (){

//beginJob

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

bool EfficiencyFilter::DoWP80(reco::GsfElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
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

bool EfficiencyFilter::DoHLTMatch(reco::GsfElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
  //Check the electrons which have been triggered by the HLT
  edm::Handle<trigger::TriggerEvent> trgEvent;
  iEvent.getByLabel(InputTag("hltTriggerSummaryAOD","","HLT"), trgEvent);
  //Variables to be matched after
  float HLTpt=0;
  float HLTeta=0;
  float HLTphi=0;
  
  edm::InputTag myLastFilter = edm::InputTag("hltEle17CaloIdLCaloIsoVLPixelMatchFilter","","HLT");
    const trigger::TriggerObjectCollection& TOC( trgEvent->getObjects() );

    for(int i=0; i != trgEvent->sizeFilters(); ++i) {
      std::string label(trgEvent->filterTag(i).label());
      if (Debug) cout<<label<<endl;
      if (Debug) if( label == myLastFilter.label() ) cout<<"HT FIlter matched ->"<<label<<endl;;
    }

    if ( trgEvent->filterIndex(myLastFilter) < trgEvent->sizeFilters() ) {
      const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ) );
      for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) {
	int hltf = keys[hlto];
	const trigger::TriggerObject& L3obj(TOC[hltf]);
	HLTpt=L3obj.pt();
	HLTeta=L3obj.eta();
	HLTphi=L3obj.phi();
	if (Debug) cout<<"The matched HLT electron has pt,eta,phi ->"<<L3obj.pt()<<" "<<L3obj.eta()<<" "<<L3obj.phi()<<endl;
      }
    }
    if (Debug) cout<<"this electron has pt,eta,phi->"<<recoElectron->pt()<<" "<<recoElectron->eta()<<" "<<recoElectron->phi()<<endl;

    //Difference between mtached HLT electron and *Electron
    float diffpt=fabs(HLTpt-recoElectron->pt());
    float diffeta=fabs(HLTeta-recoElectron->eta());
    float diffphi=fabs(HLTphi-recoElectron->phi());
    cout<<"this electron has difference wrt HLT ele of pt,eta,phi->"<<diffpt<<" "<<diffeta<<" "<<diffphi<<endl;
    if (diffphi<0.1 && diffeta<0.1) {
      cout<<"This electron is the one who triggers the HLT!"<<endl;
      return true;
    }
    return false;
}


// ------------ method called once each job just after ending the event loop  ------------
void
EfficiencyFilter::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (EfficiencyFilter);
