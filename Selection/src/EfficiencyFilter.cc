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
// $Id$
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


using namespace std;
using namespace edm;
using namespace reco;
bool Debug=false; //Activate with true if you wonna have verbosity for Debug


//
// constructors and destructor
//
EfficiencyFilter::EfficiencyFilter (const edm::ParameterSet & parameters)
{
  theElectronCollectionLabel =
    parameters.getParameter < InputTag > ("electronCollection");
  std::string outputfile_D = parameters.getUntrackedParameter<std::string>("filename");
  outputfile_ = parameters.getUntrackedParameter<std::string>("outputfile", outputfile_D);
  triggerCollection_=parameters.getUntrackedParameter<edm::InputTag>("triggerCollectionTag");
  useCombinedPrescales_ = parameters.getParameter<bool>("UseCombinedPrescales");
  triggerNames_         = parameters.getParameter< std::vector<std::string> > ("TriggerNames");
  useAllTriggers_       = (triggerNames_.size()==0);
  removePU_             = parameters.getParameter<bool>("removePU");
  electronIsolatedProducer_ = parameters.getParameter< edm::InputTag > ("electronIsolatedProducer");
  candTag_ = parameters.getParameter< edm::InputTag > ("candTag");
}



EfficiencyFilter::~EfficiencyFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}




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

  //Check the electrons which have been triggered by the HLT
    edm::Handle<trigger::TriggerEvent> trgEvent;
    iEvent.getByLabel(InputTag("hltTriggerSummaryAOD","","HLT"), trgEvent);
    int index=0;

    edm::InputTag myLastFilter = edm::InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt20PixelMatchFilter","","HLT");
    const trigger::TriggerObjectCollection& TOC( trgEvent->getObjects() );
    //const TriggerObjectCollection TOC(trgEvent->getObjects());
    // filterIndex must be less than the size of trgEvent or you get a CMSException: _M_range_check
    for(int i=0; i != trgEvent->sizeFilters(); ++i) {
      std::string label(trgEvent->filterTag(i).label());
      if( label == myLastFilter.label() ) index = i;
    }

    if ( trgEvent->filterIndex(myLastFilter) < trgEvent->sizeFilters() ) {
      const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ) );
      for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) {
	int hltf = keys[hlto];
	const trigger::TriggerObject& L3obj(TOC[hltf]);
	float pL3obj = L3obj.p(); 
      }
    }
	//   using namespace trigger;
	//std::auto_ptr<trigger::TriggerFilterObjectWithRefs> filterproduct (new trigger::TriggerFilterObjectWithRefs(path(),module()));
	//filterproduct->addCollectionTag(electronIsolatedProducer_);
   //will be a collection of Ref<reco::ElectronCollection> ref;

	//edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
	//iEvent.getByLabel (candTag_,PrevFilterOutput);

	//std::vector<edm::Ref<reco::RecoEcalCandidateCollection> > recoecalcands;
	//PrevFilterOutput->getObjects(TriggerCluster, recoecalcands);

	//  edm::Handle<ElectronCollection> electronIsolatedHandle;
	//iEvent.getByLabel(electronIsolatedProducer_,electronIsolatedHandle);

  for (reco::GsfElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {

    //loop over the electrons to find the matching one
    //for(reco::ElectronCollection::const_iterator iElectron = electronIsolatedHandle->begin(); iElectron != electronIsolatedHandle->end(); iElectron++){
    //reco::ElectronRef electronref(reco::ElectronRef(electronIsolatedHandle,iElectron - electronIsolatedHandle->begin()));
    //const reco::SuperClusterRef theClus = electronref->superCluster();
    //reco::SuperClusterRef recr2 = recoElectron->superCluster();
      
    //if(&(*recr2) ==  &(*theClus)) {
    //} 
    //}

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
      TLorentzVector b_e2(recoElectron->momentum ().x (),recoElectron->momentum ().y (),recoElectron->momentum ().z (), recoElectron->p ());
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
      TLorentzVector e_e2(recoElectron->momentum ().x (),recoElectron->momentum ().y (),recoElectron->momentum ().z (), recoElectron->p ());
      LV.push_back(e_e2);
    }

  }
  if (elIsAccepted<=0)    return false;

  return true;
 
}


// ------------ method called once each job just before starting event loop  ------------
void
EfficiencyFilter::beginJob (){
  fOfile = new TFile("EfficiencyFilter.root","RECREATE");
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
  fOfile->cd();
  fOfile->Write() ;
  fOfile->Close() ;
}

DEFINE_FWK_MODULE (EfficiencyFilter);
