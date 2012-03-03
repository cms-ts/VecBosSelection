#include <iostream>
#include "VecBosSelection/Selection/interface/SelectionUtils.h"



//DO the WP80 analysis
bool SelectionUtils::DoWP80(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
  double IsoTrk = 0;
  double IsoEcal = 0;
  double IsoHcal = 0;
  double HE = 0;

  //  if (recoElectron->et () <= 25) return false;


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
  
  if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
    
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

    if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	&& fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {

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

//DO the WP80Pf analysis
bool SelectionUtils::DoWP80Pf(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
   double IsoChg = 0;
   double IsoNeut = 0;
   double IsoPhot = 0;
   double IsoTot = 0;
   double HE=0.;
   
   //  if (recoElectron->et () <= 25) return false;
   
   HE = recoElectron->hadronicOverEm();
   if (removePU_){
      double lepIsoRho;
      
      /////// Pileup density "rho" for lepton isolation subtraction /////
      
      edm::Handle<double> rhoLepIso;
      const edm::InputTag eventrhoLepIso("kt6PFJetsForIsolation", "rho");
      iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
      if( *rhoLepIso == *rhoLepIso)  lepIsoRho = *rhoLepIso;
      else  lepIsoRho =  -999999.9;      
      IsoNeut = (recoElectron->pfIsolationVariables().neutralHadronIso - lepIsoRho*0.096);
      IsoChg = (recoElectron->pfIsolationVariables().chargedHadronIso - lepIsoRho*0.096);
      IsoPhot = (recoElectron->pfIsolationVariables().photonIso  - lepIsoRho*0.096);
      IsoTot = ((IsoNeut + IsoChg + IsoPhot)/ recoElectron->pt ());
   }
   else{
      // Define Isolation variables
      IsoNeut = recoElectron->pfIsolationVariables().neutralHadronIso;
      IsoChg = recoElectron->pfIsolationVariables().chargedHadronIso;
      IsoPhot = recoElectron->pfIsolationVariables().photonIso;
      IsoTot = ((IsoNeut + IsoChg + IsoPhot)/ recoElectron->pt ());
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
  
  if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
    
    /* Isolation */
    if (IsoTot < 0.20) {
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

    if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	&& fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {

      /* Isolation */  
      if (IsoTot < 0.20)  {
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




bool SelectionUtils::DoHLTMatch(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
  bool match=false;
  //DOesnt work, after chatting with Marco M.
  //for(std::vector<std::string>::const_iterator it = triggerNames_.begin(); it<triggerNames_.end();++it) {
  //string stringa=(string) *it;
  //if (recoElectron->triggerObjectMatchesByPath(stringa).size()>0) match=true;
  //if (Debug) cout<<"electron trigger Object Size ->"<<recoElectron->triggerObjectMatches().size()<<endl;
  //} 
  if (recoElectron->triggerObjectMatches().size()>0) match=true;

  return match;
}

std::vector<bool> SelectionUtils::MakeEleIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
  std::vector<bool> rec;
  double IsoTrk = 0;
  double IsoEcal = 0;
  double IsoHcal = 0;
  double HE = 0;

  //Ripetilo due volte, che ti fai anche in NON PU Rimosso
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

  if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
    isBarrelElectrons=true;
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

  if (isBarrelElectrons) {
    rec.push_back(isIsolatedBarrel);
    rec.push_back(isIDBarrel);
    rec.push_back(isConvertedBarrel);
    return rec;
  }

  /***** Endcap WP80 Cuts *****/

  if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
      && fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {
    isEndcapElectrons=true;
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
  
  if (isEndcapElectrons) {
    rec.push_back(isIsolatedEndcap);
    rec.push_back(isIDEndcap);
    rec.push_back(isConvertedEndcap);
    return rec;
  }

  //In the unlikely event of electrons in crack regions, return 0,0,0
    rec.push_back(0);
    rec.push_back(0);
    rec.push_back(0);
  return rec;
}
