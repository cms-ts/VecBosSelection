#include <iostream>
#include "VecBosSelection/Selection/interface/SelectionUtils.h"

bool HLTmatch=true; //This one is used as a test, to compare with analyses without HLT ELE matching... Should be TRUE!

//DO the WP80 analysis
bool SelectionUtils::DoWP80(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_)
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
   double HE=0.;
   
   //  if (recoElectron->et () <= 25) return false;
   
   HE = recoElectron->hadronicOverEm();
   
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
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIDEndcap;
  bool isConvertedEndcap;
  isBarrelElectrons = false;
  isEndcapElectrons = false;
  isIDBarrel = false;
  isConvertedBarrel = false;
  isIDEndcap = false;
  isConvertedEndcap = false;
 

/***** Barrel WP80 Cuts *****/
  
  if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
    
    
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
  
  if (isIDBarrel && isConvertedBarrel) {
    return true;
  }

    /***** Endcap WP80 Cuts *****/

    if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	&& fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {


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
    
    if (isIDEndcap && isConvertedEndcap) {
      return true;
    }
    return false;
}

//DO the WP90Pf analysis
bool SelectionUtils::DoWP90Pf(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
   double HE=0.;
   
   //  if (recoElectron->et () <= 25) return false;
   
   HE = recoElectron->hadronicOverEm();
   
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
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIDEndcap;
  bool isConvertedEndcap;
  isBarrelElectrons = false;
  isEndcapElectrons = false;
  isIDBarrel = false;
  isConvertedBarrel = false;
  isIDEndcap = false;
  isConvertedEndcap = false;
 

/***** Barrel WP80 Cuts *****/
  
  if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
    
    
    /* Identification */
    if (fabs (DeltaEtaTkClu) < 0.007 && fabs (DeltaPhiTkClu) < 0.8
	&& sigmaIeIe < 0.01 && HE < 0.12) {
      isIDBarrel = true;
    }
    
    /* Conversion Rejection */
      if ((fabs (Dist) >= 0.02 || fabs (Dcot) >= 0.02)
	  && NumberOfExpectedInnerHits <= 1) {
	isConvertedBarrel = true;
      }
      //cout<<"isIsolatedBarrel "<<isIsolatedBarrel<<" isIDBarrel "<< isIDBarrel<<" isConvertedBarrel "<<isConvertedBarrel<<endl;
  }
  
  if (isIDBarrel && isConvertedBarrel) {
    return true;
  }

    /***** Endcap WP80 Cuts *****/

    if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	&& fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {


      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.009 && fabs (DeltaPhiTkClu) < 0.7
	  && sigmaIeIe < 0.03 && HE < 0.15) {
	isIDEndcap = true;
      }

      /* Conversion Rejection */
      if ((fabs (Dcot) > 0.02 || fabs (Dist) > 0.02)
	  && NumberOfExpectedInnerHits <=1) {
	isConvertedEndcap = true;
      }
      //cout<<"isIsolatedEndcap "<<isIsolatedEndcap<<" isIDEndcap "<< isIDEndcap<<" isConvertedEndcap "<<isConvertedEndcap<<endl;
    }
    
    if (isIDEndcap && isConvertedEndcap) {
      return true;
    }
    return false;
}

//DO the WP90Pf analysis
bool SelectionUtils::DoWP90PfGSF(reco::GsfElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
   double HE=0.;
   
   //  if (recoElectron->et () <= 25) return false;
   
   HE = recoElectron->hadronicOverEm();
   
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
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIDEndcap;
  bool isConvertedEndcap;
  isBarrelElectrons = false;
  isEndcapElectrons = false;
  isIDBarrel = false;
  isConvertedBarrel = false;
  isIDEndcap = false;
  isConvertedEndcap = false;
 

/***** Barrel WP80 Cuts *****/
  
  if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
    
    
    /* Identification */
    if (fabs (DeltaEtaTkClu) < 0.007 && fabs (DeltaPhiTkClu) < 0.8
	&& sigmaIeIe < 0.01 && HE < 0.12) {
      isIDBarrel = true;
    }
    
    /* Conversion Rejection */
      if ((fabs (Dist) >= 0.02 || fabs (Dcot) >= 0.02)
	  && NumberOfExpectedInnerHits <= 1) {
	isConvertedBarrel = true;
      }
      //cout<<"isIsolatedBarrel "<<isIsolatedBarrel<<" isIDBarrel "<< isIDBarrel<<" isConvertedBarrel "<<isConvertedBarrel<<endl;
  }
  
  if (isIDBarrel && isConvertedBarrel) {
    return true;
  }

    /***** Endcap WP80 Cuts *****/

    if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	&& fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {


      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.009 && fabs (DeltaPhiTkClu) < 0.7
	  && sigmaIeIe < 0.03 && HE < 0.15) {
	isIDEndcap = true;
      }

      /* Conversion Rejection */
      if ((fabs (Dcot) > 0.02 || fabs (Dist) > 0.02)
	  && NumberOfExpectedInnerHits <=1) {
	isConvertedEndcap = true;
      }
      //cout<<"isIsolatedEndcap "<<isIsolatedEndcap<<" isIDEndcap "<< isIDEndcap<<" isConvertedEndcap "<<isConvertedEndcap<<endl;
    }
    
    if (isIDEndcap && isConvertedEndcap) {
      return true;
    }
    return false;
}

//DO the WP80Pf analysis with NEW HE VARIABLE!! (testing!!!)
bool SelectionUtils::DoWP80Pf_NewHE(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_)
{
   double HE=0.;
   
   //  if (recoElectron->et () <= 25) return false;
   
   // New definition of HE variable (!!!):
   HE = recoElectron->hcalOverEcalBc();

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
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIDEndcap;
  bool isConvertedEndcap;
  isBarrelElectrons = false;
  isEndcapElectrons = false;
  isIDBarrel = false;
  isConvertedBarrel = false;
  isIDEndcap = false;
  isConvertedEndcap = false;
 

/***** Barrel WP80 Cuts *****/
  
  if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
    
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
  
  if (isIDBarrel && isConvertedBarrel) {
    return true;
  }

    /***** Endcap WP80 Cuts *****/

    if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	&& fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {

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
    
    if (isIDEndcap && isConvertedEndcap) {
      return true;
    }
    return false;
}




bool SelectionUtils::DoHLTMatch(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent)
{
  if (!HLTmatch) {
    cout<<"WARNING: HLT matched in the event, not checking if both eles have been triggered! Go in SelectionUtils.cc to change it"<<endl;
    return true;
  }
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

std::vector<bool> SelectionUtils::MakeEleIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_)
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
     
     isBarrelElectrons = true;
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
   //isIsolatedBarrel = true;
    rec.push_back(isIsolatedBarrel);
    rec.push_back(isIDBarrel);
    rec.push_back(isConvertedBarrel);
    return rec;
  }

  /***** Endcap WP80 Cuts *****/

  if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
      && fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {    
     isEndcapElectrons = true;
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
     //isIsolatedEndcap = true;
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


std::vector<bool> SelectionUtils::MakePfEleIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_)
{
  std::vector<bool> rec;
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
      //cout << "1- removing PU ..."<<endl;
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
      
      //cout << "2- NOT removing PU ..."<<endl;
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
     isBarrelElectrons=true;
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

// ================================================ 
//|                                                |
//|    N E W    S E L E C T I O N S   x   2011     |
//|                                                |
// ================================================

//DO the MEDIUM cuts VBTF 2011 analysis
bool SelectionUtils::DoMedSel2011(pat::ElectronCollection::const_iterator recoElectron,
				  edm::Event& iEvent,
				  const edm::Handle<reco::ConversionCollection> &conversions,
				  const reco::BeamSpot &beamspot,
				  const edm::Handle<reco::VertexCollection> &vtxs)
{
  

   //Define ID variables
   
   float DeltaPhiTkClu = recoElectron->deltaPhiSuperClusterTrackAtVtx ();
   float DeltaEtaTkClu = recoElectron->deltaEtaSuperClusterTrackAtVtx ();
   float sigmaIeIe     = recoElectron->sigmaIetaIeta (); 
   float HE            = recoElectron->hadronicOverEm();
   float ooemoop       = (1.0/recoElectron->ecalEnergy() - recoElectron->eSuperClusterOverP()/recoElectron->ecalEnergy());

   // impact parameter variables
   float d0vtx         = 0.0;
   float dzvtx         = 0.0;
   if (vtxs->size() > 0) {
      reco::VertexRef vtx(vtxs, 0);    
      d0vtx = recoElectron->gsfTrack()->dxy(vtx->position());
      dzvtx = recoElectron->gsfTrack()->dz(vtx->position());
   } else {
      d0vtx = recoElectron->gsfTrack()->dxy();
      dzvtx = recoElectron->gsfTrack()->dz();
   }
   
   // conversion rejection variables
   
   bool vtxFitConversion = SelectionUtils::hasMatchedConversion(*recoElectron, conversions, beamspot.position());
   float mHits           = recoElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits(); 
   
   //quality flags
   
   bool isIDBarrel;
   bool isConvertedBarrel;
   bool isIDEndcap;
   bool isConvertedEndcap;
   isIDBarrel = false;
   isConvertedBarrel = false;
   isIDEndcap = false;
   isConvertedEndcap = false;
   
   
/***** Barrel Medium 2011 Cuts *****/
  
   if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
      
      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.004 && fabs (DeltaPhiTkClu) < 0.06
	  && sigmaIeIe < 0.01 && HE < 0.12 && fabs (ooemoop) < 0.05
	  && fabs(d0vtx) < 0.02 && fabs(dzvtx) < 0.1) {
	 isIDBarrel = true;
      }
      
      /* Conversion Rejection */
      if ( !vtxFitConversion && mHits <=1) {
	 isConvertedBarrel = true;
      }
      //cout<<"isIsolatedBarrel "<<isIsolatedBarrel<<" isIDBarrel "<< isIDBarrel<<" isConvertedBarrel "<<isConvertedBarrel<<endl;
   }
   
   if (isIDBarrel && isConvertedBarrel) {
      return true;
   }

    /***** Endcap Medium 2011 Cuts *****/

    if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	&& fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {

      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.007 && fabs (DeltaPhiTkClu) < 0.03
	&& sigmaIeIe < 0.03 && HE < 0.10 && fabs (ooemoop) < 0.05
	&& fabs (d0vtx) < 0.02 && fabs (dzvtx) < 0.1) {
	isIDEndcap = true;
      }

      /* Conversion Rejection */
      if (!vtxFitConversion && mHits <=1) {
	isConvertedEndcap = true;
      }
      //cout<<"isIsolatedEndcap "<<isIsolatedEndcap<<" isIDEndcap "<< isIDEndcap<<" isConvertedEndcap "<<isConvertedEndcap<<endl;
    }
    
    if (isIDEndcap && isConvertedEndcap) {
      return true;
    }
    return false;
}


bool SelectionUtils::isGoodConversion(const reco::Conversion &conv, const math::XYZPoint &beamspot, float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax)
{
  
  //Check if a given conversion candidate passes the conversion selection cuts
  
  const reco::Vertex &vtx = conv.conversionVertex();

  //vertex validity
  if (!vtx.isValid()) return false;

  //fit probability
  if (TMath::Prob( vtx.chi2(),  vtx.ndof() )<probMin) return false;

  //compute transverse decay length
  math::XYZVector mom(conv.refittedPairMomentum()); 
  double dbsx = vtx.x() - beamspot.x();
  double dbsy = vtx.y() - beamspot.y();
  double lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();

  //transverse decay length  
  if ( lxy<lxyMin )
    return false;
    
  //loop through daughters to check nhitsbeforevtx
  for (std::vector<uint8_t>::const_iterator it = conv.nHitsBeforeVtx().begin(); it!=conv.nHitsBeforeVtx().end(); ++it) {
    if ( (*it)>nHitsBeforeVtxMax ) return false;
  }
  
  return true;
}

bool SelectionUtils::matchesConversion(const pat::Electron &ele, const reco::Conversion &conv, bool allowCkfMatch)
{

  //check if a given GsfElectron matches a given conversion (no quality cuts applied)
  //matching is always attempted through the gsf track ref, and optionally attempted through the
  //closest ctf track ref

  const std::vector<edm::RefToBase<reco::Track> > &convTracks = conv.tracks();
  for (std::vector<edm::RefToBase<reco::Track> >::const_iterator it=convTracks.begin(); it!=convTracks.end(); ++it) {
    if ( ele.gsfTrack().isNonnull() && ele.gsfTrack().id()==it->id() && ele.gsfTrack().key()==it->key()) return true;
    else if ( allowCkfMatch && ele.closestCtfTrackRef().isNonnull() && ele.closestCtfTrackRef().id()==it->id() && ele.closestCtfTrackRef().key()==it->key() ) return true;
  }

  return false;
}

bool SelectionUtils::hasMatchedConversion(const pat::Electron &ele,
					  const edm::Handle<reco::ConversionCollection> &convCol,
					  const math::XYZPoint &beamspot, bool allowCkfMatch, float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax)
{
  //check if a given electron candidate matches to at least one conversion candidate in the
  //collection which also passes the selection cuts, optionally match with the closestckf track in
  //in addition to just the gsf track (enabled in default arguments)
  
  for (reco::ConversionCollection::const_iterator it = convCol->begin(); it!=convCol->end(); ++it) {
    if (!matchesConversion(ele, *it, allowCkfMatch)) continue;
    if (!isGoodConversion(*it,beamspot,lxyMin,probMin,nHitsBeforeVtxMax)) continue;
   
    return true;
  }
  
  return false;
  
}

// ************ Isolation cuts ***************

//bool SelectionUtils::DoIso2011(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent, std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >  &IsoMap, std::vector< edm::Handle< edm::ValueMap<double> > > &IsoVal){
bool SelectionUtils::DoIso2011(pat::ElectronCollection::const_iterator recoElectron, edm::Event& iEvent, std::vector< edm::Handle< edm::ValueMap<double> > > &IsoVals){

   double lepIsoRho;
  
   /////// Pileup density "rho" for lepton isolation subtraction /////
   //cout << "1- removing PU ..."<<endl;
   edm::Handle<double> rhoLepIso;
   const edm::InputTag eventrhoLepIso("kt6PFJetsForIsolation", "rho");
   iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
   if( *rhoLepIso == *rhoLepIso)  lepIsoRho = *rhoLepIso;
   else  lepIsoRho =  -999999.9;     

   
   //const std::vector< edm::Handle< edm::ValueMap<double> > > * electronIsoVals = &IsoVals;
   
    // use the reference to the original gsfElectron from PAT::electron

    // double charged =  (*(*electronIsoVals)[0])[recoElectron->originalObjectRef()];//myElectronRef
//     double photon = (*(*electronIsoVals)[1])[recoElectron->originalObjectRef()]; //myElectronRef
//     double neutral = (*(*electronIsoVals)[2])[recoElectron->originalObjectRef()];//myElectronRef 
    double charged = (*IsoVals[0])[recoElectron->originalObjectRef()];//myElectronRef
    double photon  = (*IsoVals[1])[recoElectron->originalObjectRef()]; //myElectronRef
    double neutral = (*IsoVals[2])[recoElectron->originalObjectRef()];//myElectronRef 

    // Effective area for 2011 data (Delta_R=0.3) (taken from https://twiki.cern.ch/twiki/bin/view/Main/HVVElectronId2012 )
    double A_eff_PH, A_eff_NH;
    if(abs(recoElectron->eta())<=1.0){A_eff_PH=0.081; A_eff_NH=0.024;}
    else if(abs(recoElectron->eta())>1.0 && abs(recoElectron->eta())<=1.479){A_eff_PH=0.084 ; A_eff_NH=0.037;}
    else if(abs(recoElectron->eta())>1.479 && abs(recoElectron->eta())<=2.0){A_eff_PH=0.048 ; A_eff_NH=0.037;}
    else if(abs(recoElectron->eta())>2.0 && abs(recoElectron->eta())<=2.2){A_eff_PH=0.089 ; A_eff_NH=0.023;}
    else if(abs(recoElectron->eta())>2.2 && abs(recoElectron->eta())<=2.3){A_eff_PH=0.092 ; A_eff_NH=0.023;}
    else if(abs(recoElectron->eta())>2.3 && abs(recoElectron->eta())<=2.4){A_eff_PH=0.097 ; A_eff_NH=0.021;}   
    else {A_eff_PH=0.11 ; A_eff_NH=0.021;}
    
    float PFIsoPUCorrected =(charged + max(photon - lepIsoRho*A_eff_PH  , 0.) +  max(neutral - lepIsoRho * A_eff_NH, 0.))/std::max(0.5, recoElectron->pt());

    if (PFIsoPUCorrected < 0.15) return true;
    
    return false;

}

std::vector<bool> SelectionUtils::MakePfEleNewIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool useNewID_,bool doWP90_,const edm::Handle<reco::ConversionCollection> &conversions,const reco::BeamSpot &beamspot,const edm::Handle<reco::VertexCollection> &vtxs, IsoDepositVals &IsoVals)
{
  std::vector<bool> rec;

  // Define all the cuts

  float mHitsCut, distCut, dcotCut, sigmaEBCut, sigmaEECut, dPhiEBCut, dPhiEECut, dEtaEBCut, dEtaEECut, hoeEBCut, hoeEECut;
  if (!doWP90_){
     mHitsCut= 0;
     distCut= 0.02;
     dcotCut= 0.02; 
     sigmaEBCut= 0.01;
     sigmaEECut= 0.03;
     dPhiEBCut= 0.06;
     dPhiEECut= 0.03;
     dEtaEBCut= 0.004;
     dEtaEECut= 0.007;
     hoeEBCut= 0.04;
     hoeEECut= 0.15;
  } else {
     mHitsCut= 1;
     distCut= 0.02;
     dcotCut= 0.02; 
     sigmaEBCut= 0.01;
     sigmaEECut= 0.03;
     dPhiEBCut= 0.8;
     dPhiEECut= 0.7;
     dEtaEBCut= 0.007;
     dEtaEECut= 0.009;
     hoeEBCut= 0.12;
     hoeEECut= 0.15;
  } 

  // Define ISOLATION variable
   double lepIsoRho;
  
   /////// Pileup density "rho" for lepton isolation subtraction /////
   //cout << "1- removing PU ..."<<endl;
   edm::Handle<double> rhoLepIso;
   const edm::InputTag eventrhoLepIso("kt6PFJetsForIsolation", "rho");
   iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
   if( *rhoLepIso == *rhoLepIso)  lepIsoRho = *rhoLepIso;
   else  lepIsoRho =  -999999.9;     

   
   //const std::vector< edm::Handle< edm::ValueMap<double> > > * electronIsoVals = &IsoVals;
   
    // use the reference to the original gsfElectron from PAT::electron

    // double charged =  (*(*electronIsoVals)[0])[recoElectron->originalObjectRef()];//myElectronRef
//     double photon = (*(*electronIsoVals)[1])[recoElectron->originalObjectRef()]; //myElectronRef
//     double neutral = (*(*electronIsoVals)[2])[recoElectron->originalObjectRef()];//myElectronRef 
    double charged = (*IsoVals[0])[recoElectron->originalObjectRef()];//myElectronRef
    double photon  = (*IsoVals[1])[recoElectron->originalObjectRef()]; //myElectronRef
    double neutral = (*IsoVals[2])[recoElectron->originalObjectRef()];//myElectronRef 

    // Effective area for 2011 data (Delta_R=0.3) (taken from https://twiki.cern.ch/twiki/bin/view/Main/HVVElectronId2012 )
    double A_eff_PH, A_eff_NH;
    if(abs(recoElectron->eta())<=1.0){A_eff_PH=0.081; A_eff_NH=0.024;}
    else if(abs(recoElectron->eta())>1.0 && abs(recoElectron->eta())<=1.479){A_eff_PH=0.084 ; A_eff_NH=0.037;}
    else if(abs(recoElectron->eta())>1.479 && abs(recoElectron->eta())<=2.0){A_eff_PH=0.048 ; A_eff_NH=0.037;}
    else if(abs(recoElectron->eta())>2.0 && abs(recoElectron->eta())<=2.2){A_eff_PH=0.089 ; A_eff_NH=0.023;}
    else if(abs(recoElectron->eta())>2.2 && abs(recoElectron->eta())<=2.3){A_eff_PH=0.092 ; A_eff_NH=0.023;}
    else if(abs(recoElectron->eta())>2.3 && abs(recoElectron->eta())<=2.4){A_eff_PH=0.097 ; A_eff_NH=0.021;}   
    else {A_eff_PH=0.11 ; A_eff_NH=0.021;}
    
    float PFIsoPUCorrected =(charged + max(photon - lepIsoRho*A_eff_PH  , 0.) +  max(neutral - lepIsoRho * A_eff_NH, 0.))/std::max(0.5, recoElectron->pt());

  //Define ID variables                            
  float DeltaPhiTkClu = recoElectron->deltaPhiSuperClusterTrackAtVtx ();
  float DeltaEtaTkClu = recoElectron->deltaEtaSuperClusterTrackAtVtx ();
  float sigmaIeIe     = recoElectron->sigmaIetaIeta ();
  float HE            = recoElectron->hadronicOverEm();
  float ooemoop       = (1.0/recoElectron->ecalEnergy() - recoElectron->eSuperClusterOverP()/recoElectron->ecalEnergy());

  //Define Conversion Rejection Variables   
  float Dcot = recoElectron->convDcot ();
  float Dist = recoElectron->convDist ();
  int NumberOfExpectedInnerHits = recoElectron->gsfTrack ()->trackerExpectedHitsInner ().numberOfHits ();
  bool vtxFitConversion = SelectionUtils::hasMatchedConversion(*recoElectron, conversions, beamspot.position());


   // impact parameter variables
   float d0vtx         = 0.0;
   float dzvtx         = 0.0;
   if (vtxs->size() > 0) {
      reco::VertexRef vtx(vtxs, 0);    
      d0vtx = recoElectron->gsfTrack()->dxy(vtx->position());
      dzvtx = recoElectron->gsfTrack()->dz(vtx->position());
   } else {
      d0vtx = recoElectron->gsfTrack()->dxy();
      dzvtx = recoElectron->gsfTrack()->dz();
   }
  

  //quality flags                                                                                                                                            
  bool isBarrelElectrons;
  bool isEndcapElectrons;
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIDEndcap;
  bool isConvertedEndcap;  
  bool isIsolated;
  isBarrelElectrons = false;
  isEndcapElectrons = false;
  isIDBarrel = false;
  isConvertedBarrel = false;
  isIDEndcap = false;
  isConvertedEndcap = false;
  isIsolated = false;
   
  // *************** isolation *************  
  // ---------------------------------------
  if (PFIsoPUCorrected < 0.15) isIsolated= true;

  if (!useNewID_){
     /***** Barrel WP80 Cuts *****/
     
     if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
	isBarrelElectrons=true;
	
	/* Identification */
	if (fabs (DeltaEtaTkClu) < dEtaEBCut && fabs (DeltaPhiTkClu) < dPhiEBCut
	    && sigmaIeIe < sigmaEBCut && HE < hoeEBCut) {
	   isIDBarrel = true;
	}
	
	/* Conversion Rejection */
	if ((fabs (Dist) >= distCut || fabs (Dcot) >= dcotCut)
	    && NumberOfExpectedInnerHits <= mHitsCut) {
	   isConvertedBarrel = true;
	}
	//cout<<"isIsolatedBarrel "<<isIsolatedBarrel<<" isIDBarrel "<< isIDBarrel<<" isConvertedBarrel "<<isConvertedBarrel<<endl;                            
     }
     
     /***** Endcap WP80 Cuts *****/
     
     if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	 && fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {
	isEndcapElectrons=true;
	
	/* Identification */
	if (fabs (DeltaEtaTkClu) < dEtaEECut && fabs (DeltaPhiTkClu) < dPhiEECut
	    && sigmaIeIe < sigmaEECut && HE < hoeEECut) {
	   isIDEndcap = true;
	}
	
	/* Conversion Rejection */
	if ((fabs (Dcot) > dcotCut || fabs (Dist) > distCut)
	    && NumberOfExpectedInnerHits <= mHitsCut) {
	   isConvertedEndcap = true;
	}
	//cout<<"isIsolatedEndcap "<<isIsolatedEndcap<<" isIDEndcap "<< isIDEndcap<<" isConvertedEndcap "<<isConvertedEndcap<<endl;                            
     }
  } else {
     
/***** Barrel Medium 2011 Cuts *****/
     
     if (fabs (recoElectron ->superCluster()->eta()) <= 1.4442) {
	isBarrelElectrons=true;
	
	/* Identification */
	if (fabs (DeltaEtaTkClu) < 0.004 && fabs (DeltaPhiTkClu) < 0.06
	    && sigmaIeIe < 0.01 && HE < 0.12 && fabs (ooemoop) < 0.05
	    && fabs(d0vtx) < 0.02 && fabs(dzvtx) < 0.1) {
	   isIDBarrel = true;
	}
	
	/* Conversion Rejection */
	if ( !vtxFitConversion && NumberOfExpectedInnerHits <=1) {
	   isConvertedBarrel = true;
	}
	//cout<<"isIsolatedBarrel "<<isIsolatedBarrel<<" isIDBarrel "<< isIDBarrel<<" isConvertedBarrel "<<isConvertedBarrel<<endl;
     }
     

    /***** Endcap Medium 2011 Cuts *****/

    if (fabs (recoElectron ->superCluster()->eta()) >= 1.5660
	&& fabs (recoElectron ->superCluster()->eta()) <= 2.5000) {
       isEndcapElectrons=true;

      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.007 && fabs (DeltaPhiTkClu) < 0.03
	&& sigmaIeIe < 0.03 && HE < 0.10 && fabs (ooemoop) < 0.05
	&& fabs (d0vtx) < 0.02 && fabs (dzvtx) < 0.1) {
	isIDEndcap = true;
      }

      /* Conversion Rejection */
      if (!vtxFitConversion && NumberOfExpectedInnerHits <=1) {
	isConvertedEndcap = true;
      }
      //cout<<"isIsolatedEndcap "<<isIsolatedEndcap<<" isIDEndcap "<< isIDEndcap<<" isConvertedEndcap "<<isConvertedEndcap<<endl;
    }
   
  }

  if (isBarrelElectrons) {
    rec.push_back(isIsolated);
    rec.push_back(isIDBarrel);
    rec.push_back(isConvertedBarrel);
    return rec;
  }

 
  
  if (isEndcapElectrons) {
    rec.push_back(isIsolated);
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
