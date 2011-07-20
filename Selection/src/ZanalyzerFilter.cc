// -*- C++ -*-
//
// Package:    Zanalyzer
// Class:      Zanalyzer
// 
/**\class Zanalyzer Zanalyzer.cc Zmonitoring/Zanalyzer/src/Zanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vieri Candelise
//         Created:  Wed May 11 14:53:26 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/ZanalyzerFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"
#include "TMath.h"
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

using namespace std;
using namespace edm;
using namespace reco;


//
// constructors and destructor
//
ZanalyzerFilter::ZanalyzerFilter (const edm::ParameterSet & parameters)
{
  theElectronCollectionLabel =parameters.getParameter < InputTag > ("electronCollection");
}



ZanalyzerFilter::~ZanalyzerFilter ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}




//
// member functions
//

// ------------ method called for each event  ------------
bool
ZanalyzerFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{

  using namespace edm;
  Handle < GsfElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ()){
    cout<<"Sorry, no valid electron collection in this event... skip!"<<endl;
    return false;
  }

  // Find the highest and 2nd highest electron in the barrel/endcap that will be selected with WP80

  bool isBarrelElectrons;
  bool isEndcapElectrons;
  bool isIsolatedBarrel;
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIsolatedEndcap;
  bool isIDEndcap;
  bool isConvertedEndcap;
  int elIsAccepted=0;

  for (reco::GsfElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    if (recoElectron->et () <= 25)  continue;
    // Define Isolation variables
    double IsoTrk = (recoElectron->dr03TkSumPt () / recoElectron->et ());
    double IsoEcal = (recoElectron->dr03EcalRecHitSumEt () / recoElectron->et ());
    double IsoHcal = (recoElectron->dr03HcalTowerSumEt () / recoElectron->et ());
    double HE = (recoElectron->hcalOverEcal ());

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
	  && NumberOfExpectedInnerHits <= 0) {
	isConvertedBarrel = true;
      }

    }

    if (isIsolatedBarrel && isIDBarrel && isConvertedBarrel) elIsAccepted++;

  /***** Endcap WP80 Cuts *****/

    if (fabs (recoElectron->eta ()) >= 1.5660 && fabs (recoElectron->eta ()) <= 2.5000) {

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
	  && NumberOfExpectedInnerHits <= 0) {
	isConvertedEndcap = true;
      }
    }
    if (isIsolatedEndcap && isIDEndcap && isConvertedEndcap) elIsAccepted++;

  }

  int acceptGoodEvent=1; //when set to 1, it required at least 2 good electrons... Set to 0 if you want to include single electron events
  if (elIsAccepted>acceptGoodEvent) {
    return true;
  }
  return false;
  
}


// ------------ method called once each job just before starting event loop  ------------
void
ZanalyzerFilter::beginJob (){
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZanalyzerFilter::endJob ()
{
}

DEFINE_FWK_MODULE (ZanalyzerFilter);
