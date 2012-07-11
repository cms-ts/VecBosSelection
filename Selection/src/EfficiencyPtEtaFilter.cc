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
// Original Author:  superben
//         Created:  Wed May 11 14:53:26 CESDo2011
// $Id: EfficiencyPtEtaFilter.cc,v 1.10 2012/07/11 21:35:15 schizzi Exp $



// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/EfficiencyPtEtaFilter.h"
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

bool Debug_flag=false;

//
// member functions
//

// ------------ method called for each event  ------------
bool
EfficiencyPtEtaFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  if (Debug_flag) cout<<"------- NEW Event -----"<<endl;
  using namespace edm;

  ////////////////////////////////////////
  // Pick up SUPERCLUSTERS / CALOmuons: //
  ////////////////////////////////////////

  Handle<reco::SuperClusterCollection> superClusters_EB_h;
  iEvent.getByLabel(superClusterCollection_EB_,superClusters_EB_h );
  if ( ! superClusters_EB_h.isValid() ) return false;

  Handle<reco::SuperClusterCollection> superClusters_EE_h;
  iEvent.getByLabel(superClusterCollection_EE_,superClusters_EE_h );
  if ( ! superClusters_EE_h.isValid() ) return false;

  reco::SuperClusterCollection::const_iterator tag_SC;
  reco::SuperClusterCollection::const_iterator probe_SC;

  bool RECO_isPassingProbe=false;
  int l=0;

  if (!muonEfficiency_) {
    //EB superclusters:
    for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EB_h->begin(); superCluster != superClusters_EB_h->end(); superCluster++) {
      if ((superCluster->energy()/cosh(superCluster->eta()))>20.0 && fabs(superCluster->eta())<=1.4442) {
	if (l==0) tag_SC=superCluster;
	if (l==1){
	  if (tag_SC->energy()/cosh(tag_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())){
	    probe_SC=tag_SC;
	    tag_SC=superCluster;
	  }
	  else{
	    probe_SC=superCluster;
	  }
	}
	if (l>1){
	  if (tag_SC->energy()/cosh(tag_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())){
	    probe_SC=tag_SC;
	    tag_SC=superCluster;
	  }
	  else{
	    if (probe_SC->energy()/cosh(probe_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())) probe_SC=superCluster;
	  }
	}
	l++;
      }
    }
    //EE superclusters:  
    for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EE_h->begin(); superCluster != superClusters_EE_h->end(); superCluster++) {
      if ((superCluster->energy()/cosh(superCluster->eta()))>20.0 && (fabs(superCluster->eta())>=1.5660 && fabs(superCluster->eta())<=2.4000)) {
	if (l==0) tag_SC=superCluster;
	if (l==1){
	  if (tag_SC->energy()/cosh(tag_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())){
	    probe_SC=tag_SC;
	    tag_SC=superCluster;
	  }
	  else{
	    probe_SC=superCluster;
	  }
	}
	if (l>1){
	  if (tag_SC->energy()/cosh(tag_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())){
	    probe_SC=tag_SC;
	    tag_SC=superCluster;
	  }
	  else{
	    if (probe_SC->energy()/cosh(probe_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())) probe_SC=superCluster;
	  }
	}
	l++;
      }
    }
  } else {
    // CALO muons:
    if (Debug_flag) cout<<"No CALO muons yet!"<<endl;
  }


  if (Debug_flag && l<1) cout<<"No valid SuperCluster or CALOmuon!"<<endl;
  if (!muonEfficiency_ && l<2) return false;
  // Mixing TAG and PROBE SuperClusters
  if (!muonEfficiency_ && RECO_efficiency_) {
    int SCchoice = rand()%2;
    if (SCchoice == 1) {
      reco::SuperClusterCollection::const_iterator temp_SC;
      temp_SC = tag_SC;
      tag_SC = probe_SC;
      probe_SC = temp_SC;
    }
  }


  //////////////////////////////
  // Pick up ELECTRONS/MUONS: //
  //////////////////////////////

  Handle < pat::ElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ()) return false;

  Handle < pat::MuonCollection > muonCollection;
  iEvent.getByLabel (theMuonCollectionLabel, muonCollection);
  if (!muonCollection.isValid ()) return false;

  pat::ElectronCollection::const_iterator tag_ele;
  pat::ElectronCollection::const_iterator probe_ele;
  pat::MuonCollection::const_iterator tag_muon;
  pat::MuonCollection::const_iterator probe_muon;

  pat::MuonCollection::const_iterator probe_CALOmu; ///TEMPORARY!!!!!!!!!

  // iso deposits
  IsoDepositVals isoVals(isoValInputTags_.size());
  for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
     iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
  }

  int i=0;
  bool protection=false;
  double deltaPhi = 0.0;

  // Electron LOOP:
  if (!muonEfficiency_) {
    if (electronCollection->size()<=1) return false;
    for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
      protection=true;
      if (recoElectron->pt()>20.0 && (fabs(recoElectron->superCluster()->eta())<1.4442 || (fabs(recoElectron->superCluster()->eta())>1.566 && fabs(recoElectron->superCluster()->eta())<2.4)) && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_)) {
	if (i==0) tag_ele=recoElectron;
	if (i==1){
	  if (tag_ele->pt()<recoElectron->pt()){
	    probe_ele=tag_ele;
	    tag_ele=recoElectron;
	  }
	  else{
	    probe_ele=recoElectron;
	  }
	}
	if (i>1){
	  if (tag_ele->pt()<recoElectron->pt()){
	    probe_ele=tag_ele;
	    tag_ele=recoElectron;
	  }
	  else{
	    if (probe_ele->pt()<recoElectron->pt()) probe_ele=recoElectron;
	  }
	}
	i++;
      }
      if (recoElectron->pt()>20.0 && (fabs(recoElectron->superCluster()->eta())<1.4442 || (fabs(recoElectron->superCluster()->eta())>1.566 && fabs(recoElectron->superCluster()->eta())<2.4)) && RECO_efficiency_) {
	if (fabs(probe_SC->phi()-recoElectron->superCluster()->phi()) < 3.1416) {
	  deltaPhi = probe_SC->phi()-recoElectron->superCluster()->phi();
	} else {
	  deltaPhi = fabs(probe_SC->phi()-recoElectron->superCluster()->phi()) - 6.2832;
	}
	if (sqrt((probe_SC->eta()-recoElectron->superCluster()->eta())*(probe_SC->eta()-recoElectron->superCluster()->eta())+(deltaPhi*deltaPhi)) < 0.2) {
	  RECO_isPassingProbe=true;
	  continue;
	}
	if (i==0) tag_ele=recoElectron;
	if (i>0){
	  if (tag_ele->pt()<recoElectron->pt()){
	    tag_ele=recoElectron;
	  }
	}
	i++;
      }
    }
  }

  // Muon LOOP:
  if (muonEfficiency_) {
    if (muonCollection->size()<=1) return false;
    for (pat::MuonCollection::const_iterator recoMuon = muonCollection->begin (); recoMuon != muonCollection->end (); recoMuon++) {
      protection=true;
      if (recoMuon->pt()>20.0 && (fabs(recoMuon->eta())<1.4442 || (fabs(recoMuon->eta())>1.566 && fabs(recoMuon->eta())<2.4)) && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_)) {
	if (i==0) tag_muon=recoMuon;
	if (i==1){
	  if (tag_muon->pt()<recoMuon->pt()){
	    probe_muon=tag_muon;
	    tag_muon=recoMuon;
	  }
	  else{
	    probe_muon=recoMuon;
	  }
	}
	if (i>1){
	  if (tag_muon->pt()<recoMuon->pt()){
	    probe_muon=tag_muon;
	    tag_muon=recoMuon;
	  }
	  else{
	    if (probe_muon->pt()<recoMuon->pt()) probe_muon=recoMuon;
	  }
	}
	i++;
      }
      if (recoMuon->pt()>20.0 && (fabs(recoMuon->eta())<1.4442 || (fabs(recoMuon->eta())>1.566 && fabs(recoMuon->eta())<2.4)) && RECO_efficiency_) {
	if (fabs(probe_SC->phi()-recoMuon->phi()) < 3.1416) {
	  deltaPhi = probe_SC->phi()-recoMuon->phi();
	} else {
	  deltaPhi = fabs(probe_SC->phi()-recoMuon->phi()) - 6.2832;
	}
	if (sqrt((probe_SC->eta()-recoMuon->eta())*(probe_SC->eta()-recoMuon->eta())+(deltaPhi*deltaPhi)) < 0.2) {
	  RECO_isPassingProbe=true;
	  continue;
	}
	if (i==0) tag_muon=recoMuon;
	if (i>0){
	  if (tag_muon->pt()<recoMuon->pt()){
	    tag_muon=recoMuon;
	  }
	}
	i++;
      }
    }
  }

  if (!protection) {
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }
  
  //  if (Debug_flag) cout<<"tag_ele->pt() = "<<tag_ele->pt() <<endl;  
  if (Debug_flag) cout<<"Out of the SC and RECOele/mu loops!"<<endl;  

  if (i<1) return false;
  if ((!muonEfficiency_) && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) && (i<2 || tag_ele->charge() == probe_ele->charge())) return false;
  if ((muonEfficiency_)  && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) && (i<2 || tag_muon->charge() == probe_muon->charge())) return false;

  // Mixing TAG and PROBE electrons:

  if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
    int tagchoice = rand()%2;
    if (tagchoice == 1) {
      if (!muonEfficiency_) {
	pat::ElectronCollection::const_iterator temp_ele;
	temp_ele = tag_ele;
	tag_ele = probe_ele;
	probe_ele = temp_ele;
      }
      if (muonEfficiency_) {
	pat::MuonCollection::const_iterator temp_muon;
        temp_muon = tag_muon;
        tag_muon = probe_muon;
        probe_muon = temp_muon;
      }
    }
  }

  ////////////////////////////////////////////////////
  probe_CALOmu = tag_muon; // TEMPORARY!!!!!!!!!!!! //
  ////////////////////////////////////////////////////

  if (Debug_flag) cout << "Finished mixing Tag and Probe eles." << endl;


  //////////////
  // MC match //
  //////////////

  if (matchMC_) {

    Handle<reco::GenParticleCollection> genPart;
    iEvent.getByLabel(genParticleCollection_, genPart);

    if (genPart.isValid()) {

      int nMatch = 0;

      for  (reco::GenParticleCollection::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {

	// electron MC truth:

        if (!muonEfficiency_ && itgen->status() == 3 && fabs(itgen->pdgId()) == 11 && itgen->mother()->pdgId() == 23) {

          double deltaPhi1, deltaR1;
          if (fabs(itgen->phi()-tag_ele->phi()) < 3.1416) {
            deltaPhi1 = itgen->phi()-tag_ele->phi();
          } else {
            deltaPhi1 = fabs(itgen->phi()-tag_ele->phi()) - 6.2832;
          }
          deltaR1 = sqrt( pow(itgen->eta()-tag_ele->eta(), 2) + deltaPhi1*deltaPhi1 );

          double deltaPhi2, deltaR2;
          if ( RECO_efficiency_) {
            if (fabs(itgen->phi()-probe_SC->phi()) < 3.1416) {
              deltaPhi2 = itgen->phi()-probe_SC->phi();
            } else {
              deltaPhi2 = fabs(itgen->phi()-probe_SC->phi()) - 6.2832;
            }
            deltaR2 = sqrt( pow(itgen->eta()-probe_SC->eta(), 2) + deltaPhi2*deltaPhi2 );
          } else {
            if (fabs(itgen->phi()-probe_ele->phi()) < 3.1416) {
              deltaPhi2 = itgen->phi()-probe_ele->phi();
            } else {
              deltaPhi2 = fabs(itgen->phi()-probe_ele->phi()) - 6.2832;
            }
            deltaR2 = sqrt( pow(itgen->eta()-probe_ele->eta(), 2) + deltaPhi2*deltaPhi2 );
          }

          if ( deltaR1 < 0.1 && tag_ele->charge() == itgen->charge() ) {
            nMatch++;
            if (Debug_flag) std::cout << "match 1" << std::endl;
          } else if ( deltaR2 < 0.1 && RECO_efficiency_ ) {
            if (Debug_flag) std::cout << "match 2" << std::endl;
            nMatch++;
          } else if ( deltaR2 < 0.1 && !RECO_efficiency_ && probe_ele->charge() == itgen->charge() ) {
            if (Debug_flag) std::cout << "match 2" << std::endl;
            nMatch++;
          } else {
            if (Debug_flag) std::cout << "no match" << std::endl;
          }

        }

	// muon MC truth:

        if (muonEfficiency_ && itgen->status() == 3 && fabs(itgen->pdgId()) == 13 && itgen->mother()->pdgId() == 23) {

          double deltaPhi1, deltaR1;
          if (fabs(itgen->phi()-tag_muon->phi()) < 3.1416) {
            deltaPhi1 = itgen->phi()-tag_muon->phi();
          } else {
            deltaPhi1 = fabs(itgen->phi()-tag_muon->phi()) - 6.2832;
          }
          deltaR1 = sqrt( pow(itgen->eta()-tag_muon->eta(), 2) + deltaPhi1*deltaPhi1 );

          double deltaPhi2, deltaR2;
          if ( RECO_efficiency_) {
            if (fabs(itgen->phi()-probe_CALOmu->phi()) < 3.1416) {
              deltaPhi2 = itgen->phi()-probe_CALOmu->phi();
            } else {
              deltaPhi2 = fabs(itgen->phi()-probe_CALOmu->phi()) - 6.2832;
            }
            deltaR2 = sqrt( pow(itgen->eta()-probe_CALOmu->eta(), 2) + deltaPhi2*deltaPhi2 );
          } else {
            if (fabs(itgen->phi()-probe_muon->phi()) < 3.1416) {
              deltaPhi2 = itgen->phi()-probe_muon->phi();
            } else {
              deltaPhi2 = fabs(itgen->phi()-probe_muon->phi()) - 6.2832;
            }
            deltaR2 = sqrt( pow(itgen->eta()-probe_muon->eta(), 2) + deltaPhi2*deltaPhi2 );
          }

          if ( deltaR1 < 0.1 && tag_muon->charge() == itgen->charge() ) {
            nMatch++;
            if (Debug_flag) std::cout << "match 1" << std::endl;
          } else if ( deltaR2 < 0.1 && RECO_efficiency_ ) {
            if (Debug_flag) std::cout << "match 2" << std::endl;
            nMatch++;
          } else if ( deltaR2 < 0.1 && !RECO_efficiency_ && probe_muon->charge() == itgen->charge() ) {
            if (Debug_flag) std::cout << "match 2" << std::endl;
            nMatch++;
          } else {
            if (Debug_flag) std::cout << "no match" << std::endl;
          }

        }

      }

      if (Debug_flag) std::cout << "# match = " << nMatch << std::endl;

      if (nMatch > 2) std::cout << "TOO MANY MATCHED ELECTRONS" << std::endl;

      if (nMatch != 2) return false;

    }

  }


  //////////////////////////////////////////////
  // TRIGGER and ID/ISO(muons-only) MATCHING: //
  //////////////////////////////////////////////

  // TAG trigger matching:

  Handle < edm::RefToBaseVector<reco::GsfElectron> > TagHLTelectronCollection;
  iEvent.getByLabel (theTagHLTElectronCollectionLabel, TagHLTelectronCollection);
  if (!TagHLTelectronCollection.isValid ()) return false;

  Handle < edm::RefToBaseVector<reco::Muon> > TagHLTmuonCollection;
  iEvent.getByLabel (theTagHLTMuonCollectionLabel, TagHLTmuonCollection);
  if (!TagHLTmuonCollection.isValid ()) return false;

  Handle < pat::MuonCollection > tightMuonCollection;
  iEvent.getByLabel (theTightMuonCollectionLabel, tightMuonCollection);
  if (!tightMuonCollection.isValid ()) return false;


  bool HLTmatch = false;
  bool IDISOmatch = false;

  double deltaReles = 999.0;

  if (!muonEfficiency_) {
    for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator TagHLTElectron = TagHLTelectronCollection->begin (); TagHLTElectron != TagHLTelectronCollection->end (); TagHLTElectron++) {
      if ((*TagHLTElectron)->pt()>20.0) {
	if (fabs(tag_ele->phi()-(*TagHLTElectron)->phi()) < 3.1416) {
	  deltaPhi = tag_ele->phi()-(*TagHLTElectron)->phi();
	} else {
	  deltaPhi = fabs(tag_ele->phi()-(*TagHLTElectron)->phi()) - 6.2832;
	}
	deltaReles = sqrt((tag_ele->eta()-(*TagHLTElectron)->eta())*
			  (tag_ele->eta()-(*TagHLTElectron)->eta())+
			  (deltaPhi*deltaPhi));
	if (deltaReles < 0.2) HLTmatch = true;
      }
    }
    if (!(HLTmatch && SelectionUtils::DoWP80Pf(tag_ele,iEvent) && SelectionUtils::DoIso2011(tag_ele, iEvent, isoVals))) return false;
  } else {
    // MuonTag HLT matching:
    for (edm::RefToBaseVector<reco::Muon>::const_iterator TagHLTMuon = TagHLTmuonCollection->begin (); TagHLTMuon != TagHLTmuonCollection->end (); TagHLTMuon++) {
      if ((*TagHLTMuon)->pt()>20.0) {
	if (fabs(tag_muon->phi()-(*TagHLTMuon)->phi()) < 3.1416) {
	  deltaPhi = tag_muon->phi()-(*TagHLTMuon)->phi();
	} else {
	  deltaPhi = fabs(tag_muon->phi()-(*TagHLTMuon)->phi()) - 6.2832;
	}
	deltaReles = sqrt((tag_muon->eta()-(*TagHLTMuon)->eta())*
			  (tag_muon->eta()-(*TagHLTMuon)->eta())+
			  (deltaPhi*deltaPhi));
	if (deltaReles < 0.2) HLTmatch = true;
      }
    }
    // MuonTag ID/ISO
    for (pat::MuonCollection::const_iterator tightMuon = tightMuonCollection->begin (); tightMuon != tightMuonCollection->end (); tightMuon++) {
      if (tightMuon->pt()>20.0) {
	if (fabs(tag_muon->phi()-tightMuon->phi()) < 3.1416) {
	  deltaPhi = tag_muon->phi()-tightMuon->phi();
	} else {
	  deltaPhi = fabs(tag_muon->phi()-tightMuon->phi()) - 6.2832;
	}
	deltaReles = sqrt((tag_muon->eta()-tightMuon->eta())*
			  (tag_muon->eta()-tightMuon->eta())+
			  (deltaPhi*deltaPhi));
	if (deltaReles < 0.2) IDISOmatch = true;
      }
    }
    if (!(HLTmatch && IDISOmatch)) return false; // stop if tag is not HLT matched and ID/ISO!!!
  }

  // PROBE trigger matching and ID/ISO:

  Handle < edm::RefToBaseVector<reco::GsfElectron> > ProbeHLTelectronCollection;
  iEvent.getByLabel (theProbeHLTElectronCollectionLabel, ProbeHLTelectronCollection);
  if (!ProbeHLTelectronCollection.isValid ()) return false;

  Handle < edm::RefToBaseVector<reco::Muon> > ProbeHLTmuonCollection;
  iEvent.getByLabel (theProbeHLTMuonCollectionLabel, ProbeHLTmuonCollection);
  if (!ProbeHLTmuonCollection.isValid ()) return false;
  
  if (HLTele17_efficiency_ || HLTele8_efficiency_) {
    HLTmatch = false;
    if (!muonEfficiency_) { //ElectronProbe trg matching
      for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator ProbeHLTElectron = ProbeHLTelectronCollection->begin (); ProbeHLTElectron != ProbeHLTelectronCollection->end (); ProbeHLTElectron++) {
	if ((*ProbeHLTElectron)->pt()>20.0) {
	  if (fabs(probe_ele->phi()-(*ProbeHLTElectron)->phi()) < 3.1416) {
	    deltaPhi = probe_ele->phi()-(*ProbeHLTElectron)->phi();
	  } else {
	    deltaPhi = fabs(probe_ele->phi()-(*ProbeHLTElectron)->phi()) - 6.2832;
	  }
	  deltaReles = sqrt((probe_ele->eta()-(*ProbeHLTElectron)->eta())*
			    (probe_ele->eta()-(*ProbeHLTElectron)->eta())+
			    (deltaPhi*deltaPhi));
	  if (deltaReles < 0.2) HLTmatch = true;
	}
      }
    } else { //MuonProbe trg matching
      for (edm::RefToBaseVector<reco::Muon>::const_iterator ProbeHLTMuon = ProbeHLTmuonCollection->begin (); ProbeHLTMuon != ProbeHLTmuonCollection->end (); ProbeHLTMuon++) {
	if ((*ProbeHLTMuon)->pt()>20.0) {
	  if (fabs(probe_muon->phi()-(*ProbeHLTMuon)->phi()) < 3.1416) {
	    deltaPhi = probe_muon->phi()-(*ProbeHLTMuon)->phi();
	  } else {
	    deltaPhi = fabs(probe_muon->phi()-(*ProbeHLTMuon)->phi()) - 6.2832;
	  }
	  deltaReles = sqrt((probe_muon->eta()-(*ProbeHLTMuon)->eta())*
			    (probe_muon->eta()-(*ProbeHLTMuon)->eta())+
			    (deltaPhi*deltaPhi));
	  if (deltaReles < 0.2) HLTmatch = true;
	}
      }
    }
  }

  // MuonProbe ID/ISO:
  if (WP80_efficiency_ && muonEfficiency_) {
    IDISOmatch = false;
    for (pat::MuonCollection::const_iterator tightMuon = tightMuonCollection->begin (); tightMuon != tightMuonCollection->end (); tightMuon++) {
      if (tightMuon->pt()>20.0) {
	if (fabs(probe_muon->phi()-tightMuon->phi()) < 3.1416) {
	  deltaPhi = probe_muon->phi()-tightMuon->phi();
	} else {
	  deltaPhi = fabs(probe_muon->phi()-tightMuon->phi()) - 6.2832;
	}
	deltaReles = sqrt((probe_muon->eta()-tightMuon->eta())*
			  (probe_muon->eta()-tightMuon->eta())+
			  (deltaPhi*deltaPhi));
	if (deltaReles < 0.2) IDISOmatch = true;
      }
    }
  }

  if (Debug_flag) cout<<"Trigger matching finished: HLTmatch = "<<HLTmatch<<endl;  


  ///////////////////////////////////
  // Calculating tag & probe stuff //
  ///////////////////////////////////

  TLorentzVector tagv;
  TLorentzVector probev;

  if (!muonEfficiency_) {
    tagv.SetPtEtaPhiM(tag_ele->pt(),tag_ele->eta(),tag_ele->phi(), 0.0);
    if (RECO_efficiency_)  {
      probev.SetXYZT(probe_SC->energy()*cos(probe_SC->phi())/cosh(probe_SC->eta()),
		     probe_SC->energy()*sin(probe_SC->phi())/cosh(probe_SC->eta()),
		     probe_SC->energy()*tanh(probe_SC->eta()),
		     probe_SC->energy());
    } else {
      probev.SetPtEtaPhiM(probe_ele->pt(),probe_ele->eta(),probe_ele->phi(), 0.0);
    }
  } else {
    tagv.SetPtEtaPhiM(tag_muon->pt(),tag_muon->eta(),tag_muon->phi(), 0.0);
    if (RECO_efficiency_)  {
    tagv.SetPtEtaPhiM(probe_CALOmu->pt(),probe_CALOmu->eta(),probe_CALOmu->phi(), 0.0);
    } else {
      probev.SetPtEtaPhiM(probe_muon->pt(),probe_muon->eta(),probe_muon->phi(), 0.0);
    }
  }


  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();

  //cut on the tag and probe mass...
  if (e_ee_invMass>120.0 || e_ee_invMass<60.0) return false;

  if (Debug_flag) cout<<"Start filling TAP histos..."<<endl;
  

  ////////////////////////////////
  // Filling TAP distributions: //
  ////////////////////////////////

  // WP80 & HLT:
  if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
    if (!muonEfficiency_) {  // Start filling electron histograms
      probeall_pt->Fill(probe_ele->pt());
      probeall_eta->Fill(probe_ele->superCluster()->eta());
      probeall_mee->Fill(e_ee_invMass);
      tagall_pt->Fill(tag_ele->pt());
      tagall_eta->Fill(tag_ele->superCluster()->eta());
      if ((WP80_efficiency_ && SelectionUtils::DoWP90Pf(probe_ele,iEvent) && SelectionUtils::DoIso2011(probe_ele, iEvent, isoVals)) || ((HLTele17_efficiency_ || HLTele8_efficiency_) && HLTmatch)) {
	probepass_pt->Fill(probe_ele->pt());
	probepass_eta->Fill(probe_ele->superCluster()->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (fabs(probe_ele->superCluster()->eta()) >= 0.0    && fabs(probe_ele->superCluster()->eta()) < 0.8) {
	  if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass1eta2pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass1eta3pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass1eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_ele->superCluster()->eta()) >= 0.8    && fabs(probe_ele->superCluster()->eta()) < 1.4442) {
	  if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass2eta1pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass2eta3pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass2eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_ele->superCluster()->eta()) >= 1.566  && fabs(probe_ele->superCluster()->eta()) < 2.0) {
	  if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass3eta1pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass3eta2pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass3eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_ele->superCluster()->eta()) >= 2.0    && fabs(probe_ele->superCluster()->eta()) < 2.5) {
	  if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass4eta1pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass4eta2pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass4eta3pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
	}
      } else {
	probefail_pt->Fill(probe_ele->pt());
	probefail_eta->Fill(probe_ele->superCluster()->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (fabs(probe_ele->superCluster()->eta()) >= 0.0    && fabs(probe_ele->superCluster()->eta()) < 0.8) {
	  if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail1eta1pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail1eta2pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail1eta3pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail1eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_ele->superCluster()->eta()) >= 0.8    && fabs(probe_ele->superCluster()->eta()) < 1.4442) {
	  if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail2eta1pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail2eta3pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail2eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_ele->superCluster()->eta()) >= 1.566  && fabs(probe_ele->superCluster()->eta()) < 2.0) {
	  if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail3eta1pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail3eta2pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail3eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_ele->superCluster()->eta()) >= 2.0    && fabs(probe_ele->superCluster()->eta()) < 2.5) {
	  if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail4eta1pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail4eta2pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail4eta3pt->Fill(e_ee_invMass);
	  if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
	}
      }
    } else { // Start filling muon histograms
      probeall_pt->Fill(probe_muon->pt());
      probeall_eta->Fill(probe_muon->eta());
      probeall_mee->Fill(e_ee_invMass);
      tagall_pt->Fill(tag_muon->pt());
      tagall_eta->Fill(tag_muon->eta());
      if ((WP80_efficiency_ && IDISOmatch) || ((HLTele17_efficiency_ || HLTele8_efficiency_) && HLTmatch)) {
	probepass_pt->Fill(probe_muon->pt());
	probepass_eta->Fill(probe_muon->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (fabs(probe_muon->eta()) >= 0.0    && fabs(probe_muon->eta()) < 0.8) {
	  if (probe_muon->pt() >= 20.0 && probe_muon->pt() <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 30.0 && probe_muon->pt() <  40.0) probepass1eta2pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 40.0 && probe_muon->pt() <  50.0) probepass1eta3pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 50.0 && probe_muon->pt() < 999.0) probepass1eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_muon->eta()) >= 0.8    && fabs(probe_muon->eta()) < 1.4442) {
	  if (probe_muon->pt() >= 20.0 && probe_muon->pt() <  30.0) probepass2eta1pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 30.0 && probe_muon->pt() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 40.0 && probe_muon->pt() <  50.0) probepass2eta3pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 50.0 && probe_muon->pt() < 999.0) probepass2eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_muon->eta()) >= 1.566  && fabs(probe_muon->eta()) < 2.0) {
	  if (probe_muon->pt() >= 20.0 && probe_muon->pt() <  30.0) probepass3eta1pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 30.0 && probe_muon->pt() <  40.0) probepass3eta2pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 40.0 && probe_muon->pt() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 50.0 && probe_muon->pt() < 999.0) probepass3eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_muon->eta()) >= 2.0    && fabs(probe_muon->eta()) < 2.5) {
	  if (probe_muon->pt() >= 20.0 && probe_muon->pt() <  30.0) probepass4eta1pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 30.0 && probe_muon->pt() <  40.0) probepass4eta2pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 40.0 && probe_muon->pt() <  50.0) probepass4eta3pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 50.0 && probe_muon->pt() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
	}
      } else {
	probefail_pt->Fill(probe_muon->pt());
	probefail_eta->Fill(probe_muon->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (fabs(probe_muon->eta()) >= 0.0    && fabs(probe_muon->eta()) < 0.8) {
	  if (probe_muon->pt() >= 20.0 && probe_muon->pt() <  30.0) probefail1eta1pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 30.0 && probe_muon->pt() <  40.0) probefail1eta2pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 40.0 && probe_muon->pt() <  50.0) probefail1eta3pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 50.0 && probe_muon->pt() < 999.0) probefail1eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_muon->eta()) >= 0.8    && fabs(probe_muon->eta()) < 1.4442) {
	  if (probe_muon->pt() >= 20.0 && probe_muon->pt() <  30.0) probefail2eta1pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 30.0 && probe_muon->pt() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 40.0 && probe_muon->pt() <  50.0) probefail2eta3pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 50.0 && probe_muon->pt() < 999.0) probefail2eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_muon->eta()) >= 1.566  && fabs(probe_muon->eta()) < 2.0) {
	  if (probe_muon->pt() >= 20.0 && probe_muon->pt() <  30.0) probefail3eta1pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 30.0 && probe_muon->pt() <  40.0) probefail3eta2pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 40.0 && probe_muon->pt() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 50.0 && probe_muon->pt() < 999.0) probefail3eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_muon->eta()) >= 2.0    && fabs(probe_muon->eta()) < 2.5) {
	  if (probe_muon->pt() >= 20.0 && probe_muon->pt() <  30.0) probefail4eta1pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 30.0 && probe_muon->pt() <  40.0) probefail4eta2pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 40.0 && probe_muon->pt() <  50.0) probefail4eta3pt->Fill(e_ee_invMass);
	  if (probe_muon->pt() >= 50.0 && probe_muon->pt() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
	}
      }
    }
  }

  // RECO:
  if (RECO_efficiency_) {
    if (!muonEfficiency_) { //Start filling electron hisyograms:
      probeall_pt->Fill((probe_SC->energy()/cosh(probe_SC->eta())));
      probeall_eta->Fill(probe_SC->eta());
      probeall_mee->Fill(e_ee_invMass);
      tagall_pt->Fill(tag_ele->pt());
      tagall_eta->Fill(tag_ele->superCluster()->eta());
      if (RECO_isPassingProbe){
	probepass_pt->Fill((probe_SC->energy()/cosh(probe_SC->eta())));
	probepass_eta->Fill(probe_SC->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (fabs(probe_SC->eta()) >= 0.0    && fabs(probe_SC->eta()) < 0.8) {
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 20.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 30.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  40.0) probepass1eta2pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 40.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  50.0) probepass1eta3pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 50.0 && (probe_SC->energy()/cosh(probe_SC->eta())) < 999.0) probepass1eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_SC->eta()) >= 0.8    && fabs(probe_SC->eta()) < 1.4442) {
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 20.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  30.0) probepass2eta1pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 30.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 40.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  50.0) probepass2eta3pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 50.0 && (probe_SC->energy()/cosh(probe_SC->eta())) < 999.0) probepass2eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_SC->eta()) >= 1.566  && fabs(probe_SC->eta()) < 2.0) {
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 20.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  30.0) probepass3eta1pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 30.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  40.0) probepass3eta2pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 40.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 50.0 && (probe_SC->energy()/cosh(probe_SC->eta())) < 999.0) probepass3eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_SC->eta()) >= 2.0    && fabs(probe_SC->eta()) < 2.5) {
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 20.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  30.0) probepass4eta1pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 30.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  40.0) probepass4eta2pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 40.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  50.0) probepass4eta3pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 50.0 && (probe_SC->energy()/cosh(probe_SC->eta())) < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
	}
      } else {
	probefail_pt->Fill((probe_SC->energy()/cosh(probe_SC->eta())));
	probefail_eta->Fill(probe_SC->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (fabs(probe_SC->eta()) >= 0.0    && fabs(probe_SC->eta()) < 0.8) {
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 20.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  30.0) probefail1eta1pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 30.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  40.0) probefail1eta2pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 40.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  50.0) probefail1eta3pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 50.0 && (probe_SC->energy()/cosh(probe_SC->eta())) < 999.0) probefail1eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_SC->eta()) >= 0.8    && fabs(probe_SC->eta()) < 1.4442) {
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 20.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  30.0) probefail2eta1pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 30.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 40.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  50.0) probefail2eta3pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 50.0 && (probe_SC->energy()/cosh(probe_SC->eta())) < 999.0) probefail2eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_SC->eta()) >= 1.566  && fabs(probe_SC->eta()) < 2.0) {
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 20.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  30.0) probefail3eta1pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 30.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  40.0) probefail3eta2pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 40.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 50.0 && (probe_SC->energy()/cosh(probe_SC->eta())) < 999.0) probefail3eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_SC->eta()) >= 2.0    && fabs(probe_SC->eta()) < 2.5) {
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 20.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  30.0) probefail4eta1pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 30.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  40.0) probefail4eta2pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 40.0 && (probe_SC->energy()/cosh(probe_SC->eta())) <  50.0) probefail4eta3pt->Fill(e_ee_invMass);
	  if ((probe_SC->energy()/cosh(probe_SC->eta())) >= 50.0 && (probe_SC->energy()/cosh(probe_SC->eta())) < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
	}
      } 
    } else { // Start filling muon histograms:
      probeall_pt->Fill(probe_CALOmu->pt());
      probeall_eta->Fill(probe_CALOmu->eta());
      probeall_mee->Fill(e_ee_invMass);
      tagall_pt->Fill(tag_muon->pt());
      tagall_eta->Fill(tag_muon->superCluster()->eta());
      if (RECO_isPassingProbe){
	probepass_pt->Fill(probe_CALOmu->pt());
	probepass_eta->Fill(probe_CALOmu->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (fabs(probe_CALOmu->eta()) >= 0.0    && fabs(probe_CALOmu->eta()) < 0.8) {
	  if (probe_CALOmu->pt() >= 20.0 && probe_CALOmu->pt() <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 30.0 && probe_CALOmu->pt() <  40.0) probepass1eta2pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 40.0 && probe_CALOmu->pt() <  50.0) probepass1eta3pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 50.0 && probe_CALOmu->pt() < 999.0) probepass1eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_CALOmu->eta()) >= 0.8    && fabs(probe_CALOmu->eta()) < 1.4442) {
	  if (probe_CALOmu->pt() >= 20.0 && probe_CALOmu->pt() <  30.0) probepass2eta1pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 30.0 && probe_CALOmu->pt() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 40.0 && probe_CALOmu->pt() <  50.0) probepass2eta3pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 50.0 && probe_CALOmu->pt() < 999.0) probepass2eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_CALOmu->eta()) >= 1.566  && fabs(probe_CALOmu->eta()) < 2.0) {
	  if (probe_CALOmu->pt() >= 20.0 && probe_CALOmu->pt() <  30.0) probepass3eta1pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 30.0 && probe_CALOmu->pt() <  40.0) probepass3eta2pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 40.0 && probe_CALOmu->pt() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 50.0 && probe_CALOmu->pt() < 999.0) probepass3eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_CALOmu->eta()) >= 2.0    && fabs(probe_CALOmu->eta()) < 2.5) {
	  if (probe_CALOmu->pt() >= 20.0 && probe_CALOmu->pt() <  30.0) probepass4eta1pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 30.0 && probe_CALOmu->pt() <  40.0) probepass4eta2pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 40.0 && probe_CALOmu->pt() <  50.0) probepass4eta3pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 50.0 && probe_CALOmu->pt() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
	}
      } else {
	probefail_pt->Fill(probe_CALOmu->pt());
	probefail_eta->Fill(probe_CALOmu->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (fabs(probe_CALOmu->eta()) >= 0.0    && fabs(probe_CALOmu->eta()) < 0.8) {
	  if (probe_CALOmu->pt() >= 20.0 && probe_CALOmu->pt() <  30.0) probefail1eta1pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 30.0 && probe_CALOmu->pt() <  40.0) probefail1eta2pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 40.0 && probe_CALOmu->pt() <  50.0) probefail1eta3pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 50.0 && probe_CALOmu->pt() < 999.0) probefail1eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_CALOmu->eta()) >= 0.8    && fabs(probe_CALOmu->eta()) < 1.4442) {
	  if (probe_CALOmu->pt() >= 20.0 && probe_CALOmu->pt() <  30.0) probefail2eta1pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 30.0 && probe_CALOmu->pt() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 40.0 && probe_CALOmu->pt() <  50.0) probefail2eta3pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 50.0 && probe_CALOmu->pt() < 999.0) probefail2eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_CALOmu->eta()) >= 1.566  && fabs(probe_CALOmu->eta()) < 2.0) {
	  if (probe_CALOmu->pt() >= 20.0 && probe_CALOmu->pt() <  30.0) probefail3eta1pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 30.0 && probe_CALOmu->pt() <  40.0) probefail3eta2pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 40.0 && probe_CALOmu->pt() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 50.0 && probe_CALOmu->pt() < 999.0) probefail3eta4pt->Fill(e_ee_invMass);
	}
	if (fabs(probe_CALOmu->eta()) >= 2.0    && fabs(probe_CALOmu->eta()) < 2.5) {
	  if (probe_CALOmu->pt() >= 20.0 && probe_CALOmu->pt() <  30.0) probefail4eta1pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 30.0 && probe_CALOmu->pt() <  40.0) probefail4eta2pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 40.0 && probe_CALOmu->pt() <  50.0) probefail4eta3pt->Fill(e_ee_invMass);
	  if (probe_CALOmu->pt() >= 50.0 && probe_CALOmu->pt() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
	}
      } 
    }
  }

  if (Debug_flag) cout<<"Finished event analysis!"<<endl;  
  return true;
}

  
// ------------ method called once each job just before starting event loop  ------------
void
EfficiencyPtEtaFilter::beginJob (){

}

// ------------ method called when starting to processes a run  ------------
bool
EfficiencyPtEtaFilter::beginRun(edm::Run &iRun, edm::EventSetup const& iSetup)
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
   if (Debug_flag) cout<<"useAllTriggers?"<<useAllTriggers_<<endl;
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

   if (Debug_flag){
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
EfficiencyPtEtaFilter::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (EfficiencyPtEtaFilter);
